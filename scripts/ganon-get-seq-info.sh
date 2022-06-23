#!/bin/bash
# Number of attempts to request data from e-utils
att=3
batch=200

retrieve_summary_xml()
{
    echo "$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id=${1}${2}")"
}

retrieve_nucleotide_fasta_xml()
{
    echo "$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=xml&id=${1}${2}")"
}

retrieve_assembly_uid_xml()
{
    echo "$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=assembly&linkname=nuccore_assembly&id=${1}${2}")"
}

retrieve_assembly_accession_xml()
{
    echo "$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=${1}${2}")"
}

get_lines()
{
    echo "$(sed -n "${2},$((${2}+${3}-1))p" ${1})"
}

sort_in_out(){ # $1 file, $2 order
    # Sort output in the same input order
    awk 'BEGIN{FS="\t"} NR==FNR{a[$1]=$0; next} $1 in a{print a[$1]}' <(echo "${1}") <(echo "${2}")
}

function showhelp {
    echo "ganon-get-seq-info.sh"
    echo
    echo "Uses NCBI E-utils http requests to get sequence information (lenght, taxid, assembly accession, assembly name) based on accessions "
    echo "outputs to STDOUT in the same order of the input (besides failed entries, which can be reported with -k)"
    echo
    echo $' -i [str] input file with one accessions per line (use - to read from STDIN)' 
    echo $' -l [str] list of accesions (comma separated)'
    echo $' -n [str] ncbi_api_key'
    echo
    echo $' -k Keep all entries even if nothing is retrieved (report "na")'
    echo $' -r Use sequence accession for unavailable asssembly accessions/names (by default report "na")'
    echo
    echo $' -e Get sequence length and taxid'
    echo $' -a Get assembly accession (latest)'
    echo $' -m Get assembly name'



    echo
}

input_file=""
list_acc=""
ncbi_api_key=""
keep_all=0
get_assembly_accession=0
get_assembly_name=0
get_length_taxid=0
replace_not_found_accession=0

OPTIND=1 # Reset getopts
while getopts "i:l:n:kretam" opt; do
  case ${opt} in
    i) input_file=${OPTARG} ;;
    l) list_acc=${OPTARG} ;;
    n) ncbi_api_key=${OPTARG} ;;
    k) keep_all=1 ;;
    r) replace_not_found_accession=1 ;;
    e) get_length_taxid=1 ;;
    a) get_assembly_accession=1 ;;
    m) get_assembly_name=1 ;;

    h|\?) showhelp; exit 1 ;;
    :) echo "Option -${OPTARG} requires an argument." >&2; exit 1 ;;
  esac
done
if [ ${OPTIND} -eq 1 ]; then showhelp; exit 1; fi
shift $((OPTIND-1))
[ "$1" = "--" ] && shift

if [[ "${get_length_taxid}" -eq 0 && "${get_assembly_accession}" -eq 0 && "${get_assembly_name}" -eq 0 ]]; then
    echo "At least one option has to be active (-e -a -m)"; exit 1;
fi

if [[ "${input_file}" == "-" ]]; then
    tmp_file="$(mktemp)"
    cat /dev/stdin > ${tmp_file}
    input_file=${tmp_file}
elif [[ ! -z "${list_acc}" ]]; then
    tmp_file="$(mktemp)"
    echo -e "${list_acc//,/$'\n'}" > ${tmp_file}
    input_file=${tmp_file}
fi

# test if ncbi api key is valid
if [[ ! -z "${ncbi_api_key}" ]]; then
    ncbi_api_key="&api_key=${ncbi_api_key}"
    api_key_test="$(retrieve_summary_xml "" "${ncbi_api_key}" | grep -o 'API key invalid')"
    if [[ ! -z "${api_key_test}" ]]; then
        ncbi_api_key=""
        (>&2 printf "Invalid NCBI API key\n")
    fi
fi

batch_count=1
acc="$(get_lines ${input_file} 1 ${batch})"
while [[ ! -z "${acc}" ]];
do
    out=""

    if [ "${get_length_taxid}" -eq 1 ]; then
        for i in $(seq 1 ${att});
        do
            # First try to get from summary, lighter resource 
            # replacing \n for , 
            xml_summary="$(retrieve_summary_xml "${acc//$'\n'/,}" "${ncbi_api_key}")"
            acc_summary="$(echo "${xml_summary}" | grep -oP '(?<=Name="AccessionVersion" Type="String">)[^<]+')"
            # if no accession returned, try again
            if [[ -z "${acc_summary}" ]]; then 
                #(>&2 printf "Continue summary ${i}\n")
                sleep ${i} # wait to not oveload ncbi server
                continue
            else
                len_summary="$(echo "${xml_summary}" | grep -oP '(?<=Name="Length" Type="Integer">)[^<]+')"
                taxid_summary="$(echo "${xml_summary}" | grep -oP '(?<=Name="TaxId" Type="Integer">)[^<]+')"
                out="$(paste <(echo "${acc_summary}") <(echo "${len_summary}") <(echo "${taxid_summary}") --delimiters '\t')"
                break
            fi
        done
        # if there are less output lines than input accessions, get accessions difference (or return is empty)
        if [[ "$(echo "${acc_summary}" | wc -l)" -lt "$(echo "${acc}" | wc -l)" ||  -z "${acc_summary}"  ]];
        then
            acc_diff="$(diff --changed-group-format='%>' --unchanged-group-format='' <(echo "${acc_summary}") <(echo "$acc"))"
        else
            acc_diff=""
        fi

        # If there are accessions left
        if [[ ! -z "${acc_diff}" ]]; then  
            for i in $(seq 1 ${att});
            do
                # try another method
                xml_fetch="$(retrieve_nucleotide_fasta_xml "${acc_diff//$'\n'/,}" "${ncbi_api_key}")"
                acc_fetch="$(echo "${xml_fetch}" | grep -oP '(?<=<TSeq_accver>)[^<]+')"
                # if no accession returned, try again
                if [[ -z "${acc_fetch}" ]]; then 
                    #(>&2 printf "Continue fetch ${i}\n")
                    sleep ${i} # wait to not oveload ncbi server
                    continue
                else
                    len_fetch="$(echo "${xml_fetch}" | grep -oP '(?<=<TSeq_length>)[^<]+')"
                    taxid_fetch="$(echo "${xml_fetch}" | grep -oP '(?<=<TSeq_taxid>)[^<]+')"
                    if [[ ! -z "${out}" ]]; then 
                        out="${out}"$'\n'
                    fi
                    out="${out}$(paste <(echo "${acc_fetch}") <(echo "${len_fetch}") <(echo "${taxid_fetch}") --delimiters '\t')"
                    break
                fi
            done
            # if there are less output lines than input accessions, get accessions missing (or return is empty)
            if [[ "$(echo "${acc_fetch}" | wc -l)" -lt "$(echo "${acc_diff}" | wc -l)" ||  -z "${acc_fetch}" ]];
            then
                acc_missing="$(diff --changed-group-format='%>' --unchanged-group-format='' <(echo "${acc_fetch}") <(echo "$acc_diff"))"
                failed="${failed}${acc_missing}\n"
            fi
        fi

        if [ "${keep_all}" -eq 1 ]; then
            # Keep all output lines
            out="$(join -1 1 -2 1 <(echo "${acc}" | sort -k 1,1 ) <(echo "${out}" | sort -k 1,1) -t$'\t' -o "1.1,2.2,2.3" -a 1 -e "na")"           
        fi

    else
        # print all accessions to out
        out="${acc}"
    fi

    # if should retrieve assembly accessions/names
    if [ "${get_assembly_accession}" -eq 0 ] && [ "${get_assembly_name}" -eq 0 ]; then
        # Print results sorted
        echo "$(sort_in_out "${out}" "${acc}")"
    else
        # Get assembly uid
        acc_retrieved="$(echo "${out}" | cut -f 1)"

        uid_link=""
        acc_link=""
        acc_uid_link=""
        if [[ ! -z "${acc_retrieved}" ]]; then 
            for i in $(seq 1 ${att});
            do
                xml_link="$(retrieve_assembly_uid_xml "${acc_retrieved//$'\n'/&id=}" "${ncbi_api_key}")"
                # request with several &id= instead of comma separated to get in order
                all_id_link="$(echo "${xml_link}" | tr -d '\n' | grep -oP '(?<=<LinkSet>).*?(?=</LinkSet>)' | tr -d ' ' | tr -d '\t')"
                
                # link accession with link containing id, remove containing errors, just pass containing desired tag
                # assuming that the results are returning on the same order of the input
                out_link="$(paste <(echo "$acc_retrieved") <(echo "${all_id_link}") --delimiters '\t' | grep -v "ERROR" | grep "</Id></Link></LinkSetDb>" )"

                # get accessions that were retrieved
                acc_link="$(echo "${out_link}" | cut -f 1)"
                
                # if no uids were retrieved
                if [[ -z "${acc_link}" ]]; then 
                    continue
                else
                    # Get only first id inside LinkSetDb of link nuccore_assembly (can have several assemblies)
                    # There is no way to know what is the latest or correct liked assembly accession for this sequence accession
                    # Later the latest assembly accession is chosen
                    uid_link="$(echo "${out_link}" | cut -f 2 | grep --color -oP '(?<=<LinkName>nuccore_assembly</LinkName><Link><Id>)[0-9]*?(?=</Id>)')"
                    acc_uid_link="$(paste <(echo "${acc_link}") <(echo "${uid_link}") --delimiters '\t')"
                    break
                fi
            done

            # if there are less output lines than input accessions, get accessions missing (or return is empty)
            if [[ "$(echo "${acc_link}" | wc -l)" -lt "$(echo "${acc_retrieved}" | wc -l)" || -z "${acc_link}" ]];
            then
                acc_missing="$(diff --changed-group-format='%>' --unchanged-group-format='' <(echo "${acc_link}") <(echo "$acc_retrieved"))"
                failed_assembly="${failed_assembly}${acc_missing}\n"
            fi
        fi

        # if managed to get uid for assembly
        uid_summary_assembly=""
        assemblyaccession_summary_assembly=""
        uid_assemblyaccession_summary_assembly=""
        if [[ ! -z "${uid_link}" ]]; then 
            for i in $(seq 1 ${att});
            do
                # Retrieve assembly accessions
                xml_summary_assembly="$(retrieve_assembly_accession_xml "${uid_link//$'\n'/,}" "${ncbi_api_key}")"
                uid_summary_assembly="$(echo "${xml_summary_assembly}" | grep -oP '(?<=DocumentSummary uid=\")[^\"]+')"
                if [[ -z "${uid_summary_assembly}" ]]; then 
                    continue
                else
                    if [ "${get_assembly_accession}" -eq 1 ]; then 
                        # Get latest and current assembly accession
                        latest_assemblyaccession_summary_assembly="$(echo "${xml_summary_assembly}" | grep -oP '(?<=<LatestAccession)[^<]+')" 
                        current_assemblyaccession_summary_assembly="$(echo "${xml_summary_assembly}" | grep -oP '(?<=<AssemblyAccession>)[^<]+')"
                        # choose always latest assembly if present, since there is no way to link exact assembly accession to sequence with eutils
                        assemblyaccession_summary_assembly="$(paste <(echo "${current_assemblyaccession_summary_assembly}") <(echo "${latest_assemblyaccession_summary_assembly}" | sed "s/>//g") --delimiters '\t' | awk 'BEGIN{FS="\t"}{if($2){print $2}else{print $1}}')"
                        uid_assemblyaccession_summary_assembly="$(paste <(echo "${uid_summary_assembly}") <(echo "${assemblyaccession_summary_assembly}") --delimiters '\t')"
                    fi
                    if [ "${get_assembly_name}" -eq 1 ]; then 
                        assemblyname_summary_assembly="$(echo "${xml_summary_assembly}" | grep -oP '(?<=<Organism>)[^<]+')"
                        uid_assemblyname_summary_assembly="$(paste <(echo "${uid_summary_assembly}") <(echo "${assemblyname_summary_assembly}") --delimiters '\t')"
                    fi
                    break
                fi
            done
            # if there are less output lines than input accessions, get accessions missing (or return is empty)
            if [[ "$(echo "${uid_summary_assembly}" | wc -l)" -lt "$(echo "${uid_link}" | wc -l)" ||  -z "${uid_summary_assembly}" ]];
            then
                uid_missing="$(diff --changed-group-format='%>' --unchanged-group-format='' <(echo "${uid_summary_assembly}") <(echo "$uid_link"))"
                acc_missing="$(join -1 1 -2 2 <(echo "${uid_missing}" | sort | uniq) <(echo "${acc_uid_link}" | sort -k 2,2) -t$'\t' -o "2.1")"
                failed_assembly="${failed_assembly}${acc_missing}\n"
            fi
        fi

        final="${out}"
        # define output
        if [ "${get_length_taxid}" -eq 1 ]; then
            out_fields="0,1.2,1.3,2.2";
        else
            out_fields="0,2.2";
        fi

        if [ "${get_assembly_accession}" -eq 1 ]; then 
            # link uids (not found for esummary)
            acc_assemblyaccession="$(join -1 2 -2 1 <(echo "${acc_uid_link}" | sort -k 2,2) <(echo "${uid_assemblyaccession_summary_assembly}" | sort -k 1,1 | uniq) -t$'\t' -o "1.1,2.2" -a 1 -e PLACEHOLDER_NOT_FOUND)"
            # check for entries without assembly found (not found for elink)
            final="$(join -1 1 -2 1 <(echo "${final}" | sort -k 1,1) <(echo "${acc_assemblyaccession}" | sort -k 1,1) -t$'\t' -o ${out_fields} -a 1 -e PLACEHOLDER_NOT_FOUND)"
            
            # fix output for assembly name
            if [ "${get_assembly_name}" -eq 1 ]; then 
                if [ "${get_length_taxid}" -eq 1 ]; then
                    out_fields="0,1.2,1.3,1.4,2.2";
                else
                    out_fields="0,1.2,2.2";
                fi
            fi
        fi

        if [ "${get_assembly_name}" -eq 1 ]; then 
            #echo "${uid_assemblyname_summary_assembly}"
            acc_assemblyname="$(join -1 2 -2 1 <(echo "${acc_uid_link}" | sort -k 2,2) <(echo "${uid_assemblyname_summary_assembly}" | sort -k 1,1 | uniq) -t$'\t' -o "1.1,2.2" -a 1 -e PLACEHOLDER_NOT_FOUND)"
            #echo "${acc_assemblyname}"
            final="$(join -1 1 -2 1 <(echo "${final}" | sort -k 1,1) <(echo "${acc_assemblyname}" | sort -k 1,1) -t$'\t' -o ${out_fields} -a 1 -e PLACEHOLDER_NOT_FOUND)"
        fi

        # Print final
        if [ "${replace_not_found_accession}" -eq 1 ]; then
            final="$(awk 'BEGIN{FS="\t";OFS="\t"}{for(i=2; i<=NF; i++){if($i=="PLACEHOLDER_NOT_FOUND") $i=$1}; print}' <(echo "${final}"))"
        else
            final="$(awk 'BEGIN{FS="\t";OFS="\t"}{for(i=2; i<=NF; i++){if($i=="PLACEHOLDER_NOT_FOUND") $i="na"}; print}' <(echo "${final}"))"
        fi     

        # Print results sorted
        echo "$(sort_in_out "${final}" "${acc}")"

    fi

    # read new batch
    batch_count=$((batch_count+batch))
    acc="$(get_lines ${input_file} ${batch_count} ${batch})"
done

if [[ ! -z "${tmp_file}" ]]; then
    rm ${tmp_file}
fi

# Print errors to STDERR
if [[ ! -z "${failed}" ]];
then
    (>&2 printf "Failed to get taxid and sequence length for:\n${failed}")
fi
if [[ ! -z "${failed_assembly}" ]];
then
    (>&2 printf "Failed to get assembly accession/name for:\n${failed_assembly}")
fi
if [[ ! -z "${failed}" || ! -z "${failed_assembly}" ]];
then
    exit 1
else
    exit 0
fi
