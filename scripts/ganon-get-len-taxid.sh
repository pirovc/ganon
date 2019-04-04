#!/bin/bash
# Number of attempts to request data from e-utils
att=10
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
	echo "$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=assembly&id=${1}${2}")"
}

retrieve_assembly_accession_xml()
{
	echo "$(curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id=${1}${2}")"
}

get_lines()
{
	echo "$(sed -n "${2},$((${2}+${3}-1))p" ${1})"
}

function showhelp {
	echo "get_len_taxid.sh by Vitor C. Piro (vitorpiro@gmail.com, http://github.com/pirovc)"
	echo
	echo $' -i [str] input file with accessions (use - to read from STDIN)' 
	echo $' -l [str] list of accesions (space separated)'
	echo $' -n [str] ncbi_api_key'
	echo $' -a get_assembly_accession'
	echo
}

input_file=""
list_acc=""
ncbi_api_key=""
get_assembly_accession=0

OPTIND=1 # Reset getopts
while getopts "i:l:n:a" opt; do
  case ${opt} in
    i) input_file=${OPTARG} ;;
    l) list_acc=${OPTARG} ;;
    n) ncbi_api_key=${OPTARG} ;;
    a) get_assembly_accession=1 ;;
    h|\?) showhelp; exit 1 ;;
    :) echo "Option -${OPTARG} requires an argument." >&2; exit 1 ;;
  esac
done
if [ ${OPTIND} -eq 1 ]; then showhelp; exit 1; fi
shift $((OPTIND-1))
[ "$1" = "--" ] && shift

if [[ "${input_file}" == "-" ]]; then
	tmp_file="$(mktemp)"
	cat /dev/stdin > ${tmp_file}
	input_file=${tmp_file}
elif [[ ! -z "${list_acc}" ]]; then
	tmp_file="$(mktemp)"
	echo -e "${list_acc// /$'\n'}" > ${tmp_file}
	input_file=${tmp_file}
fi

# test if ncbi api key is valid
if [[ ! -z "${ncbi_api_key}" ]]; then
	ncbi_api_key="&api_key=${ncbi_api_key}"
	api_key_test="$(retrieve_summary_xml "" "${ncbi_api_key}" | grep -o 'API key invalid')"
	if [[ ! -z "${api_key_test}" ]]; then
		ncbi_api_key=""
	fi
fi

batch_count=1
acc="$(get_lines ${input_file} 1 ${batch})"
while [[ ! -z "${acc}" ]];
do
	out=""

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
	        #(>&2 printf "Fetch\n")	
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

	# if should retrieve assembly accessions
	if [ "${get_assembly_accession}" -eq 0 ]; 
	then
		echo "${out}"
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
					# Get only last Id inside LinkSetDb (can have several assemblies)
					uid_link="$(echo "${out_link}" | cut -f 2 | grep --color -oP '(?<=<Id>)[0-9]*?(?=</Id></Link></LinkSetDb>)')"
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
					assemblyaccession_summary_assembly="$(echo "${xml_summary_assembly}" | grep -oP '(?<=<AssemblyAccession>)[^<]+')"
					uid_assemblyaccession_summary_assembly="$(paste <(echo "${uid_summary_assembly}") <(echo "${assemblyaccession_summary_assembly}") --delimiters '\t')"
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

		# link uids (not found for esummary)
		acc_assembly="$(join -1 2 -2 1 <(echo "${acc_uid_link}" | sort -k 2,2) <(echo "${uid_assemblyaccession_summary_assembly}" | sort -k 1,1 | uniq) -t$'\t' -o "1.1,2.2" -a 1 -e NOT_FOUND)"
		# print final (not found for elink)
		echo "$(join -1 1 -2 1 <(echo "${out}" | sort -k 1,1) <(echo "${acc_assembly}" | sort -k 1,1) -t$'\t' -o "0,1.2,1.3,2.2" -a 1 -e NOT_FOUND)"
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
	(>&2 printf "Failed to get taxid and sequence length:\n${failed}")
fi
if [[ ! -z "${failed_assembly}" ]];
then
	(>&2 printf "Failed to get assembly accession:\n${failed_assembly}")
fi
if [[ ! -z "${failed}" || ! -z "${failed_assembly}" ]];
then
	exit 1
else
	exit 0
fi
