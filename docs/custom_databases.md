# Custom databases

## Default NCBI assembly or sequence accession

Besides the automated download and build (`ganon build`) ganon provides a highly customizable build procedure (`ganon build-custom`) to create databases from local sequence files.

To use custom sequences, just provide them with `--input`. ganon will try to retrieve all necessary information necessary to build a database.

!!! note
    ganon expects assembly accessions in the filename like `GCA_002211645.1_ASM221164v1_genomic.fna.gz`. When using `--input-target sequence` filenames are not important but sequence headers should contain sequence accessions like `>CP022124.1 Fusobacterium nu...`. More information about building by file or sequence can be found [here](#target-file-or-sequence-input-target).

## Non-standard/custom accessions

It is also possible to use **non-standard accessions and headers** to build custom databases with `--input-file`. This file should contain the following fields (tab-separated): `file [<tab> target <tab> node <tab> specialization <tab> specialization_name].` Note that file is mandatory and additional fields not.

!!! tip
    If you just want to build a database without any taxonomic or target information, just sent the files with `--input`, use `--taxonomy skip` and choose between `--input-target file` or `sequence`.

<details>
  <summary>Examples of --input-file</summary>
  <br>

With --input-target file (default), where my_target_1 and my_target_2 are just names to assign sequences from (unique) sequence files:

```
sequences.fasta my_target_1
others.fasta my_target_2
```

With --input-target sequence, second column should match sequence headers on provided sequence files (that should be repeated for each header):

```
sequences.fasta HEADER1
sequences.fasta HEADER2
sequences.fasta HEADER3
others.fasta HEADER4
others.fasta HEADER5
```

A third column with taxonomic nodes can be provided to link the data with taxonomy. For example with --taxonomy ncbi:

```
sequences.fasta FILE_A  562
others.fasta FILE_B 623
```

```
sequences.fasta HEADER1 562
sequences.fasta HEADER2 562
sequences.fasta HEADER3 562
others.fasta HEADER4  623
others.fasta HEADER5  623
```

Further specializations can be used to create a additional classification level after the taxonomic leaves. For example (using --level custom):

```
sequences.fasta FILE_A  562 ID44444 Escherichia coli TW10119
others.fasta FILE_B 623 ID55555  Shigella flexneri 1a
```

```
sequences.fasta HEADER1 562 ID443 Escherichia coli TW10119
sequences.fasta HEADER2 562 ID297 Escherichia coli PCN079
sequences.fasta HEADER3 562 ID8873  Escherichia coli P0301867.7
others.fasta HEADER4  623 ID2241  Shigella flexneri 1a
others.fasta HEADER5  623 ID4422  Shigella flexneri 1b
```

</details>
<br>

## Examples

Some examples with download and build commands for custom ganon databases from useful and commonly used repositories and datasets for metagenomics analysis:

### HumGut

Collection of >30000 genomes from healthy human metagenomes. [Article](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-021-01114-w)/[Website](https://arken.nmbu.no/~larssn/humgut/).

```bash
# Download sequence files
wget --quiet --show-progress "http://arken.nmbu.no/~larssn/humgut/HumGut.tar.gz"
tar xf HumGut.tar.gz 

# Download taxonomy and metadata files
wget "https://arken.nmbu.no/~larssn/humgut/ncbi_nodes.dmp"
wget "https://arken.nmbu.no/~larssn/humgut/ncbi_names.dmp"
wget "https://arken.nmbu.no/~larssn/humgut/HumGut.tsv"
# Generate --input-file from metadata
tail -n+2 HumGut.tsv | awk -F"\t" '{print "fna/"$21"\t"$1"\t"$2}' > HumGut_ganon_input_file.tsv

# Build ganon database
ganon build-custom --input-file HumGut_ganon_input_file.tsv --taxonomy-files ncbi_nodes.dmp ncbi_names.dmp --db-prefix HumGut --level strain --threads 32
```

Similarly using GTDB taxonomy files:

```bash
# Download taxonomy files
wget "https://arken.nmbu.no/~larssn/humgut/gtdb_nodes.dmp"
wget "https://arken.nmbu.no/~larssn/humgut/gtdb_names.dmp"

# Build ganon database
ganon build-custom --input-file HumGut_ganon_input_file.tsv --taxonomy-files gtdb_nodes.dmp gtdb_names.dmp --db-prefix HumGut_gtdb --level strain --threads 32
```

!!! note
    There is no need to use ganon's gtdb integration here since GTDB files in NCBI format are available

### Plasmid, Plastid and Mitochondrion from RefSeq

Extra repositories from RefSeq release not included as default databases. [Website](https://www.ncbi.nlm.nih.gov/refseq/).

```bash
# Download sequence files
wget -A genomic.fna.gz -m -nd --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/"
wget -A genomic.fna.gz -m -nd --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/"
wget -A genomic.fna.gz -m -nd --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/"

# Build ganon database
ganon build-custom --input plasmid.* plastid.* mitochondrion.* --db-prefix ppm --input-target sequence --level leaves --threads 32 
```

### UniVec, UniVec_core

"UniVec is a non-redundant database of sequences commonly attached to cDNA or genomic DNA during the cloning process." [Website](https://ftp.ncbi.nlm.nih.gov/pub/UniVec/README.uv). Useful to screen for vector and linker/adapter contamination. UniVec_core is a sub-set of the UniVec selected to reduce the false positive hits from real biological sources.

```bash
# UniVec
wget -O "UniVec.fasta" --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec" 
grep -o '^>[^ ]*' UniVec.fasta | sed 's/^>//' | awk '{print "UniVec.fasta\t"$1"\t81077"}' > UniVec_ganon_input_file.tsv
ganon build-custom --input-file UniVec_ganon_input_file.tsv --db-prefix UniVec --input-target sequence --level leaves --threads 8

# UniVec_Core
wget -O "UniVec_Core.fasta" --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core" 
grep -o '^>[^ ]*' UniVec_Core.fasta | sed 's/^>//' | awk '{print "UniVec_Core.fasta\t"$1"\t81077"}' > UniVec_Core_ganon_input_file.tsv
ganon build-custom --input-file UniVec_Core_ganon_input_file.tsv --db-prefix UniVec_Core --input-target sequence --level leaves --threads 8
```

!!! note
    All UniVec entries in the examples are categorized as [Artificial Sequence](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=81077) (NCBI txid:81077). Some are completely artificial but others may be derived from real biological sources. More information in this [link](https://ftp.ncbi.nlm.nih.gov/pub/UniVec/README.vector.origins).

### MGnify genome catalogues (MAGs)

"Genome catalogues are biome-specific collections of metagenomic-assembled and isolate genomes". [Article](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/)/[Website](https://www.ebi.ac.uk/metagenomics/browse/genomes)/[FTP](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/).

There are currently (2023-05-04) 8 genome catalogues available: chicken-gut, human-gut, human-oral, marine, non-model-fish-gut, pig-gut and zebrafish-fecal. An example below how to download and build the human-oral catalog:


```bash
# Download metadata
wget "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-oral/v1.0/genomes-all_metadata.tsv"

# Download sequence files with 12 threads
tail -n+2 genomes-all_metadata.tsv | cut -f 1,20 | xargs -P 12 -n2 sh -c 'curl --silent ${1}| gzip -d | sed -e "1,/##FASTA/ d" | gzip > ${0}.fna.gz'

# Generate ganon input file
tail -n+2 genomes-all_metadata.tsv | cut -f 1,15  | tr ';' '\t' | awk -F"\t" '{tax="1";for(i=NF;i>1;i--){if(length($i)>3){tax=$i;break;}};print $1".fna.gz\t"$1"\t"tax}' > ganon_input_file.tsv

# Build ganon database
ganon build-custom --input-file ganon_input_file.tsv --db-prefix mgnify_human_oral_v1 --taxonomy gtdb --level leaves --threads 32
```

!!! note
    MGnify genomes catalogues will be build with GTDB taxonomy.

### Pathogen detection FDA-ARGOS

A collection of >1400 "microbes that include biothreat microorganisms, common clinical pathogens and closely related species". [Article](https://www.ncbi.nlm.nih.gov/bioproject/231221)/[Website](https://www.ncbi.nlm.nih.gov/bioproject/231221)/[BioProject](https://www.ncbi.nlm.nih.gov/bioproject/231221).

```bash
# Download sequence files
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt
grep "strain=FDAARGOS" assembly_summary_refseq.txt > fdaargos_assembly_summary.txt
genome_updater.sh -e fdaargos_assembly_summary.txt -f "genomic.fna.gz" -o download -m -t 12

# Build ganon database
ganon build-custom --input download/ --input-recursive --db-prefix fdaargos --ncbi-file-info download/assembly_summary.txt --level assembly --threads 32
```

!!! note
    The example above uses [genome_updater](https://github.com/pirovc/genome_updater) to download files

### BLAST databases (nt, env_nt, nt_prok, ...)

BLAST databases. [Website](https://blast.ncbi.nlm.nih.gov/Blast.cgi)/[FTP](https://ftp.ncbi.nlm.nih.gov/blast/db/).

Current available nucleotide databases (2023-05-04): `16S_ribosomal_RNA` `18S_fungal_sequences` `28S_fungal_sequences` `Betacoronavirus` `env_nt` `human_genome` `ITS_eukaryote_sequences` `ITS_RefSeq_Fungi` `LSU_eukaryote_rRNA` `LSU_prokaryote_rRNA` `mito` `mouse_genome` `nt` `nt_euk` `nt_others` `nt_prok` `nt_viruses` `patnt` `pdbnt` `ref_euk_rep_genomes` `ref_prok_rep_genomes` `refseq_rna` `refseq_select_rna` `ref_viroids_rep_genomes` `ref_viruses_rep_genomes` `SSU_eukaryote_rRNA` `tsa_nt`

!!! note
    List currently available nucleotide databases `curl --silent --list-only ftp://ftp.ncbi.nlm.nih.gov/blast/db/ | grep "nucl-metadata.json" | sed 's/-nucl-metadata.json/, /g' | sort`

!!! warning
    Some BLAST databases are very big and may require extreme computational resources to build. You may need to use some [reduction strategies](../default_databases/#reducing-database-size)

The example below extracts sequences and information from a BLAST db to build a ganon database:

```bash
# Define BLAST db
db="16S_ribosomal_RNA"
threads=8

# Download BLAST db - re-run this command many times until all finish (no more output)
curl --silent --list-only ftp://ftp.ncbi.nlm.nih.gov/blast/db/ | grep "^${db}\..*tar.gz$" | xargs -P ${threads:-1} -I{} wget --continue -nd --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/blast/db/{}"

# OPTIONAL Download and check MD5
wget -O - -nd --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/blast/db/${db}\.*tar.gz.md5" > "${db}.md5"
md5sum "${db}".*tar.gz > "${db}_downloaded.md5"
diff -s <(sort -k 2,2 "${db}.md5") <(sort -k 2,2 "${db}_downloaded.md5")  # Should print "Files /dev/fd/xx and /dev/fd/xx are identical"

# Extract BLAST db files, if successful, remove .tar.gz
ls "${db}"*.tar.gz | xargs -P ${threads} -I{} sh -c 'gzip -dc {} | tar --overwrite -vxf - && rm {}' > "${db}_extracted_files.txt"

# Create folder to write sequence files (split into 10 sub-folders)
seq 0 9 | xargs -i mkdir -p "${db}"/{}

# This command extracts sequences from the blastdb and writes them into taxid specific files
# It also generates the --input-file for ganon
blastdbcmd -entry all -db "${db}" -outfmt "%a %T %s" | \
awk -v db="${db}" '{file=db"/"substr($2,1,1)"/"$2".fna"; print ">"$1"\n"$3 >> file; print file"\t"$2"\t"$2}' | \
sort | uniq > "${db}_ganon_input_file.tsv"

# Build ganon database
ganon build-custom --input-file "${db}_ganon_input_file.tsv" --db-prefix "${db}" --level species --threads 12

# Delete extracted files and sequences
cat "${db}_extracted_files.txt" | xargs rm
rm "${db}_extracted_files.txt" "${db}_ganon_input_file.tsv" "${db}.md5" "${db}_downloaded.md5"
rm -rf "${db}"
```

!!! note
    `blastdbcmd` is a command from BLAST software suite and should be installed separately

### Files from genome_updater

To create a ganon database from files previosly downloaded with [genome_updater](https://github.com/pirovc/genome_updater):

```bash
ganon build-custom --input output_folder_genome_updater/version/ --input-recursive --db-prefix mydb --ncbi-file-info  output_folder_genome_updater/assembly_summary.txt --level assembly --threads 32
```

## Parameter details

### False positive and size (--max-fp, --filter-size)

ganon indices are based on bloom filters and can have false positive matches. This can be controlled with `--max-fp` parameter. The lower the `--max-fp`, the less chances of false positives matches on classification, but the larger the database size will be. For example, with `--max-fp 0.01` the database will be build so any target (defined by `--level`) will have 1 in a 100 change of reporting a false k-mer match. The false positive of the query (all k-mers of the reads) is higher but directly affected.

Alternatively, one can set a specific size for the final index with `--filter-size`. When using this option, please observe the theoretic false positive of the index reported at the end of the building process.

### minimizers (--window-size, --kmer-size)

in `ganon build`, when `--window-size` > `--kmer-size` minimizers are used. That means that for a every window, a single k-mer will be selected. It produces smaller database files and requires substantially less memory overall. It may increase building times but will have a huge benefit for classification times. Sensitivity and precision can be reduced by small margins. If `--window-size` = `--kmer-size`, all k-mers are going to be used to build the database.

### Target file or sequence (--input-target) 

Customized builds can be done either by file or sequence. `--input-target file` will consider every file provided with `--input` a single unit. `--input-target sequence` will use every sequence as a unit.

`--input-target file` is the default behavior and most efficient way to build databases. `--input-target sequence` should only be used when the input sequences are stored in a single file or when classification at sequence level is desired.

### Build level (--level)

The `--level` parameter defines the max. depth of the database for classification. This parameter is relevant because the `--max-fp` is going to be guaranteed at the `--level` chosen. By default, the level will be the same as `--input-target`, meaning that classification will be done either at file or sequence level.

Alternatively, `--level assembly` will link the file or sequence target information with assembly accessions retrieved from NCBI servers. `--level leaves` or `--level species` (or genus, family, ...) will link the targets with taxonomic information and prune the tree at the chosen level. `--level custom` will use specialization level define in the `--input-file`.

### Genome sizes (--genome-size-files)

Ganon will automatically download auxiliary files to define an approximate genome size for each entry in the taxonomic tree. For `--taxonomy ncbi` the [species_genome_size.txt.gz](https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/) is used. For `--taxonomy gtdb` the [\*_metadata.tar.gz](https://data.gtdb.ecogenomic.org/releases/latest/) files are used. Those files can be directly provided with the `--genome-size-files` argument.

Genome sizes of parent nodes are calculated as the average of the respective children nodes. Other nodes without direct assigned genome sizes will use the closest parent with a pre-calculated genome size. The genome sizes are stored in the [ganon database](#buildupdate).

### Retrieving info (--ncbi-sequence-info, --ncbi-file-info)

Further taxonomy and assembly linking information has to be collected to properly build the database. `--ncbi-sequence-info` and `--ncbi-file-info` allow customizations on this step.

When `--input-target sequence`, `--ncbi-sequence-info` argument allows the use of NCBI e-utils webservices (`eutils`) or downloads accession2taxid files to extract target information (options `nucl_gb` `nucl_wgs` `nucl_est` `nucl_gss` `pdb` `prot` `dead_nucl` `dead_wgs` `dead_prot`). By default, ganon uses `eutils` up-to 50000 input sequences, otherwise it downloads `nucl_gb` `nucl_wgs` from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/. Previously downloaded files can be directly provided with this argument.

When `--input-target file`, `--ncbi-file-info` uses `assembly_summary.txt` from https://ftp.ncbi.nlm.nih.gov/genomes/ to extract target information (options `refseq` `genbank` `refseq_historical` `genbank_historical`. Previously downloaded files can be directly provided with this argument.

If you are using outdated, removed or inactive assembly or sequence files and accessions from NCBI, make sure to include `dead_nucl` `dead_wgs` for `--ncbi-sequence-info` or `refseq_historical` `genbank_historical` for `--ncbi-file-info`. `eutils` option does not work with outdated accessions.
