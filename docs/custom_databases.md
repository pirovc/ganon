# Custom databases

Besides the automated download and build (`ganon build`) ganon provides a highly customizable build procedure (`ganon build-custom`) to create databases from local sequence files. The usage of this procedure depends on the configuration of your files:

- Filename like `GCA_002211645.1_ASM221164v1_genomic.fna.gz`: genomic fasta files in the NCBI standard, with assembly accession in the beginning of the filename. Provide the files with the `--input` parameter. ganon will try to retrieve all necessary information to build the database.
- Headers like `>NC_006297.1 Bacteroides fragilis YCH46 ...`: sequence headers are in the NCBI standard, with sequence accession in after `>` and with a space afterwards (or line break). Provide the files with the `--input` parameter and set `--input-target sequence`. ganon will try to retrieve all necessary information to build the database.
- For non-standard filenames and headers, [follow this](#non-standard-filesheaders-with-input-file)

!!! warning
    `--input-target sequence` will be slower to build and will use more disk space, since files have be re-written separately for each sequence. More information about building by file or sequence can be found [here](#target-file-or-sequence-input-target).

The `--level` is a important parameter that will define the (max.) classification level for the database ([more infos](#build-level-level)):

- `--level file` or `sequence` -> default behavior (depending on `--input-target `), use file/sequence as classification target
- `--level assembly` -> will retrieve assembly related to the file/sequence, use assembly as classification target
- `--level leaves` or `species`,`genus`,... -> group input by taxonomy, use tax. nodes at the rank chosen as classification target

More infos about other parameters [here](#parameter-details).

## Non-standard files/headers with `--input-file`

Alternatively to the automatic input methods, it is possible to manually define the input with either standard or **non-standard filenames, accessions and headers** to build custom databases with `--input-file`. This file should contain the following fields (tab-separated):

`file [<tab> target <tab> node <tab> specialization <tab> specialization_name].`

- `file`: relative or full path to the sequence file
- `target`: any unique text to name the file, to be used in the taxonomy
- `node`: taxonomic node (e.g. taxid) to link entry with taxonomy
- `specialization`: creates a specialized taxonomic level with a custom name, allowing files to be grouped
- `specialization_name`: a name for the specialization, to be used in the taxonomy

!!! warning
    the `target` and `specialization` fields (2nd and 4th col) cannot be the same as the `node` (3rd col)

Below you find example of `--input-file`. Note they are slightly different depending on the `--input-target` chosen. They need to be *tab-separated* to be properly parsed (tsv).

### Examples of `--input-file` using the default `--input-target file`

#### List of files

```
sequences.fasta
others.fasta
```

No taxonomic information is provided so `--taxonomy skip` should be set. The classification against the generated database will be performed at file level (`--level file`), since that is the only available information given.

#### List of files with alternative names
```
sequences.fasta  sequences
others.fasta     others
```

Just like above, but with a specific name to be used for each file.

#### Files and taxonomy

```
sequences.fasta  sequences  562
others.fasta     others     623
```

The classification max. level against this database will depend on the value set for `--level`:

- `--level file` -> use the file (named with target) with node as parent
- `--level leaves` or `species`,`genus`,... -> files are grouped by taxonomy

#### Files, taxonomy and specialization

```
sequences.fasta  sequences  562  ID44444  Escherichia coli TW10119
others.fasta     others     623  ID55555  Shigella flexneri 1a
```

The classification max. level against this database will depend on the value set for `--level`:

- `--level custom` -> use the specialization (named with specialization_name) with node as parent
- `--level file` -> use the file (named with target) as a tax. node as parent
- `--level leaves` or `species`,`genus`,... -> files are grouped by taxonomy

### Examples of `--input-file` using `--input-target sequence`

To provide a tabular information for every sequence in your files, you need to use the `target` field (2nd col.) of the `--input-file` to input sequence headers. For example:

#### Sequences and taxonomy

```
sequences.fasta  NZ_CP054001.1  562
sequences.fasta  NZ_CP117955.1  623
others.fasta     header1        666
others.fasta     header2        666
```

The classification max. level against this database will depend on the value set for `--level`:

- `--level sequence` -> use the sequence header with node as parent
- `--level assembly` -> will attempt to retrieve the assembly related to the sequence with node as parent
- `--level leaves` or `species`,`genus`,... -> files are grouped by taxonomy

#### Sequences, taxonomy and specialization

```
sequences.fasta  NZ_CP054001.1  562  ID44444  Escherichia coli TW10119
sequences.fasta  NZ_CP117955.1  623  ID55555  Shigella flexneri 1a
others.fasta     header1        666  StrainA  My Strain
others.fasta     header2        666  StrainA  My Strain
```

The classification max. level against this database will depend on the value set for `--level`:

- `--level custom` -> use the specialization (named with specialization_name) with node as parent
- `--level sequence` -> use the sequence header with node as parent
- `--level leaves` or `species`,`genus`,... -> files are grouped by taxonomy

## Examples

Below you will find some examples from commonly used repositories for metagenomics analysis with `ganon build-custom`: 

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
tail -n+2 HumGut.tsv | awk -F"\t" '{print "fna/"$21"\t"$1"\t"$6}' > HumGut_ganon_input_file.tsv

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

ganon build-custom --input plasmid.* plastid.* mitochondrion.* --db-prefix ppm --level species --threads 8 --input-target sequence
```

### UniVec, UniVec_core

"UniVec is a non-redundant database of sequences commonly attached to cDNA or genomic DNA during the cloning process." [Website](https://ftp.ncbi.nlm.nih.gov/pub/UniVec/README.uv). Useful to screen for vector and linker/adapter contamination. UniVec_core is a sub-set of the UniVec selected to reduce the false positive hits from real biological sources.

```bash
# UniVec
wget -O "UniVec.fasta" --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec" 
echo -e "UniVec.fasta\tUniVec\t81077" > UniVec_ganon_input_file.tsv
ganon build-custom --input-file UniVec_ganon_input_file.tsv --db-prefix UniVec --level leaves --threads 8 --skip-genome-size

# UniVec_Core
wget -O "UniVec_Core.fasta" --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core" 
echo -e "UniVec_Core.fasta\tUniVec_Core\t81077" > UniVec_Core_ganon_input_file.tsv
ganon build-custom --input-file UniVec_Core_ganon_input_file.tsv --db-prefix UniVec_Core --level leaves --threads 8 --skip-genome-size
```

!!! note
    All UniVec entries in the examples are categorized as [Artificial Sequence](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=81077) (NCBI txid:81077). Some are completely artificial but others may be derived from real biological sources. More information in this [link](https://ftp.ncbi.nlm.nih.gov/pub/UniVec/README.vector.origins).

### MGnify genome catalogues (MAGs)

"Genome catalogues are biome-specific collections of metagenomic-assembled and isolate genomes". [Article](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/)/[Website](https://www.ebi.ac.uk/metagenomics/browse/genomes)/[FTP](https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/).

Currently available genome catalogues (2024-02-09): `chicken-gut` `cow-rumen` `honeybee-gut` `human-gut` `human-oral` `human-vaginal` `marine` `mouse-gut` `non-model-fish-gut` `pig-gut` `zebrafish-fecal`

*List currently available entries `curl --silent --list-only ftp://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/`*

Example on how to download and build the `human-oral` catalog:

```bash
# Download metadata
wget "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-oral/v1.0.1/genomes-all_metadata.tsv"

# Download sequence files with 12 threads
tail -n+2 genomes-all_metadata.tsv | cut -f 1,20 | xargs -P 12 -n2 sh -c 'curl --silent ${1}| gzip -d | sed -e "1,/##FASTA/ d" | gzip > ${0}.fna.gz'

# Generate ganon input file
tail -n+2 genomes-all_metadata.tsv | cut -f 1,15  | tr ';' '\t' | awk -F"\t" '{tax="1";for(i=NF;i>1;i--){if(length($i)>3){tax=$i;break;}};print $1".fna.gz\t"$1"\t"tax}' > ganon_input_file.tsv

# Build ganon database
ganon build-custom --input-file ganon_input_file.tsv --db-prefix mgnify_human_oral_v101 --taxonomy gtdb --level leaves --threads 8
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

### BLAST databases (nt env_nt nt_prok ...)

BLAST databases. [Website](https://blast.ncbi.nlm.nih.gov/Blast.cgi)/[FTP](https://ftp.ncbi.nlm.nih.gov/blast/db/).

Current available nucleotide databases (2024-02-09): `16S_ribosomal_RNA` `18S_fungal_sequences` `28S_fungal_sequences` `Betacoronavirus` `env_nt` `human_genome` `ITS_eukaryote_sequences` `ITS_RefSeq_Fungi` `LSU_eukaryote_rRNA` `LSU_prokaryote_rRNA` `mito` `mouse_genome` `nt` `nt_euk` `nt_others` `nt_prok` `nt_viruses` `patnt` `pdbnt` `ref_euk_rep_genomes` `ref_prok_rep_genomes` `refseq_rna` `refseq_select_rna` `ref_viroids_rep_genomes` `ref_viruses_rep_genomes` `SSU_eukaryote_rRNA` `tsa_nt`

*List currently available entries `curl --silent --list-only ftp://ftp.ncbi.nlm.nih.gov/blast/db/ | grep "nucl-metadata.json" | sed 's/-nucl-metadata.json/, /g' | sort`*

!!! warning
    Some BLAST databases are very big and may require extreme computational resources to build. You may need to use some [reduction strategies](default_databases.md#reducing-database-size).

The example shows how to **download**, **parse** and **build** a ganon database from BLAST database files. It does so by splitting the database into taxonomic specific files, to speed-up the build process:

```bash
# Define BLAST db
db="16S_ribosomal_RNA"
threads=8

# Download BLAST db - re-run this command many times until all finish (no more output)
curl --silent --list-only ftp://ftp.ncbi.nlm.nih.gov/blast/db/ | grep "^${db}\..*tar.gz$" | xargs -P ${threads:-1} -I{} wget --continue -nd --quiet --show-progress "https://ftp.ncbi.nlm.nih.gov/blast/db/{}"

# OPTIONAL Download and check MD5
wget -O - -nd --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/blast/db/${db}\.*tar.gz.md5" > "${db}.md5"
find -name "${db}.*tar.gz" -type f -printf '%P\n' | xargs -P ${threads:-1} -I{} md5sum {} > "${db}_downloaded.md5"
diff -sy <(sort -k 2,2 "${db}.md5") <(sort -k 2,2 "${db}_downloaded.md5")  # Should print "Files /dev/fd/xx and /dev/fd/xx are identical"

# Extract BLAST db files, if successful, remove .tar.gz
find -name "${db}.*tar.gz" -type f -printf '%P\n' | xargs -P ${threads} -I{} sh -c 'gzip -dc {} | tar --overwrite -vxf - && rm {}' > "${db}_extracted_files.txt"

# Create folder to write sequence files (split into 10 sub-folders)
seq 0 9 | xargs -i mkdir -p "${db}"/{}

# This command extracts sequences from the blastdb and writes them into taxid specific files
# It also generates the --input-file for ganon with the fields: filepath <tab> file <tab> taxid
blastdbcmd -entry all -db "${db}" -outfmt "%a %T %s" | \
awk -v db="$(realpath ${db})" '{file=db"/"substr($2,1,1)"/"$2".fna"; print ">"$1"\n"$3 >> file; print file"\t"$2".fna\t"$2}' | \
sort | uniq > "${db}_ganon_input_file.tsv"

# Build ganon database
ganon build-custom --input-file "${db}_ganon_input_file.tsv" --db-prefix "${db}" --threads ${threads} --level leaves

# Delete extracted files and auxiliary files
cat "${db}_extracted_files.txt" | xargs rm
rm "${db}_extracted_files.txt" "${db}.md5" "${db}_downloaded.md5"
# Delete sequences and input_file
rm -rf "${db}" "${db}_ganon_input_file.tsv"
```

!!! note
    `blastdbcmd` is a command from [BLAST+ software suite](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (tested version 2.14.0) and should be installed separately.

### Files from genome_updater

To create a ganon database from files previously downloaded with [genome_updater](https://github.com/pirovc/genome_updater):

```bash
ganon build-custom --input output_folder_genome_updater/version/ --input-recursive --db-prefix mydb --ncbi-file-info  output_folder_genome_updater/assembly_summary.txt --level assembly --threads 32
```

## Parameter details

### False positive and size (--max-fp, --filter-size)

ganon indices are based on bloom filters and can have false positive matches. This can be controlled with `--max-fp` parameter. The lower the `--max-fp`, the less chances of false positives matches on classification, but the larger the database size will be. For example, with `--max-fp 0.01` the database will be build so any target (defined by `--level`) will have 1 in a 100 change of reporting a false k-mer match. [The false positive of the query](classification.md#false-positive-of-a-query-fpr-query) (all k-mers of a read) will be way lower, but directly affected by this value.

Alternatively, one can set a specific size for the final index with `--filter-size`. When using this option, please observe the theoretic false positive of the index reported at the end of the building process.

### minimizers (--window-size, --kmer-size)

in `ganon build`, when `--window-size` > `--kmer-size` minimizers are used. That means that for a every window, a single k-mer will be selected. It produces smaller database files and requires substantially less memory overall. It may increase building times but will have a huge benefit for classification times. Sensitivity and precision can be reduced by small margins. If `--window-size` = `--kmer-size`, all k-mers are going to be used to build the database.

### Target file or sequence (--input-target) 

This is a parameter that defines how ganon will parse your input files:
 - `--input-target file` (default) will consider every file provided with `--input` a single unit  (e.g. multi-fasta files are considered one input, sequence headers ignored).
 - `--input-target sequence` will use every sequence as a unit. For this, ganon will first decompose every sequence in the input files provided with `--input` into a separated file. This will take longer and use more disk space.

`--input-target file` is the default behavior and most efficient way to build databases. `--input-target sequence` should only be used when the input sequences are not separated by file (e.g. a single big FASTA file) or when classification at sequence level is desired.

### Build level (--level)

The `--level` parameter defines the max. depth of the database for classification. This parameter is relevant because the `--max-fp` is going to be guaranteed at the `--level` chosen. 

In `ganon build` the default value is `species`. In `ganon build-custom` the level will be the same as `--input-target`, meaning that classification will be done either at `file` or `sequence` level.

Alternatively, `--level assembly` will link the file or sequence target information with assembly accessions retrieved from NCBI. `--level leaves` or `--level species` (or `genus`, `family`, ...) will link the targets with taxonomic information and prune the tree at the chosen level. `--level custom` will use specialization (4th col.) defined in the `--input-file`.

### Genome sizes (--genome-size-files)

Ganon will automatically download auxiliary files to define an approximate genome size for each entry in the taxonomic tree. For `--taxonomy ncbi` the [species_genome_size.txt.gz](https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/) is used. For `--taxonomy gtdb` the [\*_metadata.tar.gz](https://data.gtdb.ecogenomic.org/releases/latest/) files are used. Those files can be directly provided with the `--genome-size-files` argument.

Genome sizes of parent nodes are calculated as the average of the respective children nodes. Other nodes without direct assigned genome sizes will use the closest parent with a pre-calculated genome size. The genome sizes are stored in the [ganon database](#buildupdate).

### Retrieving info (--ncbi-sequence-info, --ncbi-file-info)

Further taxonomy and assembly linking information has to be collected to properly build the database. `--ncbi-sequence-info` and `--ncbi-file-info` allow customizations on this step.

When `--input-target sequence`, `--ncbi-sequence-info` argument allows the use of NCBI e-utils webservices (`eutils`) or downloads accession2taxid files to extract target information (options `nucl_gb` `nucl_wgs` `nucl_est` `nucl_gss` `pdb` `prot` `dead_nucl` `dead_wgs` `dead_prot`). By default, ganon uses `eutils` up-to 50000 input sequences, otherwise it downloads `nucl_gb` `nucl_wgs` from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/. Previously downloaded files can be directly provided with this argument.

When `--input-target file`, `--ncbi-file-info` uses `assembly_summary.txt` from https://ftp.ncbi.nlm.nih.gov/genomes/ to extract target information (options `refseq` `genbank` `refseq_historical` `genbank_historical`. Previously downloaded files can be directly provided with this argument.

If you are using outdated, removed or inactive assembly or sequence files and accessions from NCBI, make sure to include `dead_nucl` `dead_wgs` for `--ncbi-sequence-info` or `refseq_historical` `genbank_historical` for `--ncbi-file-info`. `eutils` option does not work with outdated accessions.
