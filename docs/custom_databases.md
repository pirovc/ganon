# Custom databases

## Default NCBI assembly or sequence accession

Besides the automated download and build (`ganon build`) ganon provides a highly customizable build procedure (`ganon build-custom`) to create databases from local files.

To use custom sequences, just provide them with `--input`. ganon will try to retrieve all necessary information necessary to build a database.

!!! tip
    ganon expects assembly accessions in the filename like `GCA_002211645.1_ASM221164v1_genomic.fna.gz`. When using `--input-target sequence` filenames are not important but sequence headers should look like `>CP022124.1 Fusobacterium nu...`. More information about building by file or sequence can be found [here](#target-file-or-sequence---input-target).

## Non-standard/custom accessions

It is also possible to use **non-standard accessions and headers** to build databases with `--input-file`. This file should contain the following fields (tab-separated): file, [target, node, specialization, specialization_name].

<details>
  <summary>Examples of --input-file</summary>
  <br>

With --input-target file, where my_target_1 and my_target_2 are just names to assign sequences from (unique) sequence files:

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

Below a list of few example of custom databases from commonly used repositories:

### HumGut

TODO description

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

TODO description

```bash
# Download sequence files
wget -A genomic.fna.gz -m -nd --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plasmid/"
wget -A genomic.fna.gz -m -nd --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plastid/"
wget -A genomic.fna.gz -m -nd --quiet --show-progress "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/mitochondrion/"

# Build ganon database
ganon build-custom --input plasmid.* plastid.* mitochondrion.* --db-prefix ppm --input-target sequence --level leaves --threads 32 
```

### UniVec, UniVec_core

TODO description

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

!!! warning
    All UniVec entries in the examples are categorized as [Artificial Sequence](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=81077) (NCBI txid:81077). Some are completely artificial but others may be derived from real biological sources. More information in this [link](https://ftp.ncbi.nlm.nih.gov/pub/UniVec/README.vector.origins).

### MGnify genome collections (MAGs)

TODO description

Note that the database will be build with GTDB taxonomy.

```bash
# Download metadata
wget "https://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-oral/v1.0/genomes-all_metadata.tsv"

# Download sequence files with 12 threads
tail -n+2 genomes-all_metadata.tsv | cut -f 1,20 | xargs -P 12 -n2 sh -c 'curl --silent ${1}| gzip -d | sed -e "1,/##FASTA/ d" | gzip > ${0}.fna.gz'

# Generate ganon input file
tail -n+2 genomes-all_metadata.tsv | cut -f 1,15  | tr ';' '\t' | awk -F"\t" '{tax="1";for(i=NF;i>1;i--){if(length($i)>3){tax=$i;break;}};print $1".fna.gz\t"$1"\t"tax}' > ganon_input_file.tsv

ganon build-custom --input-file ganon_input_file.tsv --db-prefix mgnify_human_oral_v1 --taxonomy gtdb --level leaves --threads 32
```

### Pathogens NCBI

- https://ftp.ncbi.nlm.nih.gov/pathogen/
- https://www.ncbi.nlm.nih.gov/pathogens/

### Pathogen detection FDA-ARGOS

- https://www.ncbi.nlm.nih.gov/bioproject/231221
- https://pubmed.ncbi.nlm.nih.gov/31346170/

### Reference Viral DataBase (RVDB)

- https://journals.asm.org/doi/10.1128/mSphereDirect.00069-18
- https://rvdb.dbi.udel.edu/

### EuPathDB

- https://veupathdb.org/veupathdb/app/

### BLAST databases (nt, env_nt, ...)

- https://ftp.ncbi.nlm.nih.gov/blast/db/

### From genome_updater