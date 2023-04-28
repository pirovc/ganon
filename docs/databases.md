# Databases

Databases are one of the most important elements when analyzing metagenomics data and should be chosen carefuly.

ganon automates the download, update and build of databases based on NCBI RefSeq and GenBank genomes repositories wtih `ganon build` and `update` commands, for example:

```bash
ganon build -g archaea bacteria -d arc_bac -c -t 30
```

will download archaeal and bacterial complete genomes from RefSeq and build a database with 30 threads. Some day later, the database can be updated to include newest genomes with:


```bash
ganon update -d arc_bac -t 30
```

Additionally, [custom databases](#custom-databases) can be built with customized files and identifiers with the `ganon build-custom` command.

!!! info
    We DO NOT provide pre-built indices for download. ganon can build databases very efficiently. This way, you will always have up-to-date reference sequences and get most out of your data.

## RefSeq/GenBank

NCBI RefSeq and GenBank repositories are common resources to obtain reference sequences to analyze metagenomics data. They are mainly divided into domains/organism groups (e.g. archaea, bacteria, fungi, ...) but can be further filtered in many ways. The choice of those filters can be confusing and drastically change the outcome of results. Listed below are some examples of commonly used sub-sets:

| RefSeq | # assemblies | Size (GB) |  |
|---|---|---|---|
| Complete | 295219 | xx-500G | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 48 --db-prefix abfv_rs`</details> |
| One assembly for each species | 52779 | xx-98G | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 48 --genome-updater "-A 'species:1'" --db-prefix abfv_rs_t1s`</details> |
| Complete genome assemblies (higher quality) | 44121 | xx-64G | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 48 --complete-genomes --db-prefix abfv_rs_cg`</details> |
| One assembly for each species of complete genomes | 19713 | xx-27G | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 48 --complete-genomes "-A 'species:1'" --db-prefix abfv_rs_cg_t1s`</details> |
| One representative assembly for each species | 18073 | xx-35G | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 48 --genome-updater "-c 'representative genome'" --db-prefix abfv_rs_rg`</details> |

| GenBank | # assemblies | Size (GB) |  |
|---|---|---|---|
| Complete | 1595845 | | <details><summary>cmd</summary>`ganon build --source genbank --organism-group archaea bacteria fungi viral --threads 48 --db-prefix abfv_gb`</details> |
| One assembly for each species | 99505 | xx-420G | <details><summary>cmd</summary>`ganon build --source genbank --organism-group archaea bacteria fungi viral --threads 48 --genome-updater "-A 'species:1'" --db-prefix abfv_gb_t1s`</details> |
| Complete genome assemblies (higher quality) | 92917 | | <details><summary>cmd</summary>`ganon build --source genbank --organism-group archaea bacteria fungi viral --threads 48 --complete-genomes --db-prefix abfv_gb_cg`</details> |
| One assembly for each species of complete genomes | 34497 | | <details><summary>cmd</summary>`ganon build --source genbank --organism-group archaea bacteria fungi viral --threads 48 --complete-genomes "-A 'species:1'" --db-prefix abfv_gb_cg_t1s`</details> |

!!! warning
	The `# assemblies` were calculated on 2023-03-14 accounting for archaea, bacteria, fungi and viral groups only. By the time you are reading this, those numbers certainly grew a bit.

!!! info
    You can build one database for each organism group separately and use them together in `ganon classify` in any order or even stack them hierarchically. This way, several combinations of databases are possible for many use cases.

- As a rule of thumb, the more the better, so choose the most comprehensive sub-set as possible given your computational resources
- The `Size (GB)` is the final size of the filter and the approximate amount of RAM necessary to build it (calculated with default parameters).
- Databases can have a fixed size/RAM usage with the `--filter-size` parameter. Beware that smaller filters will increase the false positive rates when classifying. Other approaches [can reduce the size/RAM requirements with some trade-offs](#reducing-database-size).

!!! tip
    assemblies usually represent strains and can be used to do "strain-level classification"

### Specific organisms or taxonomic groups

It is also possible to generate databases for specific organisms or taxonomic branches with `-a/--taxid`, for example:

```bash
ganon build --source refseq --taxid 562 317 --threads 48 --db-prefix coli_syringae
```

will download and build a database for all *Escherichia coli* (taxid:562) and *Pseudomonas syringae* (taxid:317) assemblies from RefSeq.

### More filter options

ganon uses [genome_updater](https://github.com/pirovc/genome_updater) to manage downloads and further specific options and filters can be provided with the paramer `-u/--genome-updater`, for example:

```bash
ganon build -g bacteria -t 48 -d bac_refseq --genome-updater "-A 'genus:3' -E 20230101"
```

will download top 3 archaeal assemblies for each genus with date before 2023-01-01. For more information about genome_updater parameters, please check the [repository](https://github.com/pirovc/genome_updater).

## GTDB

By default, ganon will use the NCBI Taxonomy to build the database. However, [GTDB](gtdb.ecogenomic.org/) is fully supported and can be used with the parameter `--taxonomy gtdb`.

| GTDB R214 | # assemblies | Space/Mem |  |
|---|---|---|---|
| Complete | 402709 | | <details><summary>cmd</summary>`ganon build --source refseq genbank --organism-group archaea bacteria --threads 48 --taxonomy gtdb --db-prefix ab_gtdb`</details> |
| One assembly for each species | 85205 | | <details><summary>cmd</summary>`ganon build --source refseq genbank --organism-group archaea bacteria --threads 48 --taxonomy gtdb --top 1 --db-prefix ab_gtdb_t1s`</details> |

Filtering by taxonomic entries also work with GTDB, for example:

```bash
ganon build --db-prefix fuso_gtdb --taxid "f__Fusobacteriaceae" --source refseq genbank --taxonomy gtdb --threads 12
```

!!! warning
    GTDB contain only bacteria and archea groups

## Update

## Reducing database size

- Top organisms:
- Split groups
- IBF and HIBF
- k-mer and window size
- false positive rate

## Custom databases

Besides the automated download and build (`ganon build`) ganon provides a highly customizable build procedure (`ganon build-custom`) to create databases. 

To use custom sequences, just provide them with `--input`. ganon will try to retrieve all necessary information necessary to build a database.

*ganon expects assembly accessions if building by file (e.g. filename should be similar as `GCA_002211645.1_ASM221164v1_genomic.fna.gz`) or accession version if building by sequence (e.g. headers should look like `>CP022124.1 Fusobacterium nu...`).* More information about building by file or sequence can be found [here](#target-file-or-sequence---input-target).

It is also possible to use non-standard accessions and headers to build databases with `--input-file`. This file should contain the following fields (tab-separated): file, [target, node, specialization, specialization_name].

<details>
  <summary>Examples of --input-file</summary>

Using `--input-target sequence`:

```
sequences.fasta HEADER1
sequences.fasta HEADER2
sequences.fasta HEADER3
others.fasta HEADER4
others.fasta HEADER5
```

or using `--input-target file`:

```
sequences.fasta FILE_A
others.fasta FILE_B
```

Nodes can be provided to link the data with taxonomy. For example (using `--taxonomy ncbi`):

```
sequences.fasta HEADER1 562
sequences.fasta HEADER2 562
sequences.fasta HEADER3 562
others.fasta HEADER4  623
others.fasta HEADER5  623
```

Specialization can be used to create a additional classification level after the taxonomic leaves. For example (using `--level custom`):

```
sequences.fasta HEADER1 562 ID443 Escherichia coli TW10119
sequences.fasta HEADER2 562 ID297 Escherichia coli PCN079
sequences.fasta HEADER3 562 ID8873  Escherichia coli P0301867.7
others.fasta HEADER4  623 ID2241  Shigella flexneri 1a
others.fasta HEADER5  623 ID4422  Shigella flexneri 1b
```
</details>
<br>

Below a list of few example of custom built databases from commonly used repositories:

### HumGut
### Plasmid, Plastid
### UniVec


