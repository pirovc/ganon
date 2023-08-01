# Databases

ganon automates the download, update and build of databases based on NCBI RefSeq and GenBank genomes repositories wtih `ganon build` and `update` commands, for example:

```bash
ganon build -g archaea bacteria -d arc_bac -c -t 30
```

will download archaeal and bacterial complete genomes from RefSeq and build a database with 30 threads. Some day later, the database can be updated to include newest genomes with:


```bash
ganon update -d arc_bac -t 30
```

Additionally, [custom databases](../custom_databases) can be built with customized files and identifiers with the `ganon build-custom` command.

!!! info
    We DO NOT provide pre-built indices for download. ganon can build databases very efficiently. This way, you will always have up-to-date reference sequences and get most out of your data.


!!! warning
    Databases are one of the most important elements when analyzing metagenomics data and should be chosen carefully based on the study data

## RefSeq and GenBank

NCBI RefSeq and GenBank repositories are common resources to obtain reference sequences to analyze metagenomics data. They are mainly divided into domains/organism groups (e.g. archaea, bacteria, fungi, ...) but can be further filtered in many ways. The choice of those filters can drastically change the outcome of results.


### Commonly used sub-sets

| RefSeq | # assemblies | Size (GB) * |  |
|---|---|---|---|
| Complete | 295219 | 160 - 500 | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 48 --db-prefix abfv_rs`</details> |
| One assembly per species | 52779 | 40 - 98 | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 48 --genome-updater "-A 'species:1'" --db-prefix abfv_rs_t1s`</details> |
| Complete genomes (higher quality) | 44121 | 19 - 64 | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 48 --complete-genomes --db-prefix abfv_rs_cg`</details> |
| One assembly per species of complete genomes | 19713 | 8 - 27 | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 48 --complete-genomes "-A 'species:1'" --db-prefix abfv_rs_cg_t1s`</details> |
| One representative assembly per species | 18073 | 21 - 35 | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea bacteria fungi viral --threads 48 --representative-genomes --db-prefix abfv_rs_rg`</details> |

| GenBank | # assemblies | Size (GB) * |  |
|---|---|---|---|
| Complete | 1595845 | - | <details><summary>cmd</summary>`ganon build --source genbank --organism-group archaea bacteria fungi viral --threads 48 --db-prefix abfv_gb`</details> |
| One assembly per species | 99505 | 91 - 420 | <details><summary>cmd</summary>`ganon build --source genbank --organism-group archaea bacteria fungi viral --threads 48 --genome-updater "-A 'species:1'" --db-prefix abfv_gb_t1s`</details> |
| Complete genomes (higher quality) | 92917 | 24 - 132 | <details><summary>cmd</summary>`ganon build --source genbank --organism-group archaea bacteria fungi viral --threads 48 --complete-genomes --db-prefix abfv_gb_cg`</details> |
| One assembly per species of complete genomes | 34497 | 10 - 34 | <details><summary>cmd</summary>`ganon build --source genbank --organism-group archaea bacteria fungi viral --threads 48 --complete-genomes "-A 'species:1'" --db-prefix abfv_gb_cg_t1s`</details> |

\*  Size (GB) is the final size of the database and the approximate amount of RAM necessary to build it (calculated with default parameters). The two values represent databases built with and without the `--hibf` parameter, respectively. The trade-offs between those two modes are explained [here](#hibf).

!!! warning
	The `# assemblies` were obtained on 2023-03-14 accounting for archaea, bacteria, fungi and viral groups only. By the time you are reading this, those numbers certainly grew a bit.

- As a rule of thumb, the more the better, so choose the most comprehensive sub-set as possible given your computational resources
- Databases can have a fixed size/RAM usage with the `--filter-size` parameter. Beware that smaller filters will increase the false positive rates when classifying. Other approaches [can reduce the size/RAM requirements with some trade-offs](#reducing-database-size).
- Further examples of commonly used database can be found [here](../custom_databases/#examples).
- Alternatively, you can build one database for each organism group separately and use them in `ganon classify` in any order or even stack them hierarchically. This way combination of multiple databases are possible, extending use cases.

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

| GTDB R214 | # assemblies | Size (GB) |  |
|---|---|---|---|
| Complete | 402709 | | <details><summary>cmd</summary>`ganon build --source refseq genbank --organism-group archaea bacteria --threads 48 --taxonomy gtdb --db-prefix ab_gtdb`</details> |
| One assembly for each species | 85205 | | <details><summary>cmd</summary>`ganon build --source refseq genbank --organism-group archaea bacteria --threads 48 --taxonomy gtdb --top 1 --db-prefix ab_gtdb_t1s`</details> |

Filtering by taxonomic entries also work with GTDB, for example:

```bash
ganon build --db-prefix fuso_gtdb --taxid "f__Fusobacteriaceae" --source refseq genbank --taxonomy gtdb --threads 12
```

!!! info
    GTDB covers only bacteria and archaea groups and has assemblies from both RefSeq and GenBank.

## Update (ganon update)

Default ganon databases generated with the `ganon build` can be updated with `ganon update`. This procedure will download new files and re-generate the ganon database with the updated entries.

For example, a database generated with the following command:

```bash
ganon build --db-prefix arc_cg_rs --source refseq --organism-group archaea --complete-genomes --threads 12
```

will contain all archaeal complete genomes from NCBI RefSeq at the time of running. Some days later, the database can be updated, fetching only new sequences added to the NCBI repository with the command:

```bash
ganon update --db-prefix arc_cg_rs --threads 12
```

!!! tip
    To not overwrite the current database and create a new one with the updated files, use the `--output-db-prefix` parameter.

## Reproducibility

If you use ganon with default databases and want to re-generate it later or keep track of the content for reproducibility purposes, you can save the `assembly_summary.txt` file located inside the `{output_prefix}_files/` directory. To re-download the exact same snapshot of files used, one could use [genome_updater](https://github.com/pirovc/genome_updater), for example:

```bash
genome_updater.sh -e assembly_summary.txt -f "genomic.fna.gz" -o recovered_files -m -t 12 
```

## Reducing database size

### False positive

A higher `--max-fp` value will generate a smaller database but with a higher number of false positive matches on classification. [More details](../custom_databases/#false-positive-and-size-max-fp-filter-size). Values between `0.001` (0.1%) and `0.3` (30%) are generally used. 

!!! hint
    When using higher `--max-fp` values, more false positive results may be generated. This can be filtered with the `--fpr-query` parameter in `ganon classify` 

### Fixed size

A fixed size for the database filter can be defined with `--filter-size`. The smaller the filter size, the higher the false positive chances on classification. When using a fixed filter size, ganon will report the max. and avg. false positive rate at the end of the build. [More details](../custom_databases/#false-positive-and-size-max-fp-filter-size).

### HIBF

The Hierarchical Interleaved Bloom Filter (HIBF) is an improvement over the default Interleaved Bloom Filter (IBF) and generates *smaller* databases with *faster* query times ([article](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02971-4)). However, the HIBF takes longer to build and has less flexibility regarding size and further options in ganon. The HIBF can be generated in `ganon build` and `ganon build-custom` with the `--hibf` parameter.

Due to differences between the default IBF used in ganon and the HIBF, it is recommended to lower the false positive when using the HIBF. A recommended value for high sensitivity is 1% (`--hibf --max-fp 0.001`).

!!! hint
    - For larger reference sets, huge amount of reads to query -> HIBF
    - For quick build and more flexibility -> IBF (default)

!!! warning
    [raptor (v3.0.0)](https://github.com/seqan/raptor/releases/tag/raptor-v3.0.0) has to be installed to build databases with `--hibf`

### Top assemblies

RefSeq and GenBank are highly biased toward some few organisms. This means that some species are highly represented in number of assemblies compared to others. This can not only bias analysis but also brings redundancy to the database. Choosing a certain number of top assemblies can mitigate those issues. Database sizes can also be drastically reduced without this redundancy, but "strain-level" analysis are then not possible. We recommend using top assemblies for larger and comprehensive reference sets (like the ones listed [above](#refseq-and-genbank)) and use the full set of assemblies for specific clade analysis.

!!! Example
    - `ganon build --top 1` will select one assembly for each taxonomic leaf (NCBI taxonomy still has strain, sub-species, ...)
    - `ganon build --genome-updater "-A 'species:1'"` will select one assembly for each species
    - `ganon build --genome-updater "-A 'genus:3'"` will select three assemblies for each genus

### Mode

`--mode` offers 5 different categories to build a database controlling the trade-off between size and classification speed.

- `avg`: Balanced mode
- `smaller` or `smallest`: create smaller databases with slower classification speed
- `fast` or `fastest`: create bigger databases with faster classification speed

!!! Warning
    If `--filter-size` is used, `smaller` and `smallest` refers to the false positive and not to the database size (which is fixed). 

### Split databases

Ganon allows classification with multiple databases in one level or in an hierarchy ([More details](../classification/#multiple-and-hierarchical-classification)). This means that databases can be built separately and used in any combination as desired. There are usually some benefits of doing so:

- Smaller databases when building by organism group, for example: one for bacteria, another for viruses, ... since average genome sizes are quite different.
- Easier to maintain and update.
- Extend use cases and avoid misclassification due to contaminated databases.
- Use databases as quality control, for example: remove reads matching one database of host or vectors (check out `ganon report --skip-hierarchy`).

### k-mer and window size

Define how much unique information is stored in the database. [More details](../custom_databases/#minimizers-window-size-kmer-size)

- The smaller the `--kmer-size`, the less unique they will be, reducing database size but also sensitivity in classification. 
- The bigger the `--window-size`, the less information needs to be stored resulting in smaller databases but with decrease classification accuracy.

### Example

Besides the benefits of using [HIBF](#hibf) and specific sub-sets of big repositories shown on the [default databases table](#commonly-used-sub-sets), examples of other reduction strategies (without `--hibf`) can be seen below:

*RefSeq archaeal complete genomes from 20230505*

| Strategy | Size (MB) | Smaller |  Trade-off  | |
|---|---|---|---|---|
| `default` | 318 | - | - | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea --threads 12 --complete-genomes --db-prefix arc_rs_cg`</details> |
| `--mode smallest` | 301 | 5% | Slower classification | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea --threads 12 --complete-genomes --mode smallest --db-prefix arc_rs_cg_smallest`</details> |
| `--filter-size 256` | 256 | 19% | Higher false positive on classification | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea --threads 12 --complete-genomes --filter-size 256 --db-prefix arc_rs_cg_fs256`</details> |
| `--window-size 35` | 249 | 21% | Less sensitive classification | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea --threads 12 --complete-genomes --window-size 35 --db-prefix arc_rs_cg_ws35`</details> |
| `--max-fp 0.2` | 190 | 40% | Higher false positive on classification | <details><summary>cmd</summary>`ganon build --source refseq --organism-group archaea --threads 12 --complete-genomes --max-fp 0.2 --db-prefix arc_rs_cg_fp0.2`</details> |

!!! note
    This is an illustrative example and the reduction proportions for different configuration may be quite different
