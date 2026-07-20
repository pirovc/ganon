# Databases

`ganon build` downloads and builds databases based on NCBI RefSeq and/or GenBank genomes repositories with [GTDB](https://gtdb.ecogenomic.org){target="_blank"} support. For example:

```bash
# All archaeal and bacterial genomes from current RefSeq
ganon build -b refseq -g archaea bacteria -d rs_arc_bac -c -t 30
```

```bash
# Some days later, the database can be synced to the latest NCBI version with
ganon update -d arc_bac -t 30
```

- Pre-built databases are not provided. `ganon build` downloads and build databases easily and efficiently. That way you get the latest and most diverse data available.
    - The [command generator](#simple-ganon-build-command-generator) can help you start with the build command.
- To build databases based on local/non-standard files, check the `ganon build-custom` command and [documentation](custom_databases.md)
    - Examples of commonly used (not standard) databases can be found [here](custom_databases.md#examples).

## Simple `ganon build` command generator

<iframe src="https://pirovc.github.io/ganon/ganon_build_generator.html" width="730" height="530" frameborder="0"></iframe>

!!! note
    [More filters](#filters) and [parameters](params.md) are available for `ganon build`

!!! tip
    To have more flexibility and extend use cases you can build database separetly for each organism group. In `ganon classify` you can [combine multiple databases in one run or stack them hierarchically](classification.md#multiple-and-hierarchical-classification).

## Commonly used sub-sets

The table below lists the resources and time needed to build commonly used sub-sets. By the time you read this, these numbers will have increased slightly. As a rule of thumb, the more the better, so choose the most comprehensive sub-set as possible given your computational resources.

- To build smaller databases with less memory, you can apply some [filters](#filters) or [tweak parameters](#reducing-database-size). Note that there will be trade-offs in every reduction.

| RefSeq ¹ | #assemblies | #species | Size ² | Time ² | `ganon build` |
|:--------:|:-----------:|:--------:|:------:|:------:|:-------------:|
| Archaea, Bacteria <br> [**complete genomes**] | 62944 | 15693 | 63 | 47m | <details><summary></summary>`ganon build --source refseq --organism-group archaea bacteria --threads 48 --complete-genomes --db-prefix rs_arc_bac_cg`</details> |
| Archaea, Bacteria <br> [**reference genomes**] | 23404 | 23401 | 77 | 28m | <details><summary></summary>`ganon build --source refseq --organism-group archaea bacteria --threads 48 --reference-genomes --db-prefix rs_arc_bac_rg`</details> |
| Archaea, Bacteria <br> [**complete + reference**] | 79268 | 30912 | 194 | 58m | <details><summary></summary>`ganon build --source refseq --organism-group archaea bacteria --threads 48 --db-prefix rs_arc_bac_cgrg --verbose --genome-updater "-F $(printf "'%s'" '$5 == "reference genome" || $12 == "Complete Genome"')"`</details> |
| Archaea, Bacteria | 509046 | 76966 | 312 | 10h | <details><summary></summary>`ganon build --source refseq --organism-group archaea bacteria --threads 48 --db-prefix rs_arc_bac`</details> |
| Fungi | 674 | 668 | 18 | 5m | <details><summary></summary>`ganon build --source refseq --organism-group fungi --threads 48 --db-prefix rs_fungi`</details> |
| Human | 2 | 1 | 2.3 | 4m | <details><summary></summary>`ganon build --source refseq --organism-group human --threads 48 --db-prefix rs_human`</details> |
| Plant | 202 | 202 | 79 | 29m | <details><summary></summary>`ganon build --source refseq --organism-group plant --threads 48 --db-prefix rs_plant`</details> |
| Protozoa | 129 | 126 | 3 | 81s| <details><summary></summary>`ganon build --source refseq --organism-group protozoa --threads 48 --db-prefix rs_protozoa`</details> |
| Viral | 15089 | 14082 | 0.41 | 32m | <details><summary></summary>`ganon build --source refseq --organism-group viral --threads 48 --db-prefix rs_viral`</details> |

| Others ¹ | #sequences | #species | Size ² | Time ² | `ganon build-custom` |
|:--------:|:----------:|:--------:|:------:|:------:|:--------------------:|
| Plasmid | 135944 | 6708 | 3 | 13m | [build-custom](custom_databases.md#plasmid-plastid-and-mitochondrion-from-refseq) |
| UniVec_Core | 3155 | 1 | 0.0004 | 13s | [build-custom](custom_databases.md#univec-univec_core) |

| GTDB              | #assemblies | #species | Size ² | Time ² | `ganon build` |
|:-----------------:|:-----------:|:--------:|:------:|:------:|:-------------:|
| R232 (2026-04-15) | 900653 | 199913 | 652 | 17h | <details><summary></summary>`ganon build --source refseq genbank --organism-group archaea bacteria --threads 48 --taxonomy gtdb --db-prefix ab_gtdb`</details> |

¹ *data from 2026-06-21*

² *"Size" (in GB) is the final ganon database size. "Time" accounts for wall time for the build process after downloading files. The memory required for the build is approximate 1.5x the database size. Your time may vary based on internet, I/O, memory and CPU speed. 64 threads were used with an AMD EPYC 9454 48-Core Processor, using ganon v2.4.2.*




<details>
  <summary>Older data for comparison</summary>

```txt

|          RefSeq (2025-11-01) *        | # assemblies | # species | Size (GB) |
|---------------------------------------|--------------|-----------|-----------|
| All genomes                           | 468399       | 83798     | 299       |
| Complete genomes (CG)                 | 68639        | 27598     | 55        |
| Reference genomes (RG)                | 22862        | 22861     | 89        | 
| CG + RG                               | 85036        | 42955     | 120       | 

* archaea, bacteria, fungi and viral
```

```txt

|          RefSeq (2024-04-20) *        | # assemblies | # species | Size (GB) |
|---------------------------------------|--------------|-----------|-----------|
| All genomes                           | 366941       | 64616     | 215       |
| Complete genomes (CG)                 | 55114        | 24238     | 42        |
| Reference genomes (RG)                | 19890        | 19888     | 77        | 
| CG + RG                               | 69600        | 37864     | 100       | 

* archaea, bacteria, fungi and viral
```

```txt
|          RefSeq (2023-03-14) *        | # assemblies | # species | Size (GB) |
|---------------------------------------|--------------|-----------|-----------|
| All genomes                           | 295219       | 52781     | 160       |
| All genomes - 1 assembly/species      | 52781        | 52781     | 128       |
| Complete genomes                      | 44121        | 19715     | 35        |
| Complete genomes - 1 assembly/species | 19715        | 19715     | 29        |
| Reference genomes                     | 18073        | 18073     | 69        |

* archaea, bacteria, fungi and viral
```

```txt
|          GenBank (2023-03-14) *       | # assemblies | # species | Size (GB) |
|---------------------------------------|--------------|-----------|-----------|
| All genomes - 1 assembly/species      | 99505        | 99505     | 300       |
| Complete genomes                      | 92917        | 34815     | 42        |
| Complete genomes - 1 assembly/species | 34815        | 34815     | 34        |

* archaea, bacteria, fungi and viral
```

```txt
|                GTDB              | # assemblies | # species | Size (GB) |
|----------------------------------|--------------|-----------|-----------|
| R226                             | 731982       | 143396    | 501       |
| R220                             | 596859       | 113104    | 338       |
| R214                             | 402709       | 85205     | 260       |
```

</details>
<br>

!!! tip
    RefSeq is preferred mainly due to its superior sequence curation and quality. In the experiments [published in the ganon2 article](https://dx.doi.org/10.1093/nargab/lqaf094){target="_blank"}, the more reference genomes used, the better the results. However, this requires significant computational resources. Combining complete and reference genomes (CG+RG) strikes a good balance, providing good results with a smaller memory footprint and faster classification. Beware that the choice of the database will drastically affect the outcome of the analysis.


## Filters

### Specific taxa

It is also possible to generate databases for one or more taxonomic branches with `-a/--taxid`, for example:

```bash
ganon build --source refseq --taxid 562 317 --threads 48 --db-prefix coli_syringae
```

will download and build a database for all *Escherichia coli* (taxid:562) and *Pseudomonas syringae* (taxid:317) assemblies from RefSeq.

This is also possible with `--taxonomy gtdb`, for example:

```bash
ganon build --db-prefix fuso_gtdb --taxid "f__Fusobacteriaceae" --source refseq genbank --taxonomy gtdb --threads 12
```

### Top genomes/taxa

Select a specific number of genomes/assemblies for each taxa in the database. For example:

- `--top 3` will select three assemblies for each taxonomic leaf
- `--genome-updater "-A 'species:1'"` will select one assembly for each species node

[More infos](#top-assemblies) about top assemblies.

### Refined filters

ganon uses [genome_updater](https://github.com/pirovc/genome_updater){target="_blank"} to manage downloads and further specific options and filters can be provided with the paramer `-u/--genome-updater`, for example:

```bash
ganon build -g bacteria -t 48 -d bac_refseq --genome-updater "-A 'genus:3' -E 20230101"
```

will download top 3 archaeal assemblies for each genus with date before 2023-01-01. For more information about genome_updater parameters, please check the [repository](https://github.com/pirovc/genome_updater){target="_blank"}.

## Update (ganon update)

Default ganon databases generated with the `ganon build` can be updated with `ganon update`. This procedure will download new files and re-generate the ganon database adding new entires and removing outdated ones. This will keep the choosen database selection in sync with the latest available data.

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

If you use ganon with default databases and want to re-generate it later or keep track of the content for reproducibility purposes, you can save the `assembly_summary.txt` file located inside the `{output_prefix}_files/` directory. To re-download the exact same snapshot of files used, one could use [genome_updater](https://github.com/pirovc/genome_updater){target="_blank"}, for example:

```bash
genome_updater.sh -e assembly_summary.txt -f "genomic.fna.gz" -o recovered_files -m -t 12 
```

## Reducing database size

### False positive

A higher `--max-fp` value will generate a smaller database but with a higher number of false positive matches on classification. [More details](custom_databases.md#false-positive-max-fp). Values between `0.001` (0.1%) and `0.3` (30%) are generally used. 

!!! hint
    When using higher `--max-fp` values, more false positive results may be generated. This can be filtered with the `--fpr-query` parameter in `ganon classify` 


### k-mer and window size

Define how much unique information is stored in the database. [More details](custom_databases.md#minimizers-window-size-kmer-size)

- The smaller the `--kmer-size`, the less unique they will be, reducing database size but also sensitivity in classification. 
- The bigger the `--window-size`, the less information needs to be stored resulting in smaller databases but with decrease classification accuracy.


### Top assemblies

RefSeq and GenBank are highly biased toward some few organisms. This means that some species are highly represented in number of assemblies compared to others. This can bias analysis towards those organisms. Choosing a certain number of top assemblies can mitigate those issues. Database sizes can also be drastically reduced without this redundancy, but "strain-level" analysis are then not possible. We recommend using top assemblies for larger and comprehensive reference sets (like the ones listed [above](#commonly-used-sub-sets)) and use the full set of assemblies for specific clade analysis.

!!! Example
    - `ganon build --top 1` will select one assembly for each taxonomic leaf (NCBI taxonomy still has strain, sub-species, ...)
    - `ganon build --genome-updater "-A 'species:1'"` will select one assembly for each species
    - `ganon build --genome-updater "-A 'genus:3'"` will select three assemblies for each genus

### Database level

With the `--level` parameter one can define the final taxonomic level of the database. It can be a taxonomic rank ['species', 'genus', ...], 'leaves' for taxonomic leaves or 'assembly' for a assembly/strain based analysis. The default value in `ganon build` is `species` but if you don't need species resolution you can set to a less specific rank (e.g. `genus`). That will generate smaller databases.

### Split databases

Ganon allows classification with multiple databases in one level or in an hierarchy ([More details](classification.md#multiple-and-hierarchical-classification)). This means that databases can be built separately and used in any combination as desired. There are usually some benefits of doing so:

- Smaller databases when building by organism group, for example: one for bacteria, another for viruses, etc.
- Easier to maintain and update.
- Extend use cases and avoid misclassification due to contaminated databases.
- Use databases as quality control, for example: remove reads matching one database of host or vectors (check out `ganon report --skip-hierarchy`).
