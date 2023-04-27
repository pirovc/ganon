# Tutorials

## Custom build

## Choosing and explaining parameters

The most important parameters and trade-offs are:

- `ganon build` `--window-size --kmer-size`: the *window* value should always be the same or larger than the *kmer* value. The larger the difference between them, the smaller the database will be. However, some sensitivity/precision loss in classification is expected with small *kmer* and/or large *window*. Larger *kmer* values (e.g. `31`) will improve classification, specially read binning, at a cost of way bigger databases.
- `ganon classify` `--rel-cutoff`: this value defines the threshold for matches between reads and database. Higher `--rel-cutoff` values will improve precision and decrease sensitivity with expected less unique matches but an increase in overall matches. For taxonomic profiling, a higher value between `0.4` and `0.8` may provide better results. For read binning, lower values between `0.2` and `0.4` are recommended.
- `ganon classify` `--rel-filter`: further filter top matches after cutoff is applied. Usually set between `0` and `0.2`.
- `ganon classify` `--reassign`: runs an EM-algorithm to reassign reads that received multiple matches. It provides a unique match for each read at the level the database was built (e.g. assembly or species). Mostly useful for read binning, with little overall impact on taxonomic profiling. Can be used independently with `ganon reassign`.
- `ganon report` `--report-type`: reports either taxonomic, sequence or matches abundances. Use `corr` or `abundance` for taxonomic profiling, `reads` or `dist` for sequence profiling and `matches` to report a summary of all matches.
- `ganon report` `--min-count`: cutoff to discard underrepresented taxa. Useful to remove the common long tail of spurious matches and false positives when performing classification. Values between `0.0001` (0.01%) and `0.001` (0.1%) improved sensitivity and precision in our evaluations. The higher the value, the more precise the outcome, with a sensitivity loss. Alternatively `--top-percentile` can be used to keep a relative amount of taxa instead a hard cutoff.

The numeric values above are averages from several experiments with different sample types and database contents. They may not work as expected for your data. If you are not sure which values to use or see something unexpected, please open an [issue](https://github.com/pirovc/ganon/issues).

Below, a more in-depth explanation of important parameters:

### ganon build 

#### filter false positive and size (--max-fp, --filter-size)

ganon indices are based on bloom filters and can have false positive matches. This can be controlled with `--max-fp` parameter. The lower the `--max-fp`, the less chances of false positives, but the larger the database size will be. For example, with `--max-fp 0.01` the database will be build so any target (e.g. assembly, specificed with `--level`) will have 1 in a 100 change of reporting a false match (between minimizer k-mers).

Alternatively, one can set a specific size for the final index with `--filer-size`. When using this option, please observe the theoretic false positive of the index reported at the end of the building process.

#### minimizers (--window-size, --kmer-size)

in `ganon build`, when `--window-size` > `--kmer-size` minimizers are used. That means that for a every window, a single k-mer will be selected. It produces smaller database files and requires substantially less memory overall. It may increase building times but will have a huge benefit for classification times. Sensitivity and precision can be reduced by small margins. If `--window-size` = `--kmer-size`, all k-mers are going to be used to build the database.

### ganon classify

#### reads (--single-reads, --paired-reads)

ganon accepts single-end and paired-end reads. Both types can be use at the same time. In paired-end mode, reads are always reported with the header of the first pair. Paired-end reads are classified in a standard forward-reverse orientation. The max. read length accepted is 65535 (accounting both reads in paired mode).

#### cutoff and filter (--rel-cutoff, --rel-filter)

ganon has two parameters to control a match between reads and references: `--rel-cutoff` and `--rel-filter`.

Every read can be classified against none, one or more references. What will be reported is the remaining matches after `cutoff` and `filter` thresholds are applied, based on the number of shared minimizers (or k-mers) between sequences.

The `cutoff` is the first. It should be set as a minimal value to consider a match between a read and a reference. Next the `filter` is applied to the remaining matches. `filter` thresholds are relative to the best scoring match and control how far from the best match further matches are allowed. `cutoff` can be interpreted as the lower bound to discard spurious matches and `filter` as the fine tuning to control what to keep.

For example, using `--kmer-size 19` (and `--window-size 19` to simplify the example), a certain read (100bp) has the following matches with the 5 references (`ref1..5`), sorted by shared k-mers:

| reference | shared k-mers |
|-----------|---------------|
| ref1      | 82            |
| ref2      | 68            |
| ref3      | 44            |
| ref4      | 25            |
| ref5      | 20            |

this read can have at most 82 shared k-mers (`100-19+1=82`). With `--rel-cutoff 0.25`, the following matches will be discarded:

| reference | shared k-mers | --rel-cutoff 0.25 |
|-----------|---------------|-------------------|
| ref1      | 82            |                   |
| ref2      | 68            |                   |
| ref3      | 44            |                   |
| ref4      | 25            |                   |
| ~~ref5~~  | ~~20~~        | X                 |

since the `--rel-cutoff` threshold is `82 * 0.25 = 21` (ceiling is applied). Further, with `--rel-filter 0.3`, the following matches will be discarded:

| reference | shared k-mers | --rel-cutoff 0.25 | --rel-filter 0.3 |
|-----------|---------------|-------------------|------------------|
| ref1      | 82            |                   |                  |
| ref2      | 68            |                   |                  |
| ~~ref3~~  | ~~44~~        |                   | X                |
| ~~ref4~~  | ~~25~~        |                   | X                |
| ~~ref5~~  | ~~20~~        | X                 |                  |


since best match is 82, the filter parameter is removing any match below `0.3 * 82 = 57` (ceiling is applied) shared k-mers. `ref1` and `ref2` are reported as matches.

For databases built with `--window-size`, the relative values are not based on the maximum number of possible shared k-mers but on the actual number of unique minimizers extracted from the read.

A different `cutoff` can be set for every database in a multiple or hierarchical database classification. A different `filter` can be set for every level of a hierarchical database classification.

Note that reads that remain with only one reference match (after `cutoff` and `filter` are applied) are considered a unique match.

### ganon build-custom

#### Target file or sequence (--input-target) 

Customized builds can be done either by file or sequence. `--input-target file` will consider every file provided with `--input` a single unit. `--input-target sequence` will use every sequence as a unit.

`--input-target file` is the default behavior and most efficient way to build databases. `--input-target sequence` should only be used when the input sequences are stored in a single file or when classification at sequence level is desired.

#### Build level (--level)

The `--level` parameter defines the max. depth of the database for classification. This parameter is relevant because the `--max-fp` is going to be guaranteed at the `--level` chosen. By default, the level will be the same as `--input-target`, meaning that classification will be done either at file or sequence level.

Alternatively, `--level assembly` will link the file or sequence target information with assembly accessions retrieved from NCBI servers. `--level leaves` or `--level species` (or genus, family, ...) will link the targets with taxonomic information and prune the tree at the chosen level. `--level custom` will use specialization level define in the `--input-file`.

#### Genome sizes (--genome-size-files)

Ganon will automatically download auxiliary files to define an approximate genome size for each entry in the taxonomic tree. For `--taxonomy ncbi` the [species_genome_size.txt.gz](https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/) is used. For `--taxonomy gtdb` the [\*_metadata.tar.gz](https://data.gtdb.ecogenomic.org/releases/latest/) files are used. Those files can be directly provided with the `--genome-size-files` argument.

Genome sizes of parent nodes are calculated as the average of the respective children nodes. Other nodes without direct assigned genome sizes will use the closest parent with a pre-calculated genome size. The genome sizes are stored in the [ganon database](#buildupdate).

#### Retrieving info (--ncbi-sequence-info, --ncbi-file-info)

Further taxonomy and assembly linking information has to be collected to properly build the database. `--ncbi-sequence-info` and `--ncbi-file-info` allow customizations on this step.

When `--input-target sequence`, `--ncbi-sequence-info` argument allows the use of NCBI e-utils webservices (`eutils`) or downloads accession2taxid files to extract target information (options `nucl_gb` `nucl_wgs` `nucl_est` `nucl_gss` `pdb` `prot` `dead_nucl` `dead_wgs` `dead_prot`). By default, ganon uses `eutils` up-to 50000 input sequences, otherwise it downloads `nucl_gb` `nucl_wgs` from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/. Previously downloaded files can be directly provided with this argument.

When `--input-target file`, `--ncbi-file-info` uses `assembly_summary.txt` from https://ftp.ncbi.nlm.nih.gov/genomes/ to extract target information (options `refseq` `genbank` `refseq_historical` `genbank_historical`. Previously downloaded files can be directly provided with this argument.

If you are using outdated, removed or inactive assembly or sequence files and accessions from NCBI, make sure to include `dead_nucl` `dead_wgs` for `--ncbi-sequence-info` or `refseq_historical` `genbank_historical` for `--ncbi-file-info`. `eutils` option does not work with outdated accessions.

### ganon report

#### report type (--report-type)

Several reports are availble with `--report-type`: `reads`, `abundance`, `dist`, `corr`, `matches`:

`reads` reports **sequence abundances** which are the basic proportion of reads classified in the sample.

`abundance` will convert sequence abundance into **taxonomic abundances** by re-distributing read counts among leaf nodes and correcting by genome size. The re-distribution applies for reads classified with a LCA assignment and it is proportional to the number of unique matches of leaf nodes available in the ganon database (relative to the LCA node). Genome size is estimated based on [NCBI or GTDB auxiliary files](#genome-sizes---genome-size-files). Genome size correction is applied by rank based on default ranks only (superkingdom phylum class order family genus species assembly). Read counts in intermediate ranks will be corrected based on the closest parent default rank and re-assigned to its original rank.

`dist` is the same of `reads` with read count re-distribution

`corr` is the same of `reads` with correction by genome size

`matches` will report the total number of matches classified, either unique or shared. *This report will output the total number of matches instead the total number of reads reported in all other reports.*


## Examples

### Commonly used reference set for metagenomics analysis (archaea+bacteria+fungi+viral from NCBI RefSeq, complete genomes, one assembly per taxa)

```bash
ganon build --db-prefix abfv_rs_cg_t1a --organism-group archaea bacteria fungi viral --source refseq --taxonomy ncbi --complete-genomes --threads 12 --top 1
```

### GTDB database, one assembly per taxa
```bash
ganon build --db-prefix complete_gtdb --organism-group archaea bacteria --source refseq genbank --taxonomy gtdb --threads 12 --top 1
```

### Database based on specific taxonomic identifiers (203492 - Fusobacteriaceae)
```bash
# NCBI
ganon build --db-prefix fuso_ncbi --taxid "203492" --source refseq genbank --taxonomy ncbi --threads 12
# GTDB
ganon build --db-prefix fuso_gtdb --taxid "f__Fusobacteriaceae" --source refseq genbank --taxonomy gtdb --threads 12
```

### Customized database at assembly level 
```bash
ganon build-custom --db-prefix my_db --input my_big_fasta_file.fasta.gz --level assembly --threads 12
```

### Customized database at species level build based on files previously downloaded with [genome_updater](https://github.com/pirovc/genome_updater)
```bash
ganon build-custom --db-prefix custom_db_gu --input myfiles/2022-06-28_10-02-14/files/ --level species --ncbi-file-info outfolder/2022-06-28_10-02-14/assembly_summary.txt --threads 12
```

### Customized database with sequence as target (to classify reads at sequence level)
```bash
ganon build-custom --db-prefix seq_target --input myfiles/2022-06-28_10-02-14/files/ --ncbi-file-info outfolder/2022-06-28_10-02-14/assembly_summary.txt --input-target sequence --threads 12
```

