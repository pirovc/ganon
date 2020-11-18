# ganon [![Build Status](https://travis-ci.org/pirovc/ganon.svg?branch=master)](https://travis-ci.org/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/downloads.svg)](https://anaconda.org/bioconda/ganon) [![Anaconda-Server Badge](https://anaconda.org/bioconda/ganon/badges/platforms.svg)](https://anaconda.org/bioconda/ganon)

a k-mer based read classification tool which uses Interleaved Bloom Filters in conjunction with a taxonomic clustering and a k-mer counting-filtering scheme. 

> **ganon: precise metagenomics classification against large and up-to-date sets of reference sequences**
> Vitor C. Piro, Temesgen H. Dadi, Enrico Seiler, Knut Reinert, and Bernhard Y. Renard
> Bioinformatics 2020 doi: [10.1101/406017](https://dx.doi.org/10.1093/bioinformatics/btaa458)

## Installation

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html)

```shh
conda install -c bioconda -c conda-forge ganon
```

* There are possible performance benefits compiling ganon from source rather than using the conda version. To do so, please follow the instructions: [Install without conda](#install-without-conda)

* ganon runs on macOS but it is not available on conda due to compiler limitations. To install ganon os osx please follow the instructions: [Install without conda](#install-without-conda)

## Running ganon with sample data

### build

Create index/database from reference genomic sequences

```shh
ganon build --db-prefix sample_bacteria \
            --input-files tests/ganon/data/build/bacteria_*.fasta.gz \
            --taxdump-file tests/ganon/data/mini_nodes.dmp tests/ganon/data/mini_names.dmp \
            --seq-info-file tests/ganon/data/build/bacteria_seqinfo.txt
```

- `ganon build` (with a space in between) is different from the `ganon-build` command.
- `--taxdump-file` and `--seq-info-file` files are optional and will be automatically downloaded/generated if not provided

### classify

Classify reads against built index/database

```shh
ganon classify --db-prefix sample_bacteria --single-reads tests/ganon-classify/data/reads/bacteria.simulated.1.fq -o sample_results
```

### report

Generate readable reports, with filter options

```
ganon report --db-prefix sample_bacteria --rep-file sample_results.rep \
             --min-percentage 0.5 --ranks all \
             --output-prefix new_report -f tsv
```

### table

Generate a rank table (e.g.: samples X species) from multiple reports, with filter options

```
ganon table -i sample_results.tre new_report.tre -o sample_table.tsv
```

### update

Update index/database adding and/or removing reference genomic sequences

```
ganon update --db-prefix sample_bacteria --output-db-prefix sample_bateria_virus --input-files tests/ganon-build/data/sequences/virus*.fasta.gz
```

## Building custom indices

To build custom indices, ganon requires one (or multiple) fasta file(s) with the standard NCBI accession.version header (e.g. `>NC_009515.1`). For every sequence, taxonomic information will be automatically retrieved. Check the parameter `--seq-info` and `--taxdump-file` for more information. 

If you want to recreate the indices used in the manuscript above, please follow the instructions https://github.com/pirovc/ganon_benchmark

### Downloading sequences and building a new index

We suggest using [genome_updater](https://github.com/pirovc/genome_updater) (`conda install -c bioconda genome_updater`) to download sequences from RefSeq/Genbank. genome_updater can download and keep a subset (organism or taxonomic groups) of those repositories updated. For example:

Downloading archaeal and bacterial complete genomes from RefSeq: 

	genome_updater.sh -g "archaea,bacteria" \
	                  -d "refseq" \
	                  -l "Complete Genome" \
	                  -f "genomic.fna.gz,assembly_report.txt" \
                      -o "RefSeqCG_arc_bac" -b "v1" \
	                  -a -m -u -r -p -t 24

Where `-a` will additionally download the taxdump.tar.gz, `-m` will force the MD5 check for every file downloaded, `-u -r -p` will generate reports which can be used later, `-t` set the number of parallel downloads, `-b` will name the version and `-o` set the output working directory. If you want to download a set of defined species instead of the whole organism group, use `-g "species:562,623"` or any taxonomic group(s) `-g "taxids:620,1643685"`.

Building the index based on the files of the example above:

	ganon build --db-prefix ganon_db \
	            --input-files RefSeqCG_arc_bac/v1/files/*genomic.fna.gz

If you are getting the bash error `Argument list too long` use `--input-directory "RefSeqCG_arc_bac/v1/files/" --input-extension "genomic.fna.gz"` instead of `--input-files`

### Updating the index with new sequences

To update the folder with the most recent releases (after some days), just re-run the same download command:

	genome_updater.sh -g "archaea,bacteria" \
	                  -d "refseq" \
	                  -l "Complete Genome" \
	                  -f "genomic.fna.gz,assembly_report.txt" \
                      -o "RefSeqCG_arc_bac" -b "v2" \
	                  -a -m -u -r -p -t 24

This is going to download the new files (and remove the outdated ones) into the new label `v2`. New log and report files with the current version will be generated. It also checks if the last downloaded version is complete and downloads any missing file.

Updating the index based on the files of the example above:

	ganon update --db-prefix ganon_db --output-db-prefix ganon_db_updated \
	             --input-files $(find RefSeqCG_arc_bac/v2/files/ -name *genomic.fna.gz -type f)

If `--output-db-prefix` is not set, the database files `ganon_db.*` will be overwritten with the updated version. The `find` command will only look for files `-type f` inside the `v2/files/`, ignoring symbolic links from sequences which were not changing in this version.

By default, `ganon update` will only add new sequences provided to the index. To perform a full update, also removing sequences for the index, use the option `--update-complete` and provide the full set of updated references instead of only the new sequences, as below:

	ganon update --db-prefix ganon_db --output-db-prefix ganon_db_updated \
				 --update-complete \
	             --input-files RefSeqCG_arc_bac/v2/files/*genomic.fna.gz

### Using extra files to build and update

<details>
  <summary>Example</summary>

Optionally, some extra files generated by genome_updater can be further used to speed-up the building process:

	# Extract taxonomic information from genome_updater reports
	awk 'BEGIN {FS="\t";OFS="\t"}{if($4!="na"){ print $4,$5,$6,$2 }}' RefSeqCG_arc_bac/v1/updated_sequence_accession.txt > RefSeqCG_arc_bac/v1/seqinfo.txt
 
	# Use generated files from genome_updater on ganon build
	ganon build --db-prefix ganon_db \
	            --input-files RefSeqCG_arc_bac/v1/files/*genomic.fna.gz \
	            --seq-info-file RefSeqCG_arc_bac/v1/seqinfo.txt \
	            --taxdump-file RefSeqCG_arc_bac/v1/{TIMESTAMP}_taxdump.tar.gz

The same goes for the update process:

	# Extract taxonomic information from genome_updater reports
	awk 'BEGIN {FS="\t";OFS="\t"}{if($1=="A" && $4!="na"){ print $4,$5,$6,$2 }}' RefSeqCG_arc_bac/v2/updated_sequence_accession.txt > RefSeqCG_arc_bac/v2/seqinfo.txt

	# Use generated files on ganon update
	ganon update --db-prefix ganon_db \
	             --output-db-prefix ganon_db_updated \
	             --input-files $(find RefSeqCG_arc_bac/v2/files/ -name *genomic.fna.gz -type f) \
	             --seq-info-file RefSeqCG_arc_bac/v2/seqinfo.txt \
	             --taxdump-file RefSeqCG_arc_bac/v2/{TIMESTAMP}_taxdump.tar.gz

Obs:
-  If `-d genbank` was used in genome_updater, change the occurrences of `$4` to `$3` in the `awk` commands above.
- `{TIMESTAMP}` should be the timestamp (YYYY-MM-DD_HH-MM-SS) automatically generated when the files were downloaded.

</details>

## Output files

### build/update

Every run on `ganon build` or `ganon update` will generate the following database files:

 - {prefix}**.ibf**: main interleaved bloom filter file
 - {prefix}**.map**: tab-separated mapping between targets and bin identifiers. Targets should be present in the .tax file as a node *(fields: target, bin id)*
 - {prefix}**.tax**: taxonomic tree *(fields: node, parent, rank, name)*
 - {prefix}**.gnn**: gzipped pickled file (python) with information about clustering and parameters used

Obs:
-  Database files from version `0.1.X`, `0.2.X` and `0.3.X` are **NOT** compatible. If you want to convert a database to a newer version, please use the script `ganon-convert-db-0.1-0.2.py` or `ganon-convert-db-0.2-0.3.py`

### classify

 - {prefix}**.lca**: output with one match for each classified read after LCA. If multiple hierarchy levels are set, one file for each level will be created: {prefix}.{hierachy}.lca *(fields: read identifier, target, (max) k-mer count)*
 - {prefix}**.all**: output with all matches for each read. Only generated with --output-all/-a active. If multiple hierarchy levels are set, one file for each level will be created: {prefix}.{hierachy}.all. **Warning: file can be very large** *(fields: 1) read identifier, 2) target, 3) k-mer count)*
  - {prefix}**.rep**: plain report of the run with only targets that received a match *(fields: 1) hierarchy_label, 2) target, 3) total matches, 4) unique reads, 5) lca reads, 6) rank, 7) name)*. At the end prints 2 extra lines with `#total_classified` and `#total_unclassified`
  - {prefix}**.tre**: report file (see below)

### report

 - {prefix}**.tre**: tab-separated tree-like report with cumulative counts and taxonomic lineage. By default, this is a read-based report where each read classified is counted once. It is possible to generate this for all read matches (`ganon report --report-type matches`). In this case, single and shared matches are reported to their target. Each line in this report is a taxonomic entry, with the following fields: 

	1) taxonomic rank *(e.g. phylum, species, ...)*
	2) target *(e.g. taxid/assemblyid)*
	3) target lineage *(e.g 1|2|1224|...)*
	4) target name *(e.g. Paenibacillus polymyxa)*
	5) \# unique assignments *(number of reads that matched exclusively to this target)*
	6) \# reads/matches assignments *(number of reads/matches directly assigned to this target. Besides unique assignments, this number also includes lca assignments (in case of `--report-type reads`) or shared assignments (in case of `--report-type matches`))*
	7) \# cumulative assignments *(cumulative number of reads/matches assigned up-to this taxa)*
	8) \% cumulative assignments

- Using `--report-type reads` the first line of the file will show the number of unclassified reads

- When `--report-type reads` only taxa that received direct read matches, either unique or through lca, are considered. In cases where the classification is very ambiguous, some reads may have only shared matches and will not be reported. To look at those matches you can create a report with `--report-type matches` or look at the file {prefix}**.rep**.

### table

 - {output_file}: a tab-separated file with counts/percentages of taxa for multiple samples (rows) for a specific rank (cols)
 
### Examples of output files

The main output file is the `{prefix}.tre` which will summarize the results:

```
unclassified	      	                                     	unclassified           	0 	0 	1 	1.00000 
root        	1     	1                                    	root                   	0 	0 	99	99.00000
superkingdom	2     	1|2                                  	Bacteria               	0 	0 	99	99.00000
phylum      	1239  	1|2|1239                             	Firmicutes             	0 	0 	57	57.00000
phylum      	1224  	1|2|1224                             	Proteobacteria         	0 	0 	42	42.00000
class       	91061 	1|2|1239|91061                       	Bacilli                	0 	0 	57	57.00000
class       	28211 	1|2|1224|28211                       	Alphaproteobacteria    	0 	0 	29	29.00000
class       	1236  	1|2|1224|1236                        	Gammaproteobacteria    	0 	0 	13	13.00000
order       	1385  	1|2|1239|91061|1385                  	Bacillales             	0 	0 	57	57.00000
order       	204458	1|2|1224|28211|204458                	Caulobacterales        	0 	0 	29	29.00000
order       	72274 	1|2|1224|1236|72274                  	Pseudomonadales        	0 	0 	13	13.00000
family      	186822	1|2|1239|91061|1385|186822           	Paenibacillaceae       	0 	0 	57	57.00000
family      	76892 	1|2|1224|28211|204458|76892          	Caulobacteraceae       	0 	0 	29	29.00000
family      	468   	1|2|1224|1236|72274|468              	Moraxellaceae          	0 	0 	13	13.00000
genus       	44249 	1|2|1239|91061|1385|186822|44249     	Paenibacillus          	0 	0 	57	57.00000
genus       	75    	1|2|1224|28211|204458|76892|75       	Caulobacter            	0 	0 	29	29.00000
genus       	469   	1|2|1224|1236|72274|468|469          	Acinetobacter          	0 	0 	13	13.00000
species     	1406  	1|2|1239|91061|1385|186822|44249|1406	Paenibacillus polymyxa 	57	57	57	57.00000
species     	366602	1|2|1224|28211|204458|76892|75|366602	Caulobacter sp. K31    	29	29	29	29.00000
species     	470   	1|2|1224|1236|72274|468|469|470      	Acinetobacter baumannii	13	13	13	13.00000
```

running `ganon classify` or `ganon report` with `--ranks all`, the output will show all ranks used for classification and presented sorted by lineage:

```
unclassified 	       	                                                    	unclassified                                 	0 	0 	1 	1.00000 
root         	1      	1                                                   	root                                         	0 	0 	99	99.00000
no rank      	131567 	1|131567                                            	cellular organisms                           	0 	0 	99	99.00000
superkingdom 	2      	1|131567|2                                          	Bacteria                                     	0 	0 	99	99.00000
phylum       	1224   	1|131567|2|1224                                     	Proteobacteria                               	0 	0 	42	42.00000
class        	1236   	1|131567|2|1224|1236                                	Gammaproteobacteria                          	0 	0 	13	13.00000
order        	72274  	1|131567|2|1224|1236|72274                          	Pseudomonadales                              	0 	0 	13	13.00000
family       	468    	1|131567|2|1224|1236|72274|468                      	Moraxellaceae                                	0 	0 	13	13.00000
genus        	469    	1|131567|2|1224|1236|72274|468|469                  	Acinetobacter                                	0 	0 	13	13.00000
species group	909768 	1|131567|2|1224|1236|72274|468|469|909768           	Acinetobacter calcoaceticus/baumannii complex	0 	0 	13	13.00000
species      	470    	1|131567|2|1224|1236|72274|468|469|909768|470       	Acinetobacter baumannii                      	13	13	13	13.00000
class        	28211  	1|131567|2|1224|28211                               	Alphaproteobacteria                          	0 	0 	29	29.00000
order        	204458 	1|131567|2|1224|28211|204458                        	Caulobacterales                              	0 	0 	29	29.00000
family       	76892  	1|131567|2|1224|28211|204458|76892                  	Caulobacteraceae                             	0 	0 	29	29.00000
genus        	75     	1|131567|2|1224|28211|204458|76892|75               	Caulobacter                                  	0 	0 	29	29.00000
species      	366602 	1|131567|2|1224|28211|204458|76892|75|366602        	Caulobacter sp. K31                          	29	29	29	29.00000
no rank      	1783272	1|131567|2|1783272                                  	Terrabacteria group                          	0 	0 	57	57.00000
phylum       	1239   	1|131567|2|1783272|1239                             	Firmicutes                                   	0 	0 	57	57.00000
class        	91061  	1|131567|2|1783272|1239|91061                       	Bacilli                                      	0 	0 	57	57.00000
order        	1385   	1|131567|2|1783272|1239|91061|1385                  	Bacillales                                   	0 	0 	57	57.00000
family       	186822 	1|131567|2|1783272|1239|91061|1385|186822           	Paenibacillaceae                             	0 	0 	57	57.00000
genus        	44249  	1|131567|2|1783272|1239|91061|1385|186822|44249     	Paenibacillus                                	0 	0 	57	57.00000
species      	1406   	1|131567|2|1783272|1239|91061|1385|186822|44249|1406	Paenibacillus polymyxa                       	57	57	57	57.00000

```

## Hierarchical classification

Ganon classification can be performed in one or more databases at the same time. The databases can be provided in a hierarchical order. Multiple database classification can be performed providing several inputs for `--db-prefix`. They are required to be built with the same `k` size. To classify reads in a hierarchical order, `--hierarchy-labels` should be provided. `--max-error/--min-kmers` can be provided for each database. `--max-error-unique` can be provided for each hierarchy level. When using multiple hierarchical levels, output files will be generated for each level (use `--output-single` to generate a single output from multiple hierarchical levels).

For example: 

### Classifying reads against multiple databases:

	ganon classify --db-prefix db1 db2 db3 \
	               --min-kmers 0.75 \
	               -r reads.fq.gz

Classification against 3 database (as if they were one) using the same error rate.

### Classifying reads against multiple databases with different error rates:

	ganon classify --db-prefix db1 db2 db3 \
	               --max-error  0    1   4 \
	               -r reads.fq.gz

Classification against 3 database (as if they were one) using different error rates for each.

### Classifying reads against multiple databases hierarchically:

	ganon classify --db-prefix            db1     db2      db3 \ 
	               --hierarchy-labels 1_first 1_first 2_second \
	               -r reads.fq.gz

In this example, reads are going to be classified first against db1 and db2. Reads without a valid match will be further classified against db3. `--hierarchy-labels` are strings and are going to be sorted to define the hierarchy order, disregarding input order.

### Classifying reads against multiple databases hierarchically with different error rates:

	ganon classify --db-prefix            db1     db2      db3 \
	               --hierarchy-labels 1_first 1_first 2_second \
	               --min-kmers              1     0.5     0.25 \
	               --max-error-unique               0        1 \
	               -r reads.fq.gz

In this example, classification will be performed with different error rates for each database. For each hierarchy (`1_first` and `2_second`) a different `--max-error-unique` will be used.

## Choosing parameters

### ganon classify

#### --single-reads and --paired-reads

ganon accepts single-end or paired-end reads. In paired-end mode, reads are always reported with the header of the first pair. Paired-end reads are classified in a forward-reverse orientation.

#### --min-kmers and --max-error

Both parameters are used to define the similarity threshold between reads and references. `--max-error` will calculate the minimum amount of k-mers matches based on the q-gram lemma. `--min-kmers` will directly tell how many k-mers (in %) are necessary to consider a match. Note that the strata filter will always select reads with the best error range first (no error, 1 error, 2 errors, ...) and those parameters are controlling the lower bound of the threshold.

`--max-error` is recomended when working with precise classification. For example, most of your reads are represented in the index. `--max-error 0` means that all k-mers of a read should match a reference to be classified (very strict). Values used here are usually low (1, 3, 5), but not necessarily. Note that if you have reads of different lenghts you may want to use `--min-kmers` (with high values, e.g. 0.75) to apply roughly the same criteria for your data.

`--min-kmers` is recommended in more exploratory cases where exact matches are not possible. For example, analysing a sample with very few known species. `--min-kmers 0`  means that any read with one or more k-mers will be considered. This does not mean that any read will be classified, but that the threshold is very low for reads with few matches. Using low `--min-kmers` will mostly introduce false positives with a chance of increasing your sensitivity. Values here are usually low (0.25, 0.1, 0.05).

#### --max-error-unique

Exclusive error rate to define the similarity threshold of reads with unique matches. This is applied after filtering and only if a read is exclusively assigned to one target. If the classified read has less than `--max-error-unique`, the match is not excluded but assigned to its parent node. This is useful in a scenario when a read is poorly matched against a specific target (species, assembly) due to lack of representativity (just one references for a species, for example). Usually set to lower values (0, 1, 2).

#### --strata-filter

The strata filter is active by default. For every read, the best match - meaning most k-mers against the same target - is selected. An error value is calculated based on the q-gram lemma. All matches below this error rate are discarded (`--strata-filter 0` by default). For example: the best match of a read has 0 errors (all k-mers matched a target). With `--strata-filter 1` all matches with 0 + 1 errors will be reported. `--strata-filter -1` will disable the strata filtering and report everything up-to `--max-error`/`--min-kmers`. To use ganon as a not optimized k-mer counter use: `--min-kmers 0 --strata-filter -1`.

#### --offset

`--offset` can be used to speed-up analysis by skipping k-mers. `--offset 1` will check every k-mer of the sequences to be classified. `--offset n` will only evaluate every nth k-mer of the input sequences. For `--offset 1` there are possible performance improvements by disabling this function in compilation time with `-DGANON_OFFSET=OFF` (default is `ON`). Note that higher offset values will affect the sensitivity and precision of your classsification, specially when using 
`--min-kmers 0`.

### ganon build 

#### --max-bloom-size and --bin-length 

The most useful variable to define the IBF size (.ibf file) is the `--max-bloom-size`. It will set an approximate upper limit size for the file and estimate the `--bin-length` size based on it (ps: there is a minimum size necessary to generate the filter given a set of references and chosen parameters. Ganon will tell you if your value is too low.).

The IBF size is defined mainly on the amount of the input reference sequences (`-i`) but also can also be adjusted by a combination of parameters. Ganon will try to find the best `--bin-length` given `--max-fp`, `--kmer-size` and `--hash-functions`. Increasing `--max-fp` will generate smaller filters, but will generate more false positives in the classification step. If you know what you are doing, you can also directly set the size of the IBF with `--fixed-bloom-size` (ganon will tell you what's the resulting max. false positive).

`--bin-length` is the size in base pairs of each group for the taxonomic clustering (with TaxSBP). By default,  `--fragment-length` will be the size of `--bin-length` - `--overlap-length`, meaning that sequences will be split with overlap to fit into the bins. For example: species X has 2 sequences of 120bp each. Considering `--bin-length 50` and `--overlap-length 10` (`--fragment-length 40` consequently) each of the sequences will be split into 50bp and put into a bin with overlap of 10bp, resulting in 3 bins for each sequence (6 in total for species X).

Such adjustment is necessary to equalize the size of each bin, since the IBF requires the individual bloom filters to be of the same size by definition. Building the IBF based on the biggest sequence group in your references will generate the lowest number of bins but a very sparse and gigantic IBF. Building the IBF based on the smallest sequence group in your references will generate the smallest IBF but with too many bins. A balance between those two is necessary to achieve small and fast filters.

### ganon update

#### --update-complete

By default `ganon update` will only add sequences provided with `--input-files` to an previosly generated index. Using `--update-complete` it is possible to add and remove sequences from an index. When activating this option, ganon will consider that the files provided in `--input-files` are an actual representation of the index to build. It will automatically detect sequences that should be kept, inserted or removed given the input files and the information contained on the index to be updated.

### ganon report

#### --report-type

By default, `ganon classify` and  `ganon report` generate a read-based report (`ganon report --report-type reads`) where each read classified is counted once, either to its unique or lca assignment. It is possible to generate the same report for all read matches (`ganon report --report-type matches`). In this case, multiple matches for each reads are reported to their targets (single or shared matches). Using `--report-type matches` will not show the unclassified number of reads and it will always sum up to 100% in the root node, since this reports the overall distribution of matches and not the amount of reads classified.

#### --split-hierarchy, --keep-hierarchy or --skip-hierarchy

When using multiple databases in different hierarchical levels to classify reads, it is possible to report them separetly using `--split-hierarchy`. Once activated, one report will be generated for each hierarchical label. In each hierarchical specific report, the counts of other hierarchies are going to be reported at the root level to keep consistency in the reports.
It is also possible to select or ignore specific hierarchical labels (e.g. for a label use for host removal) using `--keep-hierarchy` or `--skip-hierarchy`.

## Install without conda

### build dependencies

System packages:
- gcc >=7 (check [gcc7 with conda](#installing-gcc7-in-a-separate-environment-with-conda)) or clang**
- cmake >=3.10
- zlib

** clang>=7 [linux] and AppleClang>=10.0.1 [osx: xcode-10.2/macOS 10.14]

### run dependencies

System packages:
- python >=3.4
- pandas >=0.22.0
- gawk
- grep
- tar
- curl
- wget
- coreutils (zcat)

** Please make sure that the system packages are supported/installed in your environment. All other packages are installed in the next steps.

### Installing taxsbp, binpacking and pylca

```shh
git clone https://github.com/pirovc/pylca.git
cd pylca
python3 setup.py install
```

```shh
pip3 install binpacking==1.4.3
binpacking -h
```

```shh
git clone https://github.com/pirovc/taxsbp.git
cd taxsbp
python3 setup.py install
taxsbp -h
```

### Downloading and building ganon

```shh
git clone --recurse-submodules https://github.com/pirovc/ganon.git # ganon, catch2, cxxopts, sdsl-lite, seqan
```
	
```shh
cd ganon
python3 setup.py install --record files.txt #optional
mkdir build_cpp
cd build_cpp
cmake -DCMAKE_BUILD_TYPE=Release -DVERBOSE_CONFIG=ON -DGANON_OFFSET=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCONDA=OFF ..
make
sudo make install #optional
```

- to change install location (e.g. `/myprefix/bin/`), set the installation prefix in the cmake command with `-DCMAKE_INSTALL_PREFIX=/myprefix/ `

- in the cmake command, set `-DGANON_OFFSET=ON` to be able to use the offset functionality

- use `-DINCLUDE_DIRS` to set alternative paths to cxxopts and Catch2 libs.

If everything was properly installed, the following commands should show the help pages without errors:

```shh
ganon -h
ganon-build -h
ganon-classify -h
```

### Run tests

```shh
python3 -m unittest discover -s tests/ganon/unit/
python3 -m unittest discover -s tests/ganon/integration/
python3 -m unittest discover -s tests/ganon/integration_online/
cd build_cpp/
ctest -VV .
```

## Installing GCC7 in a separate environment with conda

```shh
conda create -n gcc7 -c gouarin gcc-7 libgcc-7 "cmake>=3.10"
source activate gcc7
```

If you are getting the following error `ganon-classify: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version 'CXXABI_1.3.11' not found` you have to set the `LD_LIBRARY_PATH`

```shh
export LD_LIBRARY_PATH=/home/user/miniconda3/envs/gcc7/lib/
```

## Parameters

	usage: ganon [-h] [-v] {build,update,classify,report,table} ...

	ganon

	positional arguments:
	  {build,update,classify,report,table}
	    build               Build ganon database
	    update              Update ganon database
	    classify            Classify reads
	    report              Generate reports
	    table               Generate table from reports

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --version         Show program's version number and exit.

<details>
  <summary>ganon build</summary>

	usage: ganon build [-h] -d DB_PREFIX [-i [INPUT_FILES [INPUT_FILES ...]]] [-r] [-k] [-n] [-f] [-m] [-l] [-t]
	                   [--fixed-bloom-size] [--fragment-length] [--overlap-length] [--seq-info-mode [[...]]]
	                   [--seq-info-file] [--taxdump-file [[...]]] [--input-directory] [--input-extension]
	                   [--write-seq-info-file] [--verbose] [--quiet]

	optional arguments:
	  -h, --help            show this help message and exit
	  -r , --rank           Target taxonomic rank for classification [assembly,taxid,species,genus,...]. Default: species
	  -k , --kmer-size      The k-mer size for the interleaved bloom filter. Default: 19
	  -n , --hash-functions 
	                        The number of hash functions for the interleaved bloom filter. Default: 3
	  -f , --max-fp         Max. false positive rate for k-mer classification. Default: 0.05
	  -m , --max-bloom-size 
	                        Approx. maximum filter size in Megabytes (MB). Will estimate best --bin-length based on --kmer-
	                        size, --hash-functions and --max-fp [Mutually exclusive --fixed-bloom-size]
	  -l , --bin-length     Maximum length (in bp) for each bin. Default: auto
	  -t , --threads        Number of subprocesses/threads to use. Default: 2
	  --fixed-bloom-size    Fixed size for filter in Megabytes (MB), will ignore --max-fp [Mutually exclusive --max-bloom-
	                        size]
	  --fragment-length     Fragment length (in bp). Set to 0 to not fragment sequences. Default: --bin-length - --overlap-
	                        length
	  --overlap-length      Fragment overlap length (in bp). Should be bigger than the read length used for classification.
	                        Default: 300
	  --seq-info-mode [ [ ...]]
	                        Mode to obtain sequence information. For each sequence entry provided, ganon requires taxonomic
	                        and seq. length information. If a small number of sequences is provided (<50000) or when --rank
	                        assembly, ganon will automatically obtain data with NCBI E-utils websevices (eutils). Offline
	                        mode will download batch files from NCBI Taxonomy and look for taxonomic ids in the order
	                        provided. Options: [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot dead_nucl dead_wgs dead_prot],
	                        eutils (force webservices) or auto (uses eutils or [nucl_gb nucl_wgs]). Default: auto [Mutually
	                        exclusive --seq-info-file]
	  --seq-info-file       Pre-generated file with sequence information (seqid <tab> seq.len <tab> taxid [<tab> assembly
	                        id]) [Mutually exclusive --seq-info]
	  --taxdump-file [ [ ...]]
	                        Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp])
	                        file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)
	  --input-directory     Directory containing input files
	  --input-extension     Extension of files to use with --input-directory (provide it without * expansion, e.g. ".fna.gz")
	  --write-seq-info-file
	                        Write sequence information to DB_PREFIX.seqinfo.txt
	  --verbose             Verbose output mode
	  --quiet               Quiet output mode (only errors and warnings to the stderr)

	required arguments:
	  -d DB_PREFIX, --db-prefix DB_PREFIX
	                        Database output prefix (.ibf, .map, .tax, .gnn will be created)
	  -i [INPUT_FILES [INPUT_FILES ...]], --input-files [INPUT_FILES [INPUT_FILES ...]]
	                        Input reference sequence fasta files [.gz]

</details>

<details>
  <summary>ganon update</summary>

	usage: ganon update [-h] -d DB_PREFIX [-i [INPUT_FILES [INPUT_FILES ...]]] [-o] [-t] [--seq-info-mode [[...]]]
	                    [--seq-info-file] [--taxdump-file [[...]]] [--input-directory] [--input-extension]
	                    [--update-complete] [--write-seq-info-file] [--verbose] [--quiet]

	optional arguments:
	  -h, --help            show this help message and exit
	  -o , --output-db-prefix 
	                        Output database prefix (.ibf, .map, .tax, .gnn). Default: overwrite current --db-prefix
	  -t , --threads        Number of subprocesses/threads to use. Default: 2
	  --seq-info-mode [ [ ...]]
	                        Mode to obtain sequence information. For each sequence entry provided, ganon requires taxonomic
	                        and seq. length information. If a small number of sequences is provided (<50000) or when --rank
	                        assembly, ganon will automatically obtained data with NCBI E-utils websevices (eutils). Offline
	                        mode will download batch files from NCBI Taxonomy and look for taxonomic ids in the order
	                        provided. Options: [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot dead_nucl dead_wgs dead_prot],
	                        eutils (force webservices) or auto (uses eutils or [nucl_gb nucl_wgs]). Default: auto [Mutually
	                        exclusive --seq-info-file]
	  --seq-info-file       Pre-generated file with sequence information (seqid <tab> seq.len <tab> taxid [<tab> assembly
	                        id]) [Mutually exclusive --seq-info]
	  --taxdump-file [ [ ...]]
	                        Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp])
	                        file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)
	  --input-directory     Directory containing input files
	  --input-extension     Extension of files to use with --input-directory (provide it without * expansion, e.g. ".fna.gz")
	  --update-complete     Update adding and removing sequences. Input files should represent the complete updated set of
	                        references, not only new sequences.
	  --write-seq-info-file
	                        Write sequence information to DB_PREFIX.seqinfo.txt
	  --verbose             Verbose output mode
	  --quiet               Quiet output mode (only errors and warnings to the stderr)

	required arguments:
	  -d DB_PREFIX, --db-prefix DB_PREFIX
	                        Database input prefix (.ibf, .map, .tax, .gnn)
	  -i [INPUT_FILES [INPUT_FILES ...]], --input-files [INPUT_FILES [INPUT_FILES ...]]
	                        Input reference sequence fasta files [.gz] to be included to the database. Complete set of
	                        updated sequences should be provided when using --update-complete

</details>

<details>
  <summary>ganon classify</summary>

	usage: ganon classify [-h] -d [DB_PREFIX [DB_PREFIX ...]] [-s [reads.fq[.gz] [reads.fq[.gz] ...]]] [-p [reads.1.fq[.gz]
	                      reads.2.fq[.gz] [reads.1.fq[.gz] reads.2.fq[.gz] ...]]] [-o] [-c [[...]]] [-k [[...]]] [-e [[...]]]
	                      [-u [[...]]] [-l [[...]]] [-f] [-t] [-r [[...]]] [--output-all] [--output-unclassified]
	                      [--output-single] [--verbose] [--quiet]

	optional arguments:
	  -h, --help            show this help message and exit
	  -o , --output-prefix 
	                        Output prefix for .lca and .rep. Empty to output to STDOUT (only .lca will be printed)
	  -c [ [ ...]], --hierarchy-labels [ [ ...]]
	                        Hierarchy definition, one for each database input. Can also be a string, but input will be sorted
	                        to define order (e.g. 1 1 2 3). Default: 1
	  -k [ [ ...]], --min-kmers [ [ ...]]
	                        Min. percentage of k-mers matching to consider a read assigned. Single value or one per database
	                        (e.g. 0.5 0.7 1 0.25). Default: 0.25 [Mutually exclusive --max-error]
	  -e [ [ ...]], --max-error [ [ ...]]
	                        Max. number of errors allowed. Single value or one per database (e.g. 3 3 4 0) [Mutually
	                        exclusive --min-kmers]
	  -u [ [ ...]], --max-error-unique [ [ ...]]
	                        Max. number of errors allowed for unique assignments after filtering. Matches below this error
	                        rate will not be discarded, but assigned to a parent taxonomic level. Single value or one per
	                        hierarchy (e.g. 0 1 2). -1 to disable. Default: -1
	  -l [ [ ...]], --strata-filter [ [ ...]]
	                        Additional errors allowed (relative to the best match) to filter and select matches. Single value
	                        or one per hierarchy (e.g. 0 1 2). -1 to disable filtering. Default: 0
	  -f , --offset         Number of k-mers to skip during classification. Can speed up analysis but may reduce recall.
	                        (e.g. 1 = all k-mers, 3 = every 3rd k-mer). Default: 2
	  -t , --threads        Number of subprocesses/threads to use. Default: 3
	  -r [ [ ...]], --ranks [ [ ...]]
	                        Ranks to show in the report (.tre). "all" for all identified ranks. empty for default ranks:
	                        superkingdom phylum class order family genus species assembly. This file can be re-generated with
	                        the ganon report command.
	  --output-all          Output an additional file with all matches (.all). File can be very large.
	  --output-unclassified
	                        Output an additional file with unclassified read headers (.unc)
	  --output-single       When using multiple hierarchical levels, output everything in one file instead of one per
	                        hierarchy
	  --verbose             Verbose output mode
	  --quiet               Quiet output mode (only errors and warnings to the stderr)

	required arguments:
	  -d [DB_PREFIX [DB_PREFIX ...]], --db-prefix [DB_PREFIX [DB_PREFIX ...]]
	                        Database input prefix[es]
	  -s [reads.fq[.gz] [reads.fq[.gz] ...]], --single-reads [reads.fq[.gz] [reads.fq[.gz] ...]]
	                        Multi-fastq[.gz] file[s] to classify
	  -p [reads.1.fq[.gz] reads.2.fq[.gz] [reads.1.fq[.gz] reads.2.fq[.gz] ...]], --paired-reads [reads.1.fq[.gz] reads.2.fq[.gz] [reads.1.fq[.gz] reads.2.fq[.gz] ...]]
	                        Multi-fastq[.gz] pairs of file[s] to classify

</details>

<details>
  <summary>ganon report</summary>

	usage: ganon report [-h] [--min-count] [--min-percentage] [--names [[...]]] [--names-with [[...]]] [--taxids [[...]]] -i
	                    [REP_FILES [REP_FILES ...]] -o OUTPUT_PREFIX [-d [[...]]] [-e] [-r [[...]]] [-s] [-k [[...]]] [-f]
	                    [--taxdump-file [[...]]] [--verbose] [--quiet]

	optional arguments:
	  -h, --help            show this help message and exit
	  -d [ [ ...]], --db-prefix [ [ ...]]
	                        Database prefix[es] used for classification (in any order). If not provided, new taxonomy will be
	                        downloaded
	  -e , --report-type    Type of report to generate [reads, matches]. Default: reads
	  -r [ [ ...]], --ranks [ [ ...]]
	                        Fixer and ordered ranks for the report ["", "all", custom list] "all" for all possible ranks.
	                        empty for default ranks (superkingdom phylum class order family genus species assembly). Default:
	                        ""
	  -s , --sort           Sort report by [rank, lineage, count, unique]. Default: rank (with custom --ranks) or lineage
	                        (with --ranks all)
	  -k [ [ ...]], --skip-hierarchy [ [ ...]]
	                        One or more hierarchies to skip in the report (from ganon classify --hierarchy-labels)
	  -f , --output-format 
	                        Output format [text, tsv, csv]. Default: text
	  --taxdump-file [ [ ...]]
	                        Force use of a specific version of the (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp])
	                        file(s) from NCBI Taxonomy (otherwise it will be automatically downloaded)
	  --verbose             Verbose output mode
	  --quiet               Quiet output mode (only errors and warnings to the stderr)

	filter arguments:
	  --min-count           Mininum number of counts to keep the taxa
	  --min-percentage      Mininum percentage of counts to keep the taxa [0-1]
	  --names [ [ ...]]     Show only entries matching exact names of the provided list
	  --names-with [ [ ...]]
	                        Show entries containing full or partial names of the provided list
	  --taxids [ [ ...]]    One or more taxids to report (including children taxa)

	required arguments:
	  -i [REP_FILES [REP_FILES ...]], --rep-files [REP_FILES [REP_FILES ...]]
	                        One or more *.rep files from ganon classify
	  -o OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
	                        Output prefix for report (output_prefix.tre). In case of multiple files output is in the format
	                        "output_prefix + filename + .tre"

</details>

<details>
  <summary>ganon table</summary>

	usage: ganon table [-h] [--min-count] [--min-percentage] [--names [[...]]] [--names-with [[...]]] [--taxids [[...]]] -i
	                   [TRE_FILES [TRE_FILES ...]] -o OUTPUT_FILE [-l] [-f] [-t] [-a] [-r] [-m] [-p] [--add-unclassified]
	                   [--add-unclassified-rank] [--add-filtered] [--skip-zeros] [--verbose] [--quiet]

	optional arguments:
	  -h, --help            show this help message and exit
	  -l , --output-value   Output value on the table [percentage, counts]. percentage values are reported between [0-1].
	                        Default: percentage
	  -f , --output-format 
	                        Output format [tsv, csv]. Default: tsv
	  -t , --top-sample     Top hits of each sample individually
	  -a , --top-all        Top hits of all samples (ranked by percentage)
	  -r , --rank           Rank to report. Default: species
	  -m , --min-occurence 
	                        Number of occurences of a taxa among reports to be kept
	  -p , --min-occurence-percentage 
	                        Percentage of occurences of a taxa among reports to be kept [0-1]
	  --add-unclassified    Add column with unclassified count/percentage
	  --add-unclassified-rank
	                        Add column with unclassified count/percentage at the chosen rank but classified at a less
	                        specific rank
	  --add-filtered        Add column with filtered count/percentage
	  --skip-zeros          Do not print lines with only zero count/percentage
	  --verbose             Verbose output mode
	  --quiet               Quiet output mode (only errors and warnings to the stderr)

	filter arguments:
	  --min-count           Mininum number of counts to keep the taxa
	  --min-percentage      Mininum percentage of counts to keep the taxa [0-1]
	  --names [ [ ...]]     Show only entries matching exact names of the provided list
	  --names-with [ [ ...]]
	                        Show entries containing full or partial names of the provided list
	  --taxids [ [ ...]]    One or more taxids to report (including children taxa)

	required arguments:
	  -i [TRE_FILES [TRE_FILES ...]], --tre-files [TRE_FILES [TRE_FILES ...]]
	                        Report files (.tre) from ganon classify/report to make the table
	  -o OUTPUT_FILE, --output-file OUTPUT_FILE
	                        Output filename for the table

</details>