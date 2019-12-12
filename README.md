# ganon [![Build Status](https://travis-ci.org/pirovc/ganon.svg?branch=master)](https://travis-ci.org/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon)

ganon is a k-mer based read classification tool which uses Interleaved Bloom Filters in conjunction with a taxonomic clustering and a k-mer counting-filtering scheme. 

> **ganon: precise metagenomics classification against large and up-to-date sets of reference sequences**
> Vitor C. Piro, Temesgen H. Dadi, Enrico Seiler, Knut Reinert, Bernhard Y. Renard
> bioRxiv 406017; doi: [10.1101/406017](https://doi.org/10.1101/406017)

## Installation

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html)

```shh
conda install -c bioconda -c conda-forge ganon
```

* There are possible performance benefits compiling ganon from source rather than using the conda version. To do so, please follow the instructions at [manual installation](#manual-installation)

* Ganon runs on osx only with a [manual installation](#manual-installation). It was tested with gcc 7 and 8, but conda does not support those compilers for osx yet.

## Running ganon with sample data

### build

```shh
ganon build --db-prefix sample_bacteria --input-files tests/ganon-build/data/sequences/bacteria*.fasta.gz
```

Obs: `ganon build` (with a space in between) is different from the `ganon-build` command.

### classify

```shh
ganon classify --db-prefix sample_bacteria --single-reads tests/ganon-classify/data/reads/bacteria.simulated.1.fq -o sample_results
```

### report

```
ganon report --db-prefix sample_bacteria --rep-file sample_bateria.rep --min-matches-perc 0.5 --ranks all --output-report new_report.tre
```

### update

```
ganon update --db-prefix sample_bacteria --output-db-prefix sample_bateria_virus --input-files tests/ganon-build/data/sequences/virus*.fasta.gz
```

## Building custom indices

To build custom indices, ganon requires one (or multiple) fasta file(s) with the standard NCBI accession.version header (e.g. `>NC_009515.1`). For every sequence, taxonomic information will be automatically retrieved. Check the parameter `--seq-info` and `--taxdump-file` for more information. 

If you want to re-create the indices used in the manuscript above, please follow the instructions https://github.com/pirovc/ganon_benchmark

### Downloading sequences

We suggest using [genome_updater](https://github.com/pirovc/genome_updater) (`conda install -c bioconda genome_updater`) to download sequences from RefSeq/Genbank. genome_updater can download and keep a subset (organism or taxonomic groups) of those repositories updated. For example:

Downloading archaeal and bacterial complete genomes from RefSeq: 

	genome_updater.sh -g "archaea,bacteria" \
	                  -d "refseq" \
	                  -l "Complete Genome" \
	                  -f "genomic.fna.gz,assembly_report.txt" \
                      -o "RefSeqCG_arc_bac" -b "v1" \
	                  -a -m -u -r -p -t 24

Where `-a` will additionaly download the taxdump.tar.gz, `-m` will force the MD5 check for every file downloaded, `-u -r -p` will generate reports which can be used later, `-t` set the number of parallel downloads, `-b` will name the version and `-o` set the output working directory. If you want to download a set of defined species instead of the whole organism group, use `-g "species:562,623"` or any taxonomic group(s) `-g "taxids:620,1643685"`.

Building the index based on the files of the example above:

	ganon build --db-prefix ganon_db \
	            --input-files RefSeqCG_arc_bac/v1/files/*.genomic.fna.gz

If you are getting the bash error `Argument list too long` use `--input-directory "RefSeqCG_arc_bac/v1/files/" --input-extension ".genomic.fna.gz"` instead of `--input-files`

### Updating sequences

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
	             --input-files $(find RefSeqCG_arc_bac/v2/files/*.genomic.fna.gz -type f)

If `--output-db-prefix` is not set, the database files `ganon_db.*` will be overwritten with the updated version. The `find` command will only look for files `-type f` inside the `v2/files/`, ignoring symbolic links from sequences which were not changing in this version.

### Using extra files to build and update

Optionally, some extra files generated by genome_updater can be further used to speed-up the building process:

	# Extract taxonomic information from genome_updater reports
	awk 'BEGIN {FS="\t";OFS="\t"}{if($4!="na"){ print $4,$5,$6,$2 }}' RefSeqCG_arc_bac/v1/updated_sequence_accession.txt > RefSeqCG_arc_bac/v1/seqinfo.txt
 
	# Use generated files from genome_updater on ganon build
	ganon build --db-prefix ganon_db \
	            --input-files RefSeqCG_arc_bac/v1/files/*.genomic.fna.gz \
	            --seq-info-file RefSeqCG_arc_bac/v1/seqinfo.txt \
	            --taxdump-file RefSeqCG_arc_bac/v1/{TIMESTAMP}_taxdump.tar.gz

The same goes for the update process:

	# Extract taxonomic information from genome_updater reports
	awk 'BEGIN {FS="\t";OFS="\t"}{if($1=="A" && $4!="na"){ print $4,$5,$6,$2 }}' RefSeqCG_arc_bac/v2/updated_sequence_accession.txt > RefSeqCG_arc_bac/v2/seqinfo.txt

	# Use generated files on ganon update
	ganon update --db-prefix ganon_db \
	             --output-db-prefix ganon_db_updated \
	             --input-files $(find RefSeqCG_arc_bac/v2/files/*.genomic.fna.gz -type f) \
	             --seq-info-file RefSeqCG_arc_bac/v2/seqinfo.txt \
	             --taxdump-file RefSeqCG_arc_bac/v2/{TIMESTAMP}_taxdump.tar.gz

Obs:
-  If `-d genbank` was used in genome_updater, change the occurances of `$4` to `$3` in the `awk` commands above.
- `{TIMESTAMP}` should be the timestamp (YYYY-MM-DD_HH-MM-SS) automatically generated when the files were downloaded.

## Output files

### build/update

Every run on `ganon build` or `ganon update` will generate the following database files:

 - {prefix}**.ibf**: main interleaved bloom filter file
 - {prefix}**.map**: tab-separated mapping between targets and bin identifiers (fields: target, bin id). Targets should be present in the .tax file as a node.
 - {prefix}**.tax**: taxonomic tree (fields: node, parent, rank, name)
 - {prefix}**.gnn**: gziped pickled file (python) with information about clustering and parameters used

Obs:
-  Database files from version `0.1.X` are **NOT** compatible with `0.2.X`. If you want to convert a database, please use the script `scripts/convert-db-0.1-0.2.py`

### classify

 - {prefix}**.lca**: output with one match for each classified read after LCA. If multiple hiearchy levels are set, one file for each level will be created: {prefix}.{hierachy}.lca (fields: read identifier, target, (max) k-mer count)
 - {prefix}**.all**: output with all matches for each read. Only generated with --output-all/-a active. If multiple hiearchy levels are set, one file for each level will be created: {prefix}.{hierachy}.all. *Warning: file can be very large* (fields: read identifier, target, k-mer count)
  - {prefix}**.rep**: plain report of the run with only target that receive an match (fields: hierarchy_label, target, total matches, unique matches, lca matches, rank, name). At the end prints 2 extra lines with `#total_classified` and `#total_unclassified`
  - {prefix}**.tre**: report file (see below)

### report

 - {prefix}**.tre**: tab-separated tree-like report with cummulative counts and lineage with the following fields: 

	1) rank (e.g. phylum, species, ...)
	2) target (e.g. taxid, assembly id)
	3) taxid lineage (e.g 1|2|1224|...)
	4) name (e.g. Paenibacillus polymyxa)
	5) # unique assignments (number of reads that matched exclusively to this target)
	6) # reads assigned (number of reads directly assigned to this target - either unique or lca)
	7) # cummulative assignments (cummulative reads assigned up-to this taxa)
	8) % cummulative assignments

#### Example of classification output files

The main output file is the `{prefix}.tre` which will sumarize the results (too see the file in a nice format, use `column -s$'\t' -t -n {prefix}.tre`):

```
unclassified  -       -                                      -                        -   -   0    0.00000
root          1       1                                      root                     0   0   100  100.00000
superkingdom  2       1|2                                    Bacteria                 0   0   100  100.00000
phylum        1239    1|2|1239                               Firmicutes               0   0   58   58.00000
phylum        1224    1|2|1224                               Proteobacteria           0   0   42   42.00000
class         91061   1|2|1239|91061                         Bacilli                  0   0   58   58.00000
class         28211   1|2|1224|28211                         Alphaproteobacteria      0   0   29   29.00000
class         1236    1|2|1224|1236                          Gammaproteobacteria      0   0   13   13.00000
order         1385    1|2|1239|91061|1385                    Bacillales               0   0   58   58.00000
order         204458  1|2|1224|28211|204458                  Caulobacterales          0   0   29   29.00000
order         72274   1|2|1224|1236|72274                    Pseudomonadales          0   0   13   13.00000
family        186822  1|2|1239|91061|1385|186822             Paenibacillaceae         0   0   58   58.00000
family        76892   1|2|1224|28211|204458|76892            Caulobacteraceae         0   0   29   29.00000
family        468     1|2|1224|1236|72274|468                Moraxellaceae            0   0   13   13.00000
genus         44249   1|2|1239|91061|1385|186822|44249       Paenibacillus            0   0   58   58.00000
genus         75      1|2|1224|28211|204458|76892|75         Caulobacter              0   0   29   29.00000
genus         469     1|2|1224|1236|72274|468|469            Acinetobacter            0   0   13   13.00000
species       1406    1|2|1239|91061|1385|186822|44249|1406  Paenibacillus polymyxa   58  58  58   58.00000
species       366602  1|2|1224|28211|204458|76892|75|366602  Caulobacter sp. K31      29  29  29   29.00000
species       470     1|2|1224|1236|72274|468|469|470        Acinetobacter baumannii  13  13  13   13.00000
```

running `ganon classify` or `ganon report` with `--ranks all`, the output will show all ranks used for classification, sorted by lineage:

```
unclassified   -        -                                                     -                                              -   -   0    0.00000
no rank        1        1                                                     root                                           0   0   100  100.00000
no rank        131567   1|131567                                              cellular organisms                             0   0   100  100.00000
superkingdom   2        1|131567|2                                            Bacteria                                       0   0   100  100.00000
phylum         1224     1|131567|2|1224                                       Proteobacteria                                 0   0   42   42.00000
class          1236     1|131567|2|1224|1236                                  Gammaproteobacteria                            0   0   13   13.00000
order          72274    1|131567|2|1224|1236|72274                            Pseudomonadales                                0   0   13   13.00000
family         468      1|131567|2|1224|1236|72274|468                        Moraxellaceae                                  0   0   13   13.00000
genus          469      1|131567|2|1224|1236|72274|468|469                    Acinetobacter                                  0   0   13   13.00000
species group  909768   1|131567|2|1224|1236|72274|468|469|909768             Acinetobacter calcoaceticus/baumannii complex  0   0   13   13.00000
species        470      1|131567|2|1224|1236|72274|468|469|909768|470         Acinetobacter baumannii                        13  13  13   13.00000
class          28211    1|131567|2|1224|28211                                 Alphaproteobacteria                            0   0   29   29.00000
order          204458   1|131567|2|1224|28211|204458                          Caulobacterales                                0   0   29   29.00000
family         76892    1|131567|2|1224|28211|204458|76892                    Caulobacteraceae                               0   0   29   29.00000
genus          75       1|131567|2|1224|28211|204458|76892|75                 Caulobacter                                    0   0   29   29.00000
no rank        2648921  1|131567|2|1224|28211|204458|76892|75|2648921         unclassified Caulobacter                       0   0   29   29.00000
species        366602   1|131567|2|1224|28211|204458|76892|75|2648921|366602  Caulobacter sp. K31                            29  29  29   29.00000
no rank        1783272  1|131567|2|1783272                                    Terrabacteria group                            0   0   58   58.00000
phylum         1239     1|131567|2|1783272|1239                               Firmicutes                                     0   0   58   58.00000
class          91061    1|131567|2|1783272|1239|91061                         Bacilli                                        0   0   58   58.00000
order          1385     1|131567|2|1783272|1239|91061|1385                    Bacillales                                     0   0   58   58.00000
family         186822   1|131567|2|1783272|1239|91061|1385|186822             Paenibacillaceae                               0   0   58   58.00000
genus          44249    1|131567|2|1783272|1239|91061|1385|186822|44249       Paenibacillus                                  0   0   58   58.00000
species        1406     1|131567|2|1783272|1239|91061|1385|186822|44249|1406  Paenibacillus polymyxa                         58  58  58   58.00000

```

## Choosing parameters

### --single-reads and --paired-reads

ganon accepts single-end or paired-end reads. In the paired-end mode, reads are always reported with the header of the first pair. The maximum number of k-mers matches a pair can have is: `length(read1) + length(read2) + 1 - k`. Paired-end reads are classified in a forward-reverse orientation.

### --min-kmers and --max-error

Both parameters are used to define the similarity threshold between reads and references. `--max-error` will work with fixed number of errors to calculate the amount of k-mers necessary to match. `--min-kmers` will directly tell how many k-mers (in %) are necessary to consider a match.

### --max-error-unique

Exclusive error rate to define the similarity threshold of reads with unique matches. This is applied after filtering and only if a read is exclusively assigned to one target. If the threshold is not achieved, the match is not excluded but assigned to the parent node. This is useful in a scenario when a read is poorly matched against a specific target (species, assembly) due to lack of representativity (just one references for a species, for example). Usually set to 0 or 1.

### .ibf size (--max-bloom-size, --bin-length)

The most useful variable to define the IBF size (.ibf file) is the `--max-bloom-size`. It will set an approximate upper limit size for the file and estimate the `--bin-length` size based on it (ps: there is a minimum size necessary to generate the filter given a set of references and chosen parameters. Ganon will tell you if your value is too low.).

The IBF size is defined mainly on the amount of the input reference sequences (`-i`) but also can also be adjusted by a combination of parameters. Ganon will try to find the best `--bin-length` given `--max-fp`, `--kmer-size` and `--hash-functions`. Increasing `--max-fp` will generate smaller filters, but will generate more false positives in the classification step. If you know what you are doing, you can also directly set the size of the IBF with `--fixed-bloom-size` (ganon will tell you what's the resulting max. false positive).

`--bin-length` is the size in base pairs of each group for the taxonomic clustering (with TaxSBP). By default,  `--fragment-length` will be the size of `--bin-length` - `--overlap-length`, meaning that sequences will be split with overlap to fit into the bins. For example: species X has 2 sequences of 120bp each. Considering `--bin-length 50` and `--overlap-length 10` (`--fragment-length 40` consequently) each of the sequences will be split into 50bp and put into a bin with overlap of 10bp, resulting in 3 bins for each sequence (6 in total for species X).

Such adjustment is necessary to equalize the size of each bin, since the IBF requires the individual bloom filters to be of the same size by definition. Building the IBF based on the biggest sequence group in your references will generate the lowest number of bins but a very sparse and gigantic IBF. Building the IBF based on the smallest sequence group in your references will generate the smallest IBF but with too many bins. A balance between those two is necessary to achieve small and fast filters.

## Manual Installation

### Dependencies

#### build

System packages:
- gcc >=7 (check [gcc7 with conda](#installing-gcc7-in-a-separate-environment-with-conda)) or clang>=7
- cmake >=3.10

Specific packages:
- Catch2 >=2.7.0 ([d63307](https://github.com/catchorg/Catch2/commit/d63307279412de3870cf97cc6802bae8ab36089e))
- cxxopts >=2.2.0 ([a0de9f](https://github.com/jarro2783/cxxopts/commit/073dd3e645fa0c853c3836f3788ca21c39af319d))
- sdsl-lite 3.0 ([d6ed14](https://github.com/xxsds/sdsl-lite/commit/d6ed14d5d731ed4a4ec12627c1ed7154b396af48))
- seqan 2.4.0 ([c308e9](https://github.com/eseiler/seqan/commit/c308e99f10d942382d4c7ed6fc91be1a889e644c))

#### run

System packages:
- python >=3.4
- pandas
- wget
- curl
- tar
- GNU core utilities (gawk, zcat)

Specific packages:
- taxsbp >=0.1.2 ([6e1481](https://github.com/pirovc/taxsbp/commit/6e14819791a960273a191b3e1e028b084ed2d945))
- binpacking >=1.4.1 ([v1.4.1](https://pypi.org/project/binpacking/1.4.1/))

** Please make sure that the system packages are supported/installed in your environment. All other packages are installed in the next steps.

### Obtaining packages

```shh
git clone --recurse-submodules https://github.com/pirovc/ganon.git # ganon, catch2, cxxopts, sdsl-lite, seqan
git clone https://github.com/pirovc/taxsbp.git # taxsbp
```

### Installing 

#### taxsbp

```shh
cd taxsbp
python3 setup.py install
taxsbp -h
```

#### binpacking

```shh
pip3 install binpacking==1.4.1
binpacking -h
```

### Building
	
```shh
cd ganon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DVERBOSE_CONFIG=ON -DGANON_OFFSET=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCONDA=OFF ..
make
```

in the cmake command, set `-DGANON_OFFSET=ON` to be able to use the offset functionality. use `-DINCLUDE_DIRS` to set alternative paths to cxxopts and Catch2 libs.

### Testing

```shh
cd ganon
./ganon -h
python3 -m unittest discover -s tests/ganon/unit/
python3 -m unittest discover -s tests/ganon/integration/

cd build
./ganon-build -h
./ganon-classify -h
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

### build

### classify

### update

