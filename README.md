# ganon [![Build Status](https://travis-ci.org/pirovc/ganon.svg?branch=master)](https://travis-ci.org/pirovc/ganon) [![codecov](https://codecov.io/gh/pirovc/ganon/branch/master/graph/badge.svg)](https://codecov.io/gh/pirovc/ganon)

ganon is a k-mer based read classification tool which uses Interleaved Bloom Filters in conjunction with a taxonomic clustering and a k-mer counting-filtering scheme. 

> **ganon: continuously up-to-date with database growth for precise short read classification in metagenomics**
> Vitor C. Piro, Temesgen H. Dadi, Enrico Seiler, Knut Reinert, Bernhard Y. Renard
> bioRxiv 406017; doi: [10.1101/406017](https://doi.org/10.1101/406017)

## Installation

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/ganon/README.html)

```shh
conda install -c bioconda ganon
```

* There are possible performance benefits compiling ganon from source rathen than using the conda version. To do so, please follow the instructions at [manual installation](#manual-installation)

* Ganon run on MacOS only with a [manual installation](#manual-installation). It was tested with gcc/clang 7 and 8, but conda does not support those compilers for mac yet.

## Running ganon with sample data

### build

```shh
ganon build --db-prefix sample_bacteria --input-files tests/ganon-build/data/sequences/bacteria*.fasta.gz
```

### classify

```shh
ganon classify --db-prefix sample_bacteria --reads tests/ganon-classify/data/reads/bacteria.simulated.1.fq -o sample_results
```

### update

```
ganon update --db-prefix sample_bacteria --output-db-prefix sample_bateria_virus --input-files tests/ganon-build/data/sequences/virus*.fasta.gz
```

## Parameters

### build

	$ ganon build --help
	usage: ganon build [-h] -d db_prefix -i [refs.fasta[.gz] [refs.fasta[.gz]
	                   ...]] [-r] [-k] [-n] [-f] [-m] [-l] [-t]
	                   [--fixed-bloom-size] [--fragment-length] [--overlap-length]
	                   [--seq-info [[...]]] [--seq-info-file]
	                   [--taxdump-file [[...]]] [--verbose]

	optional arguments:
	  -h, --help            show this help message and exit
	  -r , --rank           Lowest taxonomic rank for classification
	                        [assembly,taxid,species,genus,...]. Default: species
	  -k , --kmer-size      The k-mer size for the bloom filter [14..32]. Default:
	                        19
	  -n , --hash-functions
	                        The number of hash functions to use for the bloom
	                        filter [2..5]. Default: 3
	  -f , --max-fp         Max. false positive rate for k-mer classification.
	                        Default: 0.05
	  -m , --max-bloom-size
	                        Approx. maximum filter size in Megabytes (MB). Will
	                        estimate best --bin-length based on --kmer-size,
	                        --hash-functions and --max-fp [Mutually exclusive
	                        --fixed-bloom-size]
	  -l , --bin-length     Maximum length (in bp) for each bin. Default: auto
	  -t , --threads        Number of subprocesses/threads to use for
	                        calculations. Default: 2
	  --fixed-bloom-size    Fixed size for filter in Megabytes (MB), will ignore
	                        --max-fp [Mutually exclusive --max-bloom-size]
	  --fragment-length     Fragment length (in bp). Set to 0 to not fragment
	                        sequences. Default: --bin-length - --overlap-length
	  --overlap-length      Fragment overlap length (in bp). Should be bigger than
	                        the read length used for classification. Default: 300
	  --seq-info [ [ ...]]  Mode to obtain sequence information. For each sequence
	                        entry provided with --input-files, ganon requires
	                        taxonomic and seq. length information. If a small
	                        number of sequences is provided (<50000) or when
	                        --rank assembly, ganon will automatically obtained
	                        data with NCBI E-utils websevices (eutils). Offline
	                        mode will download batch files from NCBI Taxonomy and
	                        look for taxonomic ids in the order provided. Options:
	                        [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot dead_nucl
	                        dead_wgs dead_prot], eutils (force webservices) or
	                        auto (uses eutils or [nucl_gb nucl_wgs]). Default:
	                        auto [Mutually exclusive --seq-info-file]
	  --seq-info-file       Pre-generated file with sequence information (seqid
	                        <tab> seq.len <tab> taxid [<tab> assembly id])
	                        [Mutually exclusive --seq-info]
	  --taxdump-file [ [ ...]]
	                        Force use of a specific version of the
	                        (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp])
	                        file(s) from NCBI Taxonomy (otherwise it will be
	                        automatically downloaded)
	  --verbose             Verbose mode for ganon

	required arguments:
	  -d db_prefix, --db-prefix db_prefix
	                        Database output prefix (.filter, .nodes, .bins. .map
	                        will be created)
	  -i [refs.fasta[.gz] [refs.fasta[.gz] ...]], --input-files [refs.fasta[.gz] [refs.fasta[.gz] ...]]
	                        Multi-fasta[.gz] file[s]

### classify

	$ ganon classify --help
	usage: ganon classify [-h] -d [db_prefix [db_prefix ...]] -r [reads.fq[.gz]
	                      [reads.fq[.gz] ...]] [-c [int [int ...]]]
	                      [-m [int [int ...]]] [-e [int [int ...]]]
	                      [-u [int [int ...]]] [-f] [-o] [-n] [-s] [-l] [-p]
	                      [-k [RANKS [RANKS ...]]] [-t] [--verbose]

	optional arguments:
	  -h, --help            show this help message and exit
	  -c [int [int ...]], --db-hierarchy [int [int ...]]
	                        Hierachy definition, one for each database input. Can
	                        also be string, but input will be always sorted (e.g.
	                        1 1 2 3). Default: 1
	  -m [int [int ...]], --min-kmers [int [int ...]]
	                        Min. percentage of k-mers matching to consider a read
	                        assigned. Can be used alternatively to --max-error for
	                        reads of variable size. Single value or one per
	                        database (e.g. 0.5 0.7 1 0.25). Default: 0.25
	  -e [int [int ...]], --max-error [int [int ...]]
	                        Max. number of errors allowed. Single value or one per
	                        database (e.g. 3 3 4 0) [Mutually exclusive --min-
	                        kmers]
	  -u [int [int ...]], --max-error-unique [int [int ...]]
	                        Max. number of errors allowed for unique assignments
	                        after filtering. Matches below this error rate will
	                        not be discarded, but assigned to parent taxonomic
	                        level. Single value or one per hierachy (e.g. 0 1 2).
	                        -1 to disable. Default: -1
	  -f , --offset         Number of k-mers to skip during clasification. Can
	                        speed up analysis but may reduce recall. (e.g. 1 = all
	                        k-mers, 3 = every 3rd k-mer). Function must be enabled
	                        on compilation time with -DGANON_OFFSET=ON. Default: 2
	  -o , --output-file-prefix
	                        Output file name prefix: .out for complete results /
	                        .lca for LCA results / .rep for report. Empty to print
	                        to STDOUT (only with lca). Default: ""
	  -n , --output-unclassified-file
	                        Output file for unclassified reads headers. Empty to
	                        not output. Default: ""
	  -s, --split-output-file-hierarchy
	                        Split output in multiple files by hierarchy. Appends
	                        "_hierachy" to the --output-file definiton.
	  -l, --skip-lca        Skip LCA step and output multiple matches. --max-
	                        error-unique will not be applied
	  -p, --skip-reports    Skip reports
	  -k [RANKS [RANKS ...]], --ranks [RANKS [RANKS ...]]
	                        Ranks for the final report. "all" for all indentified
	                        ranks. empty for default ranks: superkingdom phylum
	                        class order family genus species species+ assembly
	  -t , --threads        Number of subprocesses/threads. Default: 3)
	  --verbose             Output in verbose mode for ganon-classify

	required arguments:
	  -d [db_prefix [db_prefix ...]], --db-prefix [db_prefix [db_prefix ...]]
	                        Database prefix[es]
	  -r [reads.fq[.gz] [reads.fq[.gz] ...]], --reads [reads.fq[.gz] [reads.fq[.gz] ...]]
	                        Multi-fastq[.gz] file[s] to classify

### update

	$ ganon update --help
	usage: ganon update [-h] -d db_prefix -i [refs.fasta[.gz] [refs.fasta[.gz]
	                    ...]] [-o] [-t] [--seq-info [[...]]] [--seq-info-file]
	                    [--taxdump-file [[...]]] [--verbose]

	optional arguments:
	  -h, --help            show this help message and exit
	  -o , --output-db-prefix
	                        Alternative output database prefix. Default: overwrite
	                        current --db-prefix
	  -t , --threads        set the number of subprocesses/threads to use for
	                        calculations. Default: 2
	  --seq-info [ [ ...]]  Mode to obtain sequence information. For each sequence
	                        entry provided with --input-files, ganon requires
	                        taxonomic and seq. length information. If a small
	                        number of sequences is provided (<50000) or when
	                        --rank assembly, ganon will automatically obtained
	                        data with NCBI E-utils websevices (eutils). Offline
	                        mode will download batch files from NCBI Taxonomy and
	                        look for taxonomic ids in the order provided. Options:
	                        [nucl_gb nucl_wgs nucl_est nucl_gss pdb prot dead_nucl
	                        dead_wgs dead_prot], eutils (force webservices) or
	                        auto (uses eutils or [nucl_gb nucl_wgs]). Default:
	                        auto [Mutually exclusive --seq-info-file]
	  --seq-info-file       Pre-generated file with sequence information (seqid
	                        <tab> seq.len <tab> taxid [<tab> assembly id])
	                        [Mutually exclusive --seq-info]
	  --taxdump-file [ [ ...]]
	                        Force use of a specific version of the
	                        (taxdump.tar.gz) or (nodes.dmp names.dmp [merged.dmp])
	                        file(s) from NCBI Taxonomy (otherwise it will be
	                        automatically downloaded)
	  --verbose             Verbose mode for ganon

	required arguments:
	  -d db_prefix, --db-prefix db_prefix
	                        Database prefix
	  -i [refs.fasta[.gz] [refs.fasta[.gz] ...]], --input-files [refs.fasta[.gz] [refs.fasta[.gz] ...]]
	                        Multi-fasta[.gz] file[s]

## Output files

The main output file is the `*.tre` which will sumarize the results.

Example of the sample data classification `sample_results.tre` previously generated:

```
unclassified  -        -                                              -                               1   1.00000
root          1        1                                              root                            99  99.00000
superkingdom  2        1|2                                            Bacteria                        99  99.00000
phylum        1239     1|2|1239                                       Firmicutes                      57  57.00000
phylum        1224     1|2|1224                                       Proteobacteria                  42  42.00000
class         91061    1|2|1239|91061                                 Bacilli                         57  57.00000
class         28211    1|2|1224|28211                                 Alphaproteobacteria             29  29.00000
class         1236     1|2|1224|1236                                  Gammaproteobacteria             13  13.00000
order         1385     1|2|1239|91061|1385                            Bacillales                      57  57.00000
order         204458   1|2|1224|28211|204458                          Caulobacterales                 29  29.00000
order         72274    1|2|1224|1236|72274                            Pseudomonadales                 13  13.00000
family        186822   1|2|1239|91061|1385|186822                     Paenibacillaceae                57  57.00000
family        76892    1|2|1224|28211|204458|76892                    Caulobacteraceae                29  29.00000
family        468      1|2|1224|1236|72274|468                        Moraxellaceae                   13  13.00000
genus         44249    1|2|1239|91061|1385|186822|44249               Paenibacillus                   57  57.00000
genus         75       1|2|1224|28211|204458|76892|75                 Caulobacter                     29  29.00000
genus         469      1|2|1224|1236|72274|468|469                    Acinetobacter                   13  13.00000
species       1406     1|2|1239|91061|1385|186822|44249|1406          Paenibacillus polymyxa          57  57.00000
species       366602   1|2|1224|28211|204458|76892|75|366602          Caulobacter sp. K31             29  29.00000
species       470      1|2|1224|1236|72274|468|469|470                Acinetobacter baumannii         13  13.00000
species+      1052684  1|2|1239|91061|1385|186822|44249|1406|1052684  Paenibacillus polymyxa M1       57  57.00000
species+      696749   1|2|1224|1236|72274|468|469|470|696749         Acinetobacter baumannii 1656-2  13  13.00000
```

running with `--ranks all`, the output will show all ranks used for classification, sorted by lineage:

```
unclassified   -        -                                                             -                                              1   1.00000
no rank        1        1                                                             root                                           99  99.00000
no rank        131567   1|131567                                                      cellular organisms                             99  99.00000
superkingdom   2        1|131567|2                                                    Bacteria                                       99  99.00000
phylum         1224     1|131567|2|1224                                               Proteobacteria                                 42  42.00000
class          1236     1|131567|2|1224|1236                                          Gammaproteobacteria                            13  13.00000
order          72274    1|131567|2|1224|1236|72274                                    Pseudomonadales                                13  13.00000
family         468      1|131567|2|1224|1236|72274|468                                Moraxellaceae                                  13  13.00000
genus          469      1|131567|2|1224|1236|72274|468|469                            Acinetobacter                                  13  13.00000
species group  909768   1|131567|2|1224|1236|72274|468|469|909768                     Acinetobacter calcoaceticus/baumannii complex  13  13.00000
species        470      1|131567|2|1224|1236|72274|468|469|909768|470                 Acinetobacter baumannii                        13  13.00000
no rank        696749   1|131567|2|1224|1236|72274|468|469|909768|470|696749          Acinetobacter baumannii 1656-2                 13  13.00000
class          28211    1|131567|2|1224|28211                                         Alphaproteobacteria                            29  29.00000
order          204458   1|131567|2|1224|28211|204458                                  Caulobacterales                                29  29.00000
family         76892    1|131567|2|1224|28211|204458|76892                            Caulobacteraceae                               29  29.00000
genus          75       1|131567|2|1224|28211|204458|76892|75                         Caulobacter                                    29  29.00000
species        366602   1|131567|2|1224|28211|204458|76892|75|366602                  Caulobacter sp. K31                            29  29.00000
no rank        1783272  1|131567|2|1783272                                            Terrabacteria group                            57  57.00000
phylum         1239     1|131567|2|1783272|1239                                       Firmicutes                                     57  57.00000
class          91061    1|131567|2|1783272|1239|91061                                 Bacilli                                        57  57.00000
order          1385     1|131567|2|1783272|1239|91061|1385                            Bacillales                                     57  57.00000
family         186822   1|131567|2|1783272|1239|91061|1385|186822                     Paenibacillaceae                               57  57.00000
genus          44249    1|131567|2|1783272|1239|91061|1385|186822|44249               Paenibacillus                                  57  57.00000
species        1406     1|131567|2|1783272|1239|91061|1385|186822|44249|1406          Paenibacillus polymyxa                         57  57.00000
no rank        1052684  1|131567|2|1783272|1239|91061|1385|186822|44249|1406|1052684  Paenibacillus polymyxa M1                      57  57.00000
```

### classify

#### .out

all matches (one per hierarchy with `-s, --split-output-file-hierarchy`)

	readid <tab> assignment <tab> k-mer count

* a negative k-mer count work a flag to show that such read did not have the minimum amount of matches to pass the `--max-error-unique` threshold

#### .lca

only one match / read (one per hierarchy with `-s, --split-output-file-hierarchy`)
	
	readid <tab> lca assignment <tab> max k-mer count

#### .rep

detailed report for lca matches (one per hierarchy with `-s, --split-output-file-hierarchy`)
	
	1) lca assignment <tab>
	2) reads assigned (lca) <tab>
	3) % reads assigned <tab>
	4) # reads assigned (total) <tab>
	5) # reads uniquely assigned <tab>
	6) taxonomic rank <tab>
	7) name

#### .tre

tree-like output with cummulative counts and lineage (one per run)
	
	1) rank <tab>
	2) lca assignment <tab>
	3) lineage <tab>
	4) name <tab>
	5) cummulative # reads assigned <tab>
	6) cummulative % reads assigned

### build/update

#### .filter
main bloom filter
#### .map
maps assignment to binids
#### .bins
bins generated by TaxSBP
#### .nodes
taxonomic information required for LCA and reports

## Index size

The file db_prefix.filter is the main and biggest file for ganon database and it stores the interleaved bloom filter. Its size is defined mainly on the amount of the input reference sequences (`-i`) but also can also be adjusted by a combination of parameters: `--bin-length` `--max-fp` `--kmer-size` `--hash-functions`.

Ganon will try to find the best `--bin-length` given `--max-fp` `--kmer-size` `--hash-functions`. If you have a limited amount of resources available, you can set the `--max-bloom-size` to put an approximate limit on the size of the filter (ps: there is a minimum size necessary to generate the filter given a set of references and chosen parameters). Increasing `--max-fp` will generate smaller filters, but will generate more false positives in the classification step.

`--bin-length` is the size in bp of each group for the taxonomic clustering (with TaxSBP). By default,  `--fragment-length` will be the size of `--bin-length` - `--overlap-length`, meaning that sequences will be split with overlap to fit into the bins. For example: species X has 2 sequences of 120bp each. Considering `--bin-length 50` and `--overlap-length 10` (`--fragment-length 40` consequently) each of the sequences will be split into 50bp and put into a bin with overlap of 10bp, resulting in 3 bins for each sequence (6 in total for species X).

Such adjustment is necessary to equalize the size of each bin, since the IBF requires the individual bloom filters to be of the same size. Building the IBF based on the biggest sequence group in your references will generate the lowest number of bins but a very sparse and gigantic IBF. Building the IBF based on the smallest sequence group in your references will generate the smallest IBF but with too many bins. A balance between those two is necessary to achieve small and fast filters.

## Manual Installation

### Dependencies

#### build

System packages:
- gcc >=7 (check [gcc7 with conda](#installing-gcc7-in-a-separate-environment-with-conda))
- cmake >=3.10

Specific packages:
- Catch2 >=2.7.0 ([d63307](https://github.com/catchorg/Catch2/commit/d63307279412de3870cf97cc6802bae8ab36089e))
- cxxopts >=2.1.2 ([a0de9f](https://github.com/jarro2783/cxxopts/commit/a0de9f3ba1035a3c4f5ffcd960cb94e4e12d40c5))
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
- pylca >= 1.0.0 ([d1474b](https://github.com/pirovc/pylca/commit/d1474b2ec2c028963bafce278ccb69cc21c061fa))
- binpacking >=1.4.1 ([v1.4.1](https://pypi.org/project/binpacking/1.4.1/))

** Please make sure that the system packages are supported/installed in your environment. All other packages are installed in the next steps.

### Obtaining packages

```shh
git clone --recurse-submodules https://github.com/pirovc/ganon.git # ganon, catch2, cxxopts, sdsl-lite, seqan
git clone https://github.com/pirovc/taxsbp.git # taxsbp
git clone https://github.com/pirovc/pylca.git # pylca
```

### Installing 

#### taxsbp

```shh
cd taxsbp
python3 setup.py install
taxsbp -h
```

#### pylca

```shh
cd pylca
python3 setup.py install
python3 -c 'from pylca.pylca import *; unittest.main();'
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
conda create -n gcc7 -c quantstack gcc-7 libgcc-7 cmake>=3.8.2
source activate gcc7
```

If you are getting the following error `ganon-classify: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version 'CXXABI_1.3.11' not found` you have to set the `LD_LIBRARY_PATH`

```shh
export LD_LIBRARY_PATH=/home/user/miniconda3/envs/gcc7/lib/
```
