# Parameters

```
usage: ganon [-h] [-v]
             {build,build-custom,update,classify,reassign,report,table} ...

- - - - - - - - - -
   _  _  _  _  _   
  (_|(_|| |(_)| |  
   _|   v. 2.4.2.dev9+g3c2298a32
- - - - - - - - - -

positional arguments:
  {build,build-custom,update,classify,reassign,report,table}
    build               Download and build ganon default databases
                        (refseq/genbank)
    build-custom        Build custom ganon databases
    update              Update ganon default databases
    classify            Classify reads against built databases
    reassign            Reassign reads with multiple matches with an EM
                        algorithm
    report              Generate reports from classification results
    table               Generate table from reports

options:
  -h, --help            show this help message and exit
  -v, --version         Show program's version number and exit.
```

<details>
  <summary>ganon build</summary>

```
usage: ganon build [-h] [-g [ ...]] [-a [ ...]] [-l ] [-x ] [-m [ ...]] [-b [ ...]] [-o ] [-c] [-r] [-e] [-u ]
                   [-z [ ...]] [--skip-genome-size] [--download-threads ] -d DB_PREFIX [-t ] [-p ] [-k ] [-w ] [-s ]
                   [-f ] [-j ] [-y ] [-v ] [--restart] [--verbose] [--quiet] [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -g, --organism-group [ ...]
                        One or more organism groups to download [archaea, bacteria, fungi, human, invertebrate,
                        metagenomes, other, plant, protozoa, vertebrate_mammalian, vertebrate_other, viral]. Mutually
                        exclusive --taxid (default: None)
  -a, --taxid [ ...]    One or more taxonomic identifiers to download. e.g. 562 (-x ncbi) or 's__Escherichia coli' (-x
                        gtdb). Mutually exclusive --organism-group (default: None)
  -d, --db-prefix DB_PREFIX
                        Database output prefix

database arguments:
  -l, --level           Highest level to build the database. Options: any available taxonomic rank [species, genus,
                        ...], 'leaves' for taxonomic leaves or 'assembly' for a assembly/strain based analysis (default:
                        species)

taxonomy arguments:
  -x, --taxonomy        Use taxonomy to enable taxonomic classification, lca and tax. reports [ncbi, gtdb, skip]
                        (default: ncbi)
  -m, --taxonomy-files [ ...]
                        Use local taxonomy files instead of downloading. For ncbi: taxdump.tar.gz OR nodes.dmp
                        [names.dmp merged.dmp]. For gtdb: *taxonomy.tsv.gz (default: None)

download arguments:
  -b, --source [ ...]   Source to download [refseq, genbank] (default: ['refseq'])
  -o, --top             Download limited assemblies for each taxa. 0 for all. (default: 0)
  -c, --complete-genomes
                        Download only sub-set of complete genomes. Mutually exclusive --complete-and-reference-genomes
                        (default: False)
  -r, --reference-genomes
                        Download only sub-set of reference genomes. Mutually exclusive --complete-and-reference-genomes
                        (default: False)
  -e, --complete-and-reference-genomes
                        Download union of complete and reference genomes sub-set. Mutually exclusive --complete-
                        genomes/--reference-genomes (default: False)
  -u, --genome-updater 
                        Additional genome_updater parameters (https://github.com/pirovc/genome_updater) (default: None)
  -z, --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)
  --skip-genome-size    Do not attempt to get genome sizes. Activate this option when using sequences not representing
                        full genomes. (default: False)
  --download-threads    Number of parallel sequence downloads from NCBI. (default: 8)

general arguments:
  -t, --threads         Number of sub-processes/threads to use (default: 1)
  -p, --max-fp          Max. false positive for bloom filters. Mutually exclusive --filter-size. Defaults to 0.001 with
                        --filter-type hibf or 0.05 with --filter-type ibf. (default: None)
  -k, --kmer-size       The k-mer size to split sequences. (default: 19)
  -w, --window-size     The window-size to build filter with minimizers. (default: 31)
  -s, --hash-functions 
                        The number of hash functions for the interleaved bloom filter [1-5]. With --filter-type ibf, 0
                        will try to set optimal value. (default: 4)
  -f, --filter-size     Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. Only valid for --filter-
                        type ibf. (default: 0)
  -j, --mode            Create smaller or faster filters at the cost of classification speed or database size,
                        respectively [avg, smaller, smallest, faster, fastest]. If --filter-size is used,
                        smaller/smallest refers to the false positive. By default, an average value is calculated to
                        balance classification speed and database size. Only valid for --filter-type ibf. (default: avg)
  -y, --min-length      Skip sequences smaller then value defined. 0 to not skip any sequence. Only valid for --filter-
                        type ibf. (default: 0)
  -v, --filter-type     Variant of bloom filter to use [hibf, ibf]. hibf requires raptor >= v3.0.1 installed or binary
                        path set with --raptor-path. --mode, --filter-size and --min-length will be ignored with hibf.
                        hibf will set --max-fp 0.001 as default. (default: hibf)

optional arguments:
  --restart             Restart build/update from scratch, do not try to resume from the latest possible step.
                        {db_prefix}_files/ will be deleted if present. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
  --write-info-file     Save copy of target info generated to {db_prefix}.info.tsv. Can be re-used as --input-file for
                        further attempts. (default: False)
```

</details>

<details>
  <summary>ganon build-custom</summary>

```
usage: ganon build-custom [-h] [-i [ ...]] [-e ] [-c] [-n ] [-a ] [-l ] [-z [ ...]] [--skip-genome-size] [-x ] [-b ]
                          [-m [ ...]] [-u [ ...]] [-g [ ...]] [--keep-invalid-taxa] [-r [ ...]] [-q [ ...]] -d DB_PREFIX
                          [-t ] [-p ] [-k ] [-w ] [-s ] [-f ] [-j ] [-y ] [-v ] [--restart] [--verbose] [--quiet]
                          [--write-info-file]

options:
  -h, --help            show this help message and exit

required arguments:
  -i, --input [ ...]    Input file(s) and/or folder(s). Mutually exclusive --input-file. (default: None)
  -e, --input-extension 
                        Required if --input contains folder(s). Wildcards/Shell Expansions not supported (e.g. *).
                        (default: fna.gz)
  -c, --input-recursive
                        Look for files recursively in folder(s) provided with --input (default: False)
  -d, --db-prefix DB_PREFIX
                        Database output prefix

custom arguments:
  -n, --input-file      Tab-separated file with all necessary file/sequence information. Fields: file [<tab> target
                        <tab> node <tab> specialization <tab> specialization name]. For details:
                        https://pirovc.github.io/ganon/custom_databases/. Mutually exclusive --input (default: None)
  -a, --input-target    Target to use [file, sequence]. Parse input by file or by sequence. Using 'file' is recommended
                        and will speed-up the building process (default: file)
  -l, --level           Max. level to build the database. By default, --level is the --input-target. Options: any
                        available taxonomic rank [species, genus, ...] or 'leaves' (requires --taxonomy). Further
                        specialization options [assembly, custom]. assembly will retrieve and use the assembly accession
                        and name. custom requires and uses the specialization field in the --input-file. (default: None)
  -z, --genome-size-files [ ...]
                        Specific files for genome size estimation - otherwise files will be downloaded (default: None)
  --skip-genome-size    Do not attempt to get genome sizes. Activate this option when using sequences not representing
                        full genomes. (default: False)

taxonomy arguments:
  -x, --taxonomy        Taxonomy matching the --input/--input-file. Enables taxonomic classification, lca and tax.
                        reports [ncbi, gtdb, gtdb-80, gtdb-83, gtdb-86.2, gtdb-89, gtdb-95, gtdb-202, gtdb-207,
                        gtdb-214.1, gtdb-220, gtdb-226, gtdb-232, skip] (default: ncbi)
  -b, --convert-taxonomy 
                        Convert input taxonomy nodes (--taxonomy) to [ncbi-latest, gtdb-80, gtdb-83, gtdb-86.2, gtdb-89,
                        gtdb-95, gtdb-202, gtdb-207, gtdb-214.1, gtdb-220, gtdb-226, gtdb-232]. (default: None)
  -m, --taxonomy-files [ ...]
                        Use local taxonomy files instead of downloading. For ncbi: taxdump.tar.gz OR nodes.dmp
                        [names.dmp merged.dmp]. For gtdb: *taxonomy.tsv.gz (default: None)
  -u, --convert-taxonomy-files [ ...]
                        Use local taxonomy files instead of downloading. For ncbi-latest: taxdump.tar.gz OR nodes.dmp
                        [names.dmp merged.dmp]. For gtdb-version: *taxonomy.tsv.gz (default: None)
  -g, --convert-gtdb-files [ ...]
                        Use local gtdb conversion files instead of downloading. One for each version used in --taxonomy
                        and --convert-taxonomy. Files from https://github.com/pirovc/multitax/tree/main/data/gtdb
                        (default: None)
  --keep-invalid-taxa   Keep invalid taxa in the database, will be assigned to the root of the taxonomic tree. (default:
                        False)

ncbi arguments:
  -r, --ncbi-sequence-info [ ...]
                        Uses NCBI e-utils webservices or downloads accession2taxid files to extract target information.
                        [eutils, nucl_gb, nucl_wgs, nucl_est, nucl_gss, pdb, prot, dead_nucl, dead_wgs, dead_prot or one
                        or more accession2taxid files from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/].
                        By default uses e-utils up-to 50000 sequences or downloads nucl_gb nucl_wgs otherwise. (default:
                        [])
  -q, --ncbi-file-info [ ...]
                        Downloads assembly_summary files to extract target information. [refseq, genbank,
                        refseq_historical, genbank_historical or one or more assembly_summary files from
                        https://ftp.ncbi.nlm.nih.gov/genomes/] (default: ['refseq', 'genbank'])

general arguments:
  -t, --threads         Number of sub-processes/threads to use (default: 1)
  -p, --max-fp          Max. false positive for bloom filters. Mutually exclusive --filter-size. Defaults to 0.001 with
                        --filter-type hibf or 0.05 with --filter-type ibf. (default: None)
  -k, --kmer-size       The k-mer size to split sequences. (default: 19)
  -w, --window-size     The window-size to build filter with minimizers. (default: 31)
  -s, --hash-functions 
                        The number of hash functions for the interleaved bloom filter [1-5]. With --filter-type ibf, 0
                        will try to set optimal value. (default: 4)
  -f, --filter-size     Fixed size for filter in Megabytes (MB). Mutually exclusive --max-fp. Only valid for --filter-
                        type ibf. (default: 0)
  -j, --mode            Create smaller or faster filters at the cost of classification speed or database size,
                        respectively [avg, smaller, smallest, faster, fastest]. If --filter-size is used,
                        smaller/smallest refers to the false positive. By default, an average value is calculated to
                        balance classification speed and database size. Only valid for --filter-type ibf. (default: avg)
  -y, --min-length      Skip sequences smaller then value defined. 0 to not skip any sequence. Only valid for --filter-
                        type ibf. (default: 0)
  -v, --filter-type     Variant of bloom filter to use [hibf, ibf]. hibf requires raptor >= v3.0.1 installed or binary
                        path set with --raptor-path. --mode, --filter-size and --min-length will be ignored with hibf.
                        hibf will set --max-fp 0.001 as default. (default: hibf)

optional arguments:
  --restart             Restart build/update from scratch, do not try to resume from the latest possible step.
                        {db_prefix}_files/ will be deleted if present. (default: False)
  --verbose             Verbose output mode (default: False)
  --quiet               Quiet output mode (default: False)
  --write-info-file     Save copy of target info generated to {db_prefix}.info.tsv. Can be re-used as --input-file for
                        further attempts. (default: False)
```

</details>

<details>
  <summary>ganon update</summary>

```
