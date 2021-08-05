from ganon.util import *
from ganon.report import report
from ganon.config import Config
from ganon.gnn import Gnn

def classify(cfg):
    print_log("Classifying reads (ganon-classify)", cfg.quiet)
    run_ganon_classify = " ".join([cfg.path_exec['classify'],
                                   "--single-reads " +  ",".join(cfg.single_reads) if cfg.single_reads else "",
                                   "--paired-reads " +  ",".join(cfg.paired_reads) if cfg.paired_reads else "",
                                   "--ibf " + ",".join([db_prefix+".ibf" for db_prefix in cfg.db_prefix]),
                                   "--map " + ",".join([db_prefix+".map" for db_prefix in cfg.db_prefix]), 
                                   "--tax " + ",".join([db_prefix+".tax" for db_prefix in cfg.db_prefix]),
                                   "--kmer-size " + ",".join([str(Gnn(file=db_prefix+".gnn").kmer_size) for db_prefix in cfg.db_prefix]),
                                   "--hierarchy-labels " + ",".join(cfg.hierarchy_labels) if cfg.hierarchy_labels else "",
                                   "--max-error " + ",".join([str(me) for me in cfg.max_error]) if cfg.max_error else "",
                                   "--min-kmers " + ",".join([str(mk) for mk in cfg.min_kmers]) if cfg.min_kmers else "",
                                   "--max-error-unique " + ",".join([str(meu) for meu in cfg.max_error_unique]) if cfg.max_error_unique else "",
                                   "--strata-filter " + ",".join([str(sf) for sf in cfg.strata_filter]) if cfg.strata_filter else "",
                                   "--offset " + str(cfg.offset) if cfg.offset else "",
                                   "--output-prefix " + cfg.output_prefix if cfg.output_prefix else "",
                                   "--output-lca",
                                   "--output-all" if cfg.output_all else "",
                                   "--output-unclassified" if cfg.output_unclassified else "",
                                   "--output-single" if cfg.output_single else "",
                                   "--threads " + str(cfg.threads) if cfg.threads else "",
                                   "--n-reads " + str(cfg.n_reads) if cfg.n_reads is not None else "",
                                   "--n-batches " + str(cfg.n_batches) if cfg.n_batches is not None else "",
                                   "--verbose" if cfg.verbose else "",
                                   "--quiet" if cfg.quiet else ""])
    #print(run_ganon_classify)
    stdout, stderr = run(run_ganon_classify)
    if not cfg.output_prefix:
        print(stdout)
    print_log(stderr, cfg.quiet)

    if cfg.output_prefix:
        report_params = {"db_prefix": cfg.db_prefix,
                         "rep_file": cfg.output_prefix+".rep",
                         "output_prefix": cfg.output_prefix,
                         "ranks": cfg.ranks,
                         "output_format": "tsv",
                         "verbose": cfg.verbose,
                         "quiet": cfg.quiet}
        report_cfg = Config("report", **report_params)
        ret = report(report_cfg)
        return ret
    else:
        return True
