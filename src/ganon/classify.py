from ganon.util import *
from ganon.report import report
from ganon.config import Config
from ganon.gnn import Gnn


def classify(cfg):
    print_log("Classifying reads (ganon-classify)", cfg.quiet)

    kmer_size = set()
    window_size = set()
    for i, db_prefix in enumerate(cfg.db_prefix):
        gnn = Gnn(file=db_prefix+".gnn")
        if cfg.hierarchy_labels:
            # one per hierarchical label
            kmer_size.add((cfg.hierarchy_labels[i], gnn.kmer_size))
            window_size.add((cfg.hierarchy_labels[i], gnn.window_size))
        else:
            kmer_size.add(gnn.kmer_size)
            window_size.add(gnn.window_size)

    if cfg.hierarchy_labels:
        # Sort by hierarchical label and get value
        kmer_size = [value for _, value in sorted(kmer_size, key=lambda tup: tup[0])]
        window_size = [value for _, value in sorted(window_size, key=lambda tup: tup[0])]

    if len(kmer_size) != len(window_size):
        print_log("Incompatible databases", cfg.quiet)
        return False

    run_ganon_classify = " ".join([cfg.path_exec['classify'],
                                   "--single-reads " + ",".join(cfg.single_reads) if cfg.single_reads else "",
                                   "--paired-reads " + ",".join(cfg.paired_reads) if cfg.paired_reads else "",
                                   "--ibf " + ",".join([db_prefix+".ibf" for db_prefix in cfg.db_prefix]),
                                   "--map " + ",".join([db_prefix+".map" for db_prefix in cfg.db_prefix]),
                                   "--tax " + ",".join([db_prefix+".tax" for db_prefix in cfg.db_prefix]),
                                   "--kmer-size " + ",".join(map(str, kmer_size)),
                                   "--window-size " + ",".join(map(str, window_size)),
                                   "--hierarchy-labels " + ",".join(cfg.hierarchy_labels) if cfg.hierarchy_labels else "",
                                   "--abs-cutoff " + ",".join([str(ac) for ac in cfg.abs_cutoff]) if cfg.abs_cutoff else "",
                                   "--rel-cutoff " + ",".join([str(rc) for rc in cfg.rel_cutoff]) if cfg.rel_cutoff else "",
                                   "--abs-filter " + ",".join([str(af) for af in cfg.abs_filter]) if cfg.abs_filter else "",
                                   "--rel-filter " + ",".join([str(rf) for rf in cfg.rel_filter]) if cfg.rel_filter else "",
                                   "--offset " + str(cfg.offset) if cfg.offset else "",
                                   "--output-prefix " + cfg.output_prefix if cfg.output_prefix else "",
                                   "--output-lca" if cfg.output_lca else "",
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
