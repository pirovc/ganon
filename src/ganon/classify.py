from ganon.util import run, print_log
from ganon.report import report
from ganon.config import Config


def classify(cfg):
    print_log("Classifying reads (ganon-classify)", cfg.quiet)

    run_ganon_classify = " ".join([cfg.path_exec['classify'],
                                   "--single-reads " + ",".join(cfg.single_reads) if cfg.single_reads else "",
                                   "--paired-reads " + ",".join(cfg.paired_reads) if cfg.paired_reads else "",
                                   "--ibf " + ",".join([db_prefix+".ibf" for db_prefix in cfg.db_prefix]),
                                   "--tax " + ",".join([db_prefix+".tax" for db_prefix in cfg.db_prefix]),
                                   "--hierarchy-labels " + ",".join(cfg.hierarchy_labels) if cfg.hierarchy_labels else "",
                                   "--rel-cutoff " + ",".join([str(rc) for rc in cfg.rel_cutoff]) if cfg.rel_cutoff else "",
                                   "--rel-filter " + ",".join([str(rf) for rf in cfg.rel_filter]) if cfg.rel_filter else "",
                                   "--output-prefix " + cfg.output_prefix if cfg.output_prefix else "",
                                   "--output-lca" if cfg.output_lca else "",
                                   "--output-all" if cfg.output_all else "",
                                   "--output-unclassified" if cfg.output_unclassified else "",
                                   "--output-single" if cfg.output_single else "",
                                   "--threads " + str(cfg.threads) if cfg.threads else "",
                                   "--n-reads " + str(cfg.n_reads) if cfg.n_reads is not None else "",
                                   "--n-batches " + str(cfg.n_batches) if cfg.n_batches is not None else "",
                                   "--verbose" if cfg.verbose else "",
                                   "--hibf" if cfg.hibf else "",
                                   "--quiet" if cfg.quiet else ""])
    stdout = run(run_ganon_classify, ret_stdout=True, quiet=cfg.quiet)

    if not cfg.output_prefix:
        print(stdout)

    if cfg.output_prefix:
        report_params = {"db_prefix": cfg.db_prefix,
                         "input": cfg.output_prefix+".rep",
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
