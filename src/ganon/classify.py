from ganon.util import run, print_log, check_file
from ganon.report import report
from ganon.reassign import reassign
from ganon.config import Config


def classify(cfg):

    # Set paths
    if not cfg.set_paths():
        return False

    print_log("Classifying reads (ganon-classify)", cfg.quiet)

    filter_files = []
    tax_files = []
    hibf = False
    for db_prefix in cfg.db_prefix:
        if check_file(db_prefix + ".hibf"):
            filter_files.append(db_prefix + ".hibf")
            hibf = True
        elif check_file(db_prefix + ".ibf"):
            filter_files.append(db_prefix + ".ibf")

        # No need to send .tax when reassign is requested
        if check_file(db_prefix + ".tax"):
            tax_files.append(db_prefix + ".tax")

    # Just use if all dbs have tax
    tax_files = ",".join(tax_files) if len(tax_files) == len(filter_files) else ""
    filter_files = ",".join(filter_files)

    run_ganon_classify = " ".join([cfg.path_exec['classify'],
                                   "--single-reads " + ",".join(cfg.single_reads) if cfg.single_reads else "",
                                   "--paired-reads " + ",".join(cfg.paired_reads) if cfg.paired_reads else "",
                                   "--ibf " + filter_files,
                                   "--tax " + tax_files if tax_files else "",
                                   "--hierarchy-labels " + ",".join(cfg.hierarchy_labels) if cfg.hierarchy_labels else "",
                                   "--rel-cutoff " + ",".join([str(rc) for rc in cfg.rel_cutoff]) if cfg.rel_cutoff else "",
                                   "--rel-filter " + ",".join([str(rf) for rf in cfg.rel_filter]) if cfg.rel_filter else "",
                                   "--output-prefix " + cfg.output_prefix if cfg.output_prefix else "",
                                   "--output-lca" if cfg.output_lca else "",
                                   "--output-all" if cfg.output_all or cfg.reassign else "",
                                   "--output-unclassified" if cfg.output_unclassified else "",
                                   "--output-single" if cfg.output_single else "",
                                   "--threads " + str(cfg.threads) if cfg.threads else "",
                                   "--n-reads " + str(cfg.n_reads) if cfg.n_reads is not None else "",
                                   "--n-batches " + str(cfg.n_batches) if cfg.n_batches is not None else "",
                                   "--verbose" if cfg.verbose else "",
                                   "--hibf" if hibf else "",
                                   "--quiet" if cfg.quiet else ""])
    stdout = run(run_ganon_classify, ret_stdout=True, quiet=cfg.quiet)

    if not cfg.output_prefix:
        print(stdout)

    report_input = cfg.output_prefix + ".rep"
    report_output = cfg.output_prefix
    if cfg.reassign:
        reassign_params = {"input_prefix": cfg.output_prefix,
                           "output_prefix": cfg.output_prefix,
                           "verbose": cfg.verbose,
                           "quiet": cfg.quiet}
        reassign_cfg = Config("reassign", **reassign_params)
        ret = reassign(reassign_cfg)
        if ret==False:
            return ret

    if cfg.output_prefix and tax_files:
        report_params = {"db_prefix": cfg.db_prefix,
                         "input": cfg.output_prefix + ".rep",
                         "output_prefix": cfg.output_prefix,
                         "min_count": 0.0001,
                         "ranks": cfg.ranks,
                         "output_format": "tsv",
                         "verbose": cfg.verbose,
                         "quiet": cfg.quiet}
        report_cfg = Config("report", **report_params)
        ret = report(report_cfg)
        return ret
    else:
        return True
