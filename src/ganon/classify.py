from ganon.util import run, print_log, check_file
from ganon.report import report
from ganon.reassign import reassign
from ganon.config import Config


def classify(cfg):

    # Set paths
    if not cfg.set_paths():
        return False

    print_log("Classifying reads", cfg.quiet)

    filter_files = []
    tax_files = []
    for db_prefix in cfg.db_prefix:
        if check_file(db_prefix + ".hibf"):
            filter_files.append(db_prefix + ".hibf")
        elif check_file(db_prefix + ".ibf"):
            filter_files.append(db_prefix + ".ibf")

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
                                   "--fpr-query " + ",".join([str(fq) for fq in cfg.fpr_query]) if cfg.fpr_query else "",
                                   "--output-prefix " + cfg.output_prefix if cfg.output_prefix else "",
                                   "--skip-lca" if cfg.multiple_matches != "lca" else "",
                                   "--output-lca" if cfg.multiple_matches == "lca" and cfg.output_one else "",
                                   "--output-all" if cfg.output_all or cfg.multiple_matches == "em" else "",
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
    else:
        if cfg.multiple_matches == "em":
            reassign_params = {"input_prefix": cfg.output_prefix,
                               "output_prefix": cfg.output_prefix,
                               "remove_all": False if cfg.output_all else True,
                               "skip_one": False if cfg.output_one else True,
                               "verbose": cfg.verbose,
                               "quiet": cfg.quiet}
            reassign_cfg = Config("reassign", **reassign_params)
            print_log("- - - - - - - - - -", cfg.quiet)
            ret = reassign(reassign_cfg)
            if not ret:
                return False

        if tax_files and not cfg.skip_report:
            report_params = {"db_prefix": cfg.db_prefix,
                             "input": cfg.output_prefix + ".rep",
                             "output_prefix": cfg.output_prefix,
                             "min_count": cfg.min_count,
                             "ranks": cfg.ranks,
                             "output_format": "tsv",
                             "verbose": cfg.verbose,
                             "report_type": cfg.report_type,
                             "quiet": cfg.quiet}
            report_cfg = Config("report", **report_params)
            print_log("- - - - - - - - - -", cfg.quiet)
            ret = report(report_cfg)
            if not ret:
                return False

    return True
