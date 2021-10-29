from ganon.tax import Tax
from ganon.util import *


def table(cfg):
    #validate input input files
    tre_files = validate_input_files(cfg.tre_files, cfg.input_directory, cfg.input_extension, cfg.quiet)

    print_log("Generating table", cfg.quiet)

    # Reports are parsed with cumulative counts
    reports, total_taxa = parse_reports(tre_files, cfg.rank)
    print_log(" - " + str(len(reports)) + " files parsed", cfg.quiet)
    print_log(" - " + str(total_taxa) + " taxa parsed", cfg.quiet)

    # filter reports
    filtered_total_taxa = filter_reports(reports, cfg)
    if total_taxa - filtered_total_taxa:
        print_log(" - Skipped " + str(total_taxa - filtered_total_taxa) + " taxa with filters", cfg.quiet)

    # top_sample and top_all are mutually exclusive
    if cfg.top_sample:
        top_sample_total_taxa = select_top_sample(reports, cfg.top_sample)
        print_log(" - Skipped " + str(filtered_total_taxa - top_sample_total_taxa) + " taxa (--top-sample "+ str(cfg.top_sample)+")", cfg.quiet)
        filtered_total_taxa = top_sample_total_taxa
    elif cfg.top_all:
        top_all_total_taxa = select_top_all(reports, cfg.top_all)
        print_log(" - Skipped " + str(filtered_total_taxa - top_all_total_taxa) + " taxa (--top-all "+ str(cfg.top_all)+")", cfg.quiet)
        filtered_total_taxa = top_all_total_taxa

    if cfg.min_occurrence or cfg.min_occurrence_percentage:
        if cfg.min_occurrence_percentage:
            cfg.min_occurrence = int(len(reports)*cfg.min_occurrence_percentage)
        min_occurrence_total_taxa = select_occurrence(reports, cfg.min_occurrence)
        if cfg.min_occurrence_percentage:
            print_log(" - Skipped " + str(filtered_total_taxa-min_occurrence_total_taxa) + " taxa (--min-occurrence-percentage " + str(cfg.min_occurrence_percentage)+" [" + str(cfg.min_occurrence)+ " files])", cfg.quiet)
        else:
            print_log(" - Skipped " + str(filtered_total_taxa-min_occurrence_total_taxa) + " taxa (--min-occurrence " + str(cfg.min_occurrence)+ ")", cfg.quiet)
        filtered_total_taxa = min_occurrence_total_taxa

    if not cfg.rank:
        # remove cumulative values from multiple ranks
        # done after filtration
        adjust_counts_ranks(reports)

    if not filtered_total_taxa: 
        print_log(" - No taxa left to report", cfg.quiet)
    else:
        lines, cols = write_tsv(reports, cfg)
        print_log(" - " + str(lines) + "x" + str(cols) + " table saved to " + cfg.output_file, cfg.quiet)

    return True

def parse_reports(tre_files, rank):
    reports = {}
    total_taxa = set()
    for tre_file in tre_files:
        reports[tre_file] = {}
        count, lineage, name, total, unclassified = parse_tre_rank(tre_file, rank)
        total_taxa.update(count.keys())
        reports[tre_file]["label"] = tre_file
        reports[tre_file]["count"] = count
        reports[tre_file]["lineage"] = lineage
        reports[tre_file]["name"] = name
        reports[tre_file]["total"] = total
        reports[tre_file]["unclassified"] = unclassified
        reports[tre_file]["filtered"] = 0
    return reports, len(total_taxa)


def parse_tre_rank(tre_file, selected_rank):
    count = {}
    lineage = {}
    name = {}
    unclassified = 0
    classified = 0
    total = 0

    with open(tre_file, "r") as file:
        for line in file:
            rank, taxid, lin, taxa_name, unique_assign, shared_assign, children_assign, cum_assign, cum_perc = line.rstrip().split("\t")
            if rank == "unclassified":
                unclassified = int(cum_assign)
                continue
            elif rank == "root":
                classified = int(cum_assign)
            elif selected_rank and rank != selected_rank:
                continue

            lineage[taxid] = lin.split("|")
            name[taxid] = taxa_name
            count[taxid] = int(cum_assign)

    total = unclassified + classified
    return count, lineage, name, total, unclassified


def filter_reports(reports, cfg):
    filtered_total_taxa = set()
    for file, rep in reports.items():
        for taxid in list(rep["count"]):
            count = rep["count"][taxid]
            filtered = False
            if (cfg.min_count and count < cfg.min_count) or (cfg.min_percentage and count/rep["total"] < cfg.min_percentage):
                filtered = True
            elif cfg.taxids and not any(t in cfg.taxids for t in rep["lineage"][taxid]): 
                filtered = True
            elif cfg.names and not rep["name"][taxid] in cfg.names:
                filtered = True
            elif cfg.names_with and not any(n in rep["name"][taxid] for n in cfg.names_with):
                filtered = True

            if filtered:
                rep["filtered"] += count
                del rep["count"][taxid]
                del rep["lineage"][taxid]
                del rep["name"][taxid]
            else:
                filtered_total_taxa.add(taxid)

    return len(filtered_total_taxa)


def select_top_sample(reports, top_sample):
    top_sample_total_taxa = set()
    for file, rep in reports.items():
        for i, (taxid, count) in enumerate(sorted(rep["count"].items(), key=lambda x: x[1], reverse=True)): # sorted by count
            if i < top_sample:
                top_sample_total_taxa.add(taxid)
                continue
            rep["filtered"] += count
            del rep["count"][taxid]
            del rep["lineage"][taxid]
            del rep["name"][taxid]

    return len(top_sample_total_taxa)


def select_top_all(reports, top_all):
    total_taxa = set()
    total_counts = get_total_counts(reports)
    top_taxids = []
    for i, taxid in enumerate(sorted(total_counts, key=lambda kv: total_counts[kv]["sum_percentage"], reverse=True)):
        if i < top_all:
            top_taxids.append(taxid)

    for file, rep in reports.items():
        for taxid in list(rep["count"]):
            if taxid in top_taxids:
                total_taxa.add(taxid)
                continue
            rep["filtered"] += rep["count"][taxid]
            del rep["count"][taxid]
            del rep["lineage"][taxid]
            del rep["name"][taxid]

    return len(total_taxa)


def select_occurrence(reports, min_occurrence):
    min_occ_taxids = []
    for taxid, val in get_total_counts(reports).items():
        if val["occurrence"] >= min_occurrence:
            min_occ_taxids.append(taxid)

    min_occurrence_total_taxa = set()
    for file, rep in reports.items():
        for taxid in list(rep["count"]):
            if taxid in min_occ_taxids:
                min_occurrence_total_taxa.add(taxid)
                continue
            rep["filtered"] += rep["count"][taxid]
            del rep["count"][taxid]
            del rep["lineage"][taxid]
            del rep["name"][taxid]

    return len(min_occurrence_total_taxa)


def get_total_counts(reports):
    # return sum of % of all samples
    total_counts = {}
    for d in reports.values():
        for taxid, count in d["count"].items():
            if taxid not in total_counts:
                total_counts[taxid] = {"sum_percentage": 0, "occurrence": 0}
            total_counts[taxid]["sum_percentage"] += count/d["total"]
            total_counts[taxid]["occurrence"] += 1
    return total_counts


def adjust_counts_ranks(reports):
    # If reporting multiple ranks, cumulative counts have to be adjusted
    # On the .tre report, every ranks sums up to the total with repeated counts
    # but in the table output, only single counts should be reported
    # In addition, cumulative counts include children assignments that may not be in the report (e.g. assignment to a rank not listed)
    # Going from leaf to root, all reported counts removed from parent ranks counts, and left over are unique matches for the target + unaccounted
    for file, rep in reports.items():
        # Start from higher ranks
        for t in sorted(rep["lineage"], key=lambda k: len(rep["lineage"][k]), reverse=True):
            # Remove already reported count from all parents in its lineage
            for parent in rep["lineage"][t][:-1]:
                if parent in rep["count"]:
                    rep["count"][parent] -= rep["count"][t]


def write_tsv(reports, cfg):
    total_counts = get_total_counts(reports)

    # Sort by taxid
    sorted_taxids = sorted(total_counts.keys())

    # Generate headers
    if cfg.header == "taxid":
        header = [""] + list(sorted_taxids)
    elif cfg.header == "name":
        # Merge all possible names
        names = {}
        for file in reports:
            names.update(reports[file]["name"])
        header = [""] + [names[taxid] for taxid in sorted_taxids]
    elif cfg.header == "lineage":
        # Merge all possible lineages
        lineages = {}
        for file in reports:
            lineages.update(reports[file]["lineage"])
        header = [""] + ["|".join(lineages[taxid]) for taxid in sorted_taxids]
    if cfg.add_unclassified:
        header.append("unclassified")
    if cfg.add_filtered:
        header.append("filtered")

    # generate output as a list of lists with each file in one line
    out_table = []
    out_table.append(header)
    for file in sorted(reports):
        res = reports[file]
        out_line = [res["label"]]
        for taxid in sorted_taxids:
            if taxid in res["count"]:
                v = res["count"][taxid]
                if cfg.output_value == "percentage":
                    v = v/res["total"]
            else:
                v = 0
            out_line.append(v)

        if cfg.add_unclassified:
            out_line.append(res["unclassified"]/res["total"] if cfg.output_value=="percentage" else res["unclassified"])

        if cfg.add_filtered:
            out_line.append(res["filtered"]/res["total"] if cfg.output_value=="percentage" else res["filtered"])

        if cfg.skip_zeros and (len(out_line) > 1 and max(out_line[1:]) == 0):
            continue

        out_table.append(out_line)

    if len(reports) > len(out_table):
        print_log(" - Skipped " + str(len(reports)-len(out_table)) + " files with only zero counts", cfg.quiet)

    # check only zero samples/taxa
    if cfg.skip_zeros:
        filtered_out_table = []
        # loop on transposed
        for i, line in enumerate(map(list, zip(*out_table))):
            if i == 0:
                continue  # skip header
            if len(line) > 1 and max(line[1:]) > 0:
                filtered_out_table.append(line)
        if len(out_table) > len(filtered_out_table):
            print_log(" - Skipped " + str(len(out_table)-len(filtered_out_table)) + " taxa with only zero counts", cfg.quiet)
        # transpose again
        out_table = map(list, zip(*filtered_out_table))
    else:
        taxa = len(header) - 1

    # "--transpose" table (by default is already transposed)
    if not cfg.transpose:
        out_table = map(list, zip(*out_table))

    # Write file
    out_file = open(cfg.output_file, "w")
    lines = 0
    for line in out_table:
        print(*line, sep="\t" if cfg.output_format=="tsv" else ",", file=out_file)
        lines += 1
    cols = len(line)
    out_file.close()

    return lines-1, cols-1
