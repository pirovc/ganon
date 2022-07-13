from ganon.util import validate_input_files
from ganon.util import print_log


def table(cfg):
    #validate input input files
    tre_files = validate_input_files(cfg.input, cfg.input_extension, cfg.quiet)

    print_log("Generating table", cfg.quiet)

    # Reports are parsed with cumulative counts
    reports, total_taxa = parse_reports(tre_files, cfg.rank)
    root_node = set()
    for rep in reports.values():
        root_node.add(rep["root_node"])
    if len(root_node) > 1:
        print_log("ERROR: input files should share the same root node, but " + str(len(root_node)) +
                  " root nodes were found: " ",".join(root_node), cfg.quiet)
    else:
        root_node = root_node.pop()

    print_log(" - " + str(len(reports)) + " files parsed", cfg.quiet)
    print_log(" - " + str(total_taxa) + " taxa parsed", cfg.quiet)

    # filter reports
    filtered_total_taxa = filter_reports(reports, cfg, root_node)
    if total_taxa - filtered_total_taxa:
        print_log(" - Skipped " + str(total_taxa - filtered_total_taxa) + " taxa with filters", cfg.quiet)

    # top_sample and top_all are mutually exclusive
    if cfg.top_sample:
        top_sample_total_taxa = select_top_sample(reports, cfg.top_sample, root_node)
        print_log(" - Skipped " + str(filtered_total_taxa - top_sample_total_taxa) +
                  " taxa (--top-sample " + str(cfg.top_sample)+")", cfg.quiet)
        filtered_total_taxa = top_sample_total_taxa
    elif cfg.top_all:
        top_all_total_taxa = select_top_all(reports, cfg.top_all, root_node)
        print_log(" - Skipped " + str(filtered_total_taxa - top_all_total_taxa) +
                  " taxa (--top-all " + str(cfg.top_all)+")", cfg.quiet)
        filtered_total_taxa = top_all_total_taxa

    if cfg.min_frequency:
        if cfg.min_frequency < 1:
            mf = int(len(reports)*cfg.min_frequency)
        else:
            mf = cfg.min_frequency
        min_frequency_total_taxa = select_frequency(reports, mf)
        print_log(" - Skipped " + str(filtered_total_taxa - min_frequency_total_taxa) +
                  " taxa (--min-frequency " + str(cfg.min_frequency) + ")", cfg.quiet)
        filtered_total_taxa = min_frequency_total_taxa

    if not cfg.rank:
        # remove cumulative values from multiple ranks
        # step should be at the end after filters were applied on the cumulative counts
        # if requested, counts directed to root are remove to unclassified
        adjust_counts_ranks(reports, cfg.no_root, root_node)

    # remove root from lineages
    if cfg.no_root:
        for file, rep in reports.items():
            for t in rep["count"]:
                del rep["lineage"][t][0]

    if not filtered_total_taxa:
        print_log(" - No taxa left to report", cfg.quiet)
    else:
        # build table in a list of lists
        out_table = build_table(reports, cfg)

        # "--skip-zeros" trim table on lines and cols
        if cfg.skip_zeros:
            l = len(out_table)
            # trim rows and cols
            out_table = trim_table(out_table)
            if len(out_table) < l:
                print_log(" - Skipped " + str(l-len(out_table)) + " files with only zero counts", cfg.quiet)
            c = len(out_table[0])
            out_table = transpose(trim_table(transpose(out_table)))
            if len(out_table[0]) < c:
                print_log(" - Skipped " + str(c-len(out_table[0])) + " taxa with only zero counts", cfg.quiet)

        # "--transpose" table (by default is already transposed)
        if not cfg.transpose:
            out_table = transpose(out_table)

        lines, cols = write_tsv(out_table, cfg.output_file, cfg.output_format)
        print_log(" - " + str(lines) + "x" + str(cols) + " table saved to " + cfg.output_file, cfg.quiet)

    return True


def parse_reports(tre_files, rank):
    reports = {}
    total_taxa = set()
    for tre_file in tre_files:
        reports[tre_file] = {}
        count, lineage, name, total, unclassified, root_node = parse_tre_rank(tre_file, rank)
        total_taxa.update(count.keys())
        reports[tre_file]["label"] = tre_file
        reports[tre_file]["count"] = count
        reports[tre_file]["lineage"] = lineage
        reports[tre_file]["name"] = name
        reports[tre_file]["total"] = total
        reports[tre_file]["unclassified"] = unclassified
        reports[tre_file]["filtered"] = 0
        reports[tre_file]["root_node"] = root_node
    return reports, len(total_taxa)


def parse_tre_rank(tre_file, selected_rank):
    count = {}
    lineage = {}
    name = {}
    unclassified = 0
    classified = 0
    total = 0
    root_node = "1"

    with open(tre_file, "r") as file:
        for line in file:
            rank, taxid, lin, taxa_name, unique_assign, shared_assign, children_assign, cum_assign, cum_perc = line.rstrip().split("\t")
            if rank == "unclassified":
                unclassified = int(cum_assign)
                continue
            elif rank == "root":
                classified = int(cum_assign)
                root_node = taxid
                if selected_rank:
                    continue  # do not include root to the report when using single rank
            elif selected_rank and rank != selected_rank:
                continue

            lineage[taxid] = lin.split("|")
            name[taxid] = taxa_name
            count[taxid] = int(cum_assign)

    total = unclassified + classified
    return count, lineage, name, total, unclassified, root_node


def filter_reports(reports, cfg, root_node):
    filtered_total_taxa = set()
    for file, rep in reports.items():
        for taxid in list(rep["count"]):
            count = rep["count"][taxid]
            filtered = False
            if cfg.min_count:
                if cfg.min_count > 1 and count < cfg.min_count:
                    filtered = True
                elif cfg.min_count < 1 and (count/rep["total"]) < cfg.min_count:
                    filtered = True

            if cfg.max_count:
                if cfg.max_count > 1 and count > cfg.max_count:
                    filtered = True
                elif cfg.max_count < 1 and (count/rep["total"]) > cfg.max_count:
                    filtered = True

            if cfg.taxids and not any(t in cfg.taxids for t in rep["lineage"][taxid]):
                filtered = True
            elif cfg.names and not rep["name"][taxid] in cfg.names:
                filtered = True
            elif cfg.names_with and not any(n in rep["name"][taxid] for n in cfg.names_with):
                filtered = True

            # do not filter root node
            if filtered and taxid != root_node:
                rep["filtered"] += count
                del rep["count"][taxid]
                del rep["lineage"][taxid]
                del rep["name"][taxid]
            else:
                filtered_total_taxa.add(taxid)

    return len(filtered_total_taxa)


def select_top_sample(reports, top_sample, root_node):
    top_sample_total_taxa = set(root_node)  # always keep root
    for file, rep in reports.items():
        i = 0
        for taxid, count in sorted(rep["count"].items(), key=lambda x: x[1], reverse=True):  # sorted by count
            if taxid == root_node:  # do not count root as an top entry
                continue
            if i < top_sample:
                top_sample_total_taxa.add(taxid)
                i += 1
                continue
            rep["filtered"] += count
            del rep["count"][taxid]
            del rep["lineage"][taxid]
            del rep["name"][taxid]

    return len(top_sample_total_taxa)


def select_top_all(reports, top_all, root_node):
    total_taxa = set()
    total_counts = get_total_counts(reports)
    top_taxids = set(root_node)  # always keep root
    i = 0
    for taxid in sorted(total_counts, key=lambda kv: total_counts[kv]["sum_percentage"], reverse=True):
        if taxid == root_node:  # do not count root as an top entry
            continue
        elif i < top_all:
            top_taxids.add(taxid)
            i += 1

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


def select_frequency(reports, min_frequency):
    min_occ_taxids = []
    for taxid, val in get_total_counts(reports).items():
        if val["frequency"] >= min_frequency:
            min_occ_taxids.append(taxid)

    min_frequency_total_taxa = set()
    for file, rep in reports.items():
        for taxid in list(rep["count"]):
            if taxid in min_occ_taxids:
                min_frequency_total_taxa.add(taxid)
                continue
            rep["filtered"] += rep["count"][taxid]
            del rep["count"][taxid]
            del rep["lineage"][taxid]
            del rep["name"][taxid]

    return len(min_frequency_total_taxa)


def get_total_counts(reports):
    # return sum of % of all samples
    total_counts = {}
    for d in reports.values():
        for taxid, count in d["count"].items():
            if taxid not in total_counts:
                total_counts[taxid] = {"sum_percentage": 0, "frequency": 0}
            total_counts[taxid]["sum_percentage"] += count/d["total"]
            total_counts[taxid]["frequency"] += 1
    return total_counts


def adjust_counts_ranks(reports, no_root, root_node):
    # If reporting multiple ranks, cumulative counts have to be adjusted
    # On the .tre report, every ranks sums up to the total with repeated counts
    # but in the table output, only single counts should be reported
    # In addition, cumulative counts include children assignments that may not
    # be in the report (e.g. assignment to a rank not listed)
    # Going from leaf to root, all reported counts removed from parent ranks counts,
    # and left over are unique matches for the target + unaccounted
    for file, rep in reports.items():
        # Start from higher ranks
        for t in sorted(rep["lineage"], key=lambda k: len(rep["lineage"][k]), reverse=True):
            # Remove already reported count from all parents in its lineage
            for parent in rep["lineage"][t][:-1]:
                if parent in rep["count"]:
                    rep["count"][parent] -= rep["count"][t]

        # Move left over counts at root to unclassified
        if no_root:
            rep["unclassified"] += rep["count"][root_node]
            del rep["count"][root_node]
            del rep["lineage"][root_node]
            del rep["name"][root_node]


def build_table(reports, cfg):
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
    if cfg.unclassified_label:
        header.append(cfg.unclassified_label)
    if cfg.filtered_label and cfg.filtered_label != cfg.unclassified_label:
        header.append(cfg.filtered_label)

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

        # Add unclassified/filtered at the end in the according labels
        if cfg.unclassified_label:
            unc = res["unclassified"]/res["total"] if cfg.output_value == "percentage" else res["unclassified"]
            if cfg.unclassified_label != cfg.filtered_label:
                out_line.append(unc)
        if cfg.filtered_label:
            fil = res["filtered"]/res["total"] if cfg.output_value == "percentage" else res["filtered"]
            if cfg.filtered_label == cfg.unclassified_label:
                out_line.append(unc+fil)
            else:
                out_line.append(fil)

        out_table.append(out_line)
    return out_table


def write_tsv(out_table, output_file, output_format):
    # Write file
    out_file = open(output_file, "w")
    lines = 0
    for line in out_table:
        print(*line, sep="\t" if output_format == "tsv" else ",", file=out_file)
        lines += 1
    cols = len(line)
    out_file.close()

    return lines-1, cols-1


def trim_table(table):
    trimmed_table = [table[0]]  # get header
    for line in table[1:]:
        values = line[1:]  # skip header on each line
        if max(values) > 0:
            trimmed_table.append(line)  # keep non-zero line
    return trimmed_table


def transpose(table):
    return list(map(list, zip(*table)))
