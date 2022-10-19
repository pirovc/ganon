import os
from math import floor

from ganon.util import validate_input_files
from ganon.util import print_log
from ganon.tax_util import get_genome_size

from multitax import CustomTx, NcbiTx, GtdbTx


def report(cfg):
    #validate input input files
    rep_files = validate_input_files(cfg.input, cfg.input_extension, cfg.quiet)

    # Parse taxonomy or download new
    tax_args = {"undefined_node": "",
                "undefined_rank": "na",
                "undefined_name": "na",
                "root_name": "root",
                "root_rank": "root"}

    genome_sizes = {}
    if cfg.db_prefix:
        dbp = []
        for prefix in cfg.db_prefix:
            if prefix.endswith(".tax"):
                dbp.append(prefix)
            else:
                dbp.append(prefix+".tax")

        tax = CustomTx(files=dbp,
                       cols=["node", "parent", "rank", "name"],
                       **tax_args)

        for f in dbp:
            with open(f, "r") as file:
                for line in file:
                    node, _, _, _, gsize = line.rstrip().split("\t")
                    # keep largest genome size
                    gsize = int(gsize)
                    if node in genome_sizes and genome_sizes[node] > gsize:
                        continue
                    genome_sizes[node] = int(gsize)

    else:
        if cfg.taxonomy_files:
            print_log("Parsing " + cfg.taxonomy + " taxonomy", cfg.quiet)
        else:
            print_log("Downloading and parsing " + cfg.taxonomy + " taxonomy", cfg.quiet)

        if cfg.taxonomy == "ncbi":
            tax = NcbiTx(files=cfg.taxonomy_files, **tax_args)
        elif cfg.taxonomy == "gtdb":
            tax = GtdbTx(files=cfg.taxonomy_files, **tax_args)

        # TODO get_genome_size -- how to define used nodes?


    # define fixed_ranks or leave it empty for all
    if cfg.ranks and cfg.ranks[0] == "all":
        fixed_ranks = []
    else:
        if not cfg.ranks or cfg.ranks == [""]:
            fixed_ranks = [tax.name(tax.root_node),
                           "superkingdom",
                           "phylum",
                           "class",
                           "order",
                           "family",
                           "genus",
                           "species",
                           "assembly"]
        else:
            fixed_ranks = [tax.name(tax.root_node)] + cfg.ranks

        # If reporting orphan nodes, add at the end
        if not cfg.no_orphan:
            fixed_ranks.append(tax.undefined_rank)

    any_rep = False
    print_log("Generating report(s)", cfg.quiet)
    # Parse report file
    for rep_file in rep_files:
        print_log("", cfg.quiet)

        reports, counts = parse_rep(rep_file)
        if not reports:
            print_log(" - nothing to report for " + rep_file, cfg.quiet)
            continue

        # If skipping/keeping hiearchies, remove all assignments from reports
        if cfg.skip_hierarchy or cfg.keep_hierarchy:
            reports = remove_hierarchy(reports, counts, cfg.skip_hierarchy, cfg.keep_hierarchy, cfg.quiet)

        # General output file
        if len(rep_files) == 1:
            output_file = cfg.output_prefix
        else:
            file_pre = os.path.splitext(os.path.basename(rep_file))[0]
            output_file = cfg.output_prefix+file_pre

        if cfg.split_hierarchy:
            for h in reports:
                if h not in cfg.skip_hierarchy:
                    output_file_h = output_file + "." + h + ".tre"
                    r = build_report({h: reports[h]}, counts, tax, genome_sizes, output_file_h, fixed_ranks, cfg)
                    if not r:
                        print_log(" - nothing to report for hierarchy " + h + " in the " + rep_file, cfg.quiet)
                        continue
                    else:
                        print_log(" - report saved to " + output_file_h, cfg.quiet)
                        any_rep = True

        else:
            output_file = output_file + ".tre"
            r = build_report(reports, counts, tax, genome_sizes, output_file, fixed_ranks, cfg)
            if not r:
                print_log(" - nothing to report for " + rep_file, cfg.quiet)
                continue
            else:
                print_log(" - report saved to " + output_file, cfg.quiet)
                any_rep = True

    return True if any_rep else False


def parse_rep(rep_file):
    counts = {}
    reports = {}
    total_direct_matches = 0
    with open(rep_file, 'r') as rep_file:
        for line in rep_file:
            fields = line.rstrip().split("\t")
            if fields[0] == "#total_classified":
                classified_reads = int(fields[1])
            elif fields[0] == "#total_unclassified":
                unclassified_reads = int(fields[1])
            else:
                hierarchy_name, target, direct_matches, unique_reads, lca_reads, rank, name = fields
                direct_matches = int(direct_matches)
                unique_reads = int(unique_reads)
                lca_reads = int(lca_reads)

                if hierarchy_name not in reports:
                    reports[hierarchy_name] = {}
                    counts[hierarchy_name] = {"matches": 0,
                                              "reads": 0}

                # entries should be unique by hiearchy_name
                if target not in reports[hierarchy_name]:
                    reports[hierarchy_name][target] = {"direct_matches": 0,
                                                       "unique_reads": 0,
                                                       "lca_reads": 0}

                reports[hierarchy_name][target]["direct_matches"] += direct_matches
                reports[hierarchy_name][target]["unique_reads"] += unique_reads
                reports[hierarchy_name][target]["lca_reads"] += lca_reads

                counts[hierarchy_name]["matches"] += direct_matches
                # sum up classified reads for the hiearchy
                counts[hierarchy_name]["reads"] += unique_reads + lca_reads

                total_direct_matches += direct_matches

    counts["total"] = {"matches": total_direct_matches,
                       "reads": classified_reads,
                       "unclassified": unclassified_reads}

    return reports, counts


def build_report(reports, counts, tax, genome_sizes, output_file, fixed_ranks, cfg):


    # total
    if cfg.report_type == "matches":
        total = counts["total"]["matches"]
    else:
        total = counts["total"]["reads"] + counts["total"]["unclassified"]

    if len(reports) == 1:
        merged_rep = list(reports.values())[0]
    else:
        merged_rep = merge_reports(reports)

    # Re-distribute lca reads
    if cfg.report_type == "abundance":
        merged_rep = redistribute_shared_reads(merged_rep, tax)

    # Iterate over the taxonomic tree (all ranks) and sum the entries (cummulative)
    # tree_cum_counts[node] = cum_count
    tree_cum_counts = cummulative_count_tree(merged_rep, tax, cfg.report_type)
    tree_cum_perc = cummulative_perc_tree(tree_cum_counts, total)

    # Adjust percentages (tree_cum_perc) based on estimated genome sizes
    if cfg.report_type == "abundance":
        adjust_perc_genome_size(tree_cum_perc, tree_cum_counts, merged_rep, genome_sizes, total, tax, fixed_ranks)

    # filter with fixed ranks and user parameters (names, taxid)
    # filtered_cum_counts[node] = cum_count
    filtered_cum_counts = filter_report(tree_cum_counts, tree_cum_perc, tax, fixed_ranks, cfg)

    if not filtered_cum_counts:
        return False

    # sort entries based on report or user-defined
    # sorted_nodes = ["1", "1224", ...]
    sorted_nodes = sort_report(filtered_cum_counts, cfg.sort, fixed_ranks, tax, merged_rep)

    # Output file
    tre_file = open(output_file, 'w')
    output_rows = []

    # Reporting reads, first line prints unclassified entries
    if cfg.report_type in ["reads", "abundance"]:
        unclassified_line = ["unclassified",
                             "-",
                             "-",
                             "unclassified",
                             "0",
                             "0",
                             "0",
                             str(counts["total"]["unclassified"]),
                             str("%.5f" % ((counts["total"]["unclassified"]/total)*100))]
        if cfg.output_format in ["tsv", "csv"]:
            print(*unclassified_line, file=tre_file, sep="\t" if cfg.output_format == "tsv" else ",")
        else:
            output_rows.append(unclassified_line)

    orphan_nodes = 0
    # All entries
    for node in sorted_nodes:
        cum_count = filtered_cum_counts[node]
        cum_perc = tree_cum_perc[node]*100
        unique = 0
        shared = 0
        children = cum_count
        if node in merged_rep:
            unique = merged_rep[node]['unique_reads'] 
            if cfg.report_type == "matches":
                shared = merged_rep[node]['direct_matches'] - merged_rep[node]['unique_reads']
                children = children - shared
            else:
                shared = merged_rep[node]['lca_reads']
                children = children - unique - shared

        # Orphan node (not skipped in filter), reported directly to root as parent
        if tax.latest(node) == tax.undefined_node:
            orphan_nodes += 1
            lineage = [tax.root_node, node]
        else:
            if fixed_ranks:
                r = fixed_ranks.index(tax.rank(node))
                lineage = tax.lineage(node, ranks=fixed_ranks[:r+1])
            else:
                lineage = tax.lineage(node)

        out_line = [tax.rank(node),
                    node,
                    "|".join(lineage),
                    tax.name(node),
                    str(unique),
                    str(shared),
                    str(children),
                    str(cum_count),
                    str("%.5f" % cum_perc)]
        if cfg.output_format in ["tsv", "csv"]:
            print(*out_line, file=tre_file, sep="\t" if cfg.output_format == "tsv" else ",")
        else:
            output_rows.append(out_line)

    # Print formated text
    if cfg.output_format == "text":
        # Check max width for each col
        max_width = [0]*len(output_rows[0])
        for row in output_rows:
            for i, w in enumerate(max_width):
                l = len(row[i])
                if l > w:
                    max_width[i] = l
        # apply format when printing with max_width
        for row in output_rows:
            print("\t".join(['{0: <{width}}'.format(field, width=max_width[i]) for i, field in enumerate(row)]),
                  file=tre_file)

    if output_file:
        tre_file.close()

    if orphan_nodes and not cfg.no_orphan:
        print_log(" - " + str(orphan_nodes) + " orphan nodes not found in the taxonomy. " +
                  "\n   Orphan nodes are reported as 'na' with root as direct parent. " +
                  "\n   Too ommit them from the report, use --no-orphan", cfg.quiet)
    print_log(" - " + str(len(sorted_nodes)) + " taxa reported", cfg.quiet)
    return True


def merge_reports(reports, report_type):
    merged_rep = {}
    for hierarchy_name, report in reports.items():
        for target, rep in report.items():
            if target not in merged_rep:
                merged_rep[target] = {'unique_reads': 0, 'lca_reads': 0, 'direct_matches': 0}
            merged_rep[target]['unique_reads'] += rep['unique_reads']
            merged_rep[target]['lca_reads'] += rep['lca_reads']
            merged_rep[target]['direct_matches'] += rep['direct_matches']
    return merged_rep


def redistribute_shared_reads(merged_rep, tax):
    dist_rep = {}
    for target, v in merged_rep.items():
        if target not in dist_rep:
            dist_rep[target] = {'unique_reads': 0, 'lca_reads': 0, 'direct_matches': 0}

        # always keep unique reads on their original targets
        if v["unique_reads"] > 0:
            dist_rep[target]['unique_reads'] += v["unique_reads"]

        # if there are shared reads to redistribute among leaves
        if v["lca_reads"] > 0:
            leaves_unique = set()
            total_leaves = 0
            # Get leaves or return itself in case of target is already a leaf
            leaves = tax.leaves(target)
            
            # Distribute shared reads among leaves with unique reads
            redist_field = "unique_reads"
            for leaf in leaves:
                if leaf in merged_rep and merged_rep[leaf]["unique_reads"] > 0:
                    leaves_unique.add(leaf)
                    total_leaves += merged_rep[leaf]["unique_reads"]

            # If leaves of this target got no unique assignments, use the shared matches to redistribute reads
            if len(leaves_unique) == 0:
                redist_field = "direct_matches"
                for leaf in leaves:
                    if leaf in merged_rep and merged_rep[leaf]["direct_matches"] > 0:
                        leaves_unique.add(leaf)
                        total_leaves += merged_rep[leaf]["direct_matches"]

            total_redist = 0
            for leaf in leaves_unique:
                # redistribue proportionally to the number of unique assignments (or matches) for each leaf
                red = floor(v["lca_reads"]*(merged_rep[leaf][redist_field]/total_leaves))
                total_redist += red
                if leaf not in dist_rep:
                    dist_rep[leaf] = {'unique_reads': 0, 'lca_reads': 0, 'direct_matches': 0}
                dist_rep[leaf]['lca_reads'] += red

            # If there are left overs to redistribute
            left_overs = v["lca_reads"] - total_redist
            if left_overs:
                # Distribute left_over for the top leaves (by unique, matches and follow by leaf name to keep consistency in case of tie)
                for leaf in sorted(leaves_unique, key=lambda x: (-merged_rep[x]["unique_reads"], -merged_rep[x]["direct_matches"], x))[:left_overs]:
                    dist_rep[leaf]['lca_reads'] += 1
                
    return dist_rep



def adjust_perc_genome_size(tree_cum_perc, tree_cum_counts, merged_rep, genome_sizes, total, tax, fixed_ranks):
    """
    Uses genome sizes to adjust percentage of assigned reads
    It only adjusts ranks with direct assignments (unique or lca)
    It adjusts based on percentage of assigned reads to the defined rank (not 100%)
    """
    ranks_with_counts = set()
    for leaf, val in merged_rep.items():
        if val['unique_reads'] > 0 and tax.rank(leaf) in fixed_ranks:
            ranks_with_counts.add(tax.rank(leaf))
    
    total_rank_ratio = {r:0 for r in ranks_with_counts}
    total_rank_count = {r:0 for r in ranks_with_counts} 
    for node, cum_count in tree_cum_counts.items():
        if tax.rank(node) in ranks_with_counts:
            total_rank_ratio[tax.rank(node)] += cum_count/genome_sizes[node]
            total_rank_count[tax.rank(node)] += cum_count

    for node, cum_count in tree_cum_counts.items():
        if tax.rank(node) in ranks_with_counts:
            perc_adjusted = (cum_count/genome_sizes[node])/total_rank_ratio[tax.rank(node)]
            tree_cum_perc[node] = (total_rank_count[tax.rank(node)]/total) * perc_adjusted

    lowest_rank = fixed_ranks[min(map(fixed_ranks.index, ranks_with_counts))]

    # for node, adj_perc in tree_cum_perc[node].items():
    #     if tax.rank(node) == lowest_rank:
  
    #     # Node not found in tax, account as orphan direct to root
    #     if tax.latest(node) == tax.undefined_node:
    #         lin = [tax.root_node, leaf]
    #     else:
    #         lin = tax.lineage(leaf)

    #     # Cummulative sum of all nodes up to root
    #     for t in lin:
    #         if t not in filtered_cum_counts:
    #             filtered_cum_counts[t] = 0
    #         filtered_cum_counts[t] += sum_count


def cummulative_count_tree(merged_rep, tax, report_type):
    """
    Iterate over the taxonomic tree and sum the entries (cummulative) for each leaf on the report
    """
    filtered_cum_counts = {}
    for leaf in merged_rep.keys():
        if report_type == "matches":
            sum_count = merged_rep[leaf]['direct_matches'] - merged_rep[leaf]['unique_reads']
        else:
            sum_count = merged_rep[leaf]['unique_reads'] + merged_rep[leaf]['lca_reads']

        if sum_count==0:
            continue

        # Node not found in tax, account as orphan direct to root
        if tax.latest(leaf) == tax.undefined_node:
            lin = [tax.root_node, leaf]
        else:
            lin = tax.lineage(leaf)

        # Cummulative sum of all nodes up to root
        for t in lin:
            if t not in filtered_cum_counts:
                filtered_cum_counts[t] = 0
            filtered_cum_counts[t] += sum_count

    return filtered_cum_counts


def cummulative_perc_tree(tree_cum_counts, total):
    """
    Calculate percentage based on total counts on the tree
    """
    tree_cum_perc = {}
    for node, cum_count in tree_cum_counts.items():
        tree_cum_perc[node] = cum_count/total

    return tree_cum_perc


def filter_report(tree_cum_counts, tree_cum_perc, tax, fixed_ranks, cfg):
    """
    filter with fixed ranks and user parameters (names, taxid)
    """
    filtered_cum_counts = {}
    for node, cum_count in tree_cum_counts.items():
        # always keep root
        if node == tax.root_node:
            filtered_cum_counts[node] = cum_count
            continue

        # found in tax
        if tax.latest(node) != tax.undefined_node:
            # skip if not in fixed ranks
            if fixed_ranks and tax.rank(node) not in fixed_ranks:
                continue
        elif cfg.no_orphan:
            # not found in tax
            # skip if not reporting orphan nodes
            continue

        # Filter by value
        if cfg.min_count:
            if cfg.min_count > 1 and cum_count < cfg.min_count:
                continue
            elif cfg.min_count < 1 and tree_cum_perc[node] < cfg.min_count:
                continue

        if cfg.max_count:
            if cfg.max_count > 1 and cum_count > cfg.max_count:
                continue
            elif cfg.max_count < 1 and tree_cum_perc[node] > cfg.max_count:
                continue

        if cfg.taxids and not any(t in cfg.taxids for t in tax.lineage(node)):
            continue

        if cfg.names and not tax.name(node) in cfg.names:
            continue

        if cfg.names_with and not any(n in tax.name(node) for n in cfg.names_with):
            continue

        filtered_cum_counts[node] = cum_count

    return filtered_cum_counts


def sort_report(filtered_cum_counts, sort, fixed_ranks, tax, merged_rep):

    # Always keep root at the top
    # Sort report entries, return sorted keys of the dict
    if not sort:  # not user-defined, use defaults
        if not fixed_ranks:
            sorted_nodes = sorted(filtered_cum_counts,
                                  key=lambda k: tax.lineage(k))
        else:
            sorted_nodes = sorted(filtered_cum_counts,
                                  key=lambda k: (fixed_ranks.index(tax.rank(k)), -filtered_cum_counts[k]),
                                  reverse=False)
    else:  # user-defined
        if sort == "lineage":
            sorted_nodes = sorted(filtered_cum_counts, key=lambda k: tax.lineage(k))
        elif sort == "rank":
            # if sorting by rank showing all ranks, order of the ranks is not know, sort them alphabetically
            if not fixed_ranks:
                sorted_nodes = sorted(filtered_cum_counts,
                                      key=lambda k: (tax.rank(k), -filtered_cum_counts[k]),
                                      reverse=False)
            else:
                sorted_nodes = sorted(filtered_cum_counts,
                                      key=lambda k: (fixed_ranks.index(tax.rank(k)), -filtered_cum_counts[k]),
                                      reverse=False)
        elif sort == "unique":
            sorted_nodes = sorted(filtered_cum_counts,
                                  key=lambda k: (-merged_rep[k]['unique'] if k in merged_rep else 0, -filtered_cum_counts[k]),
                                  reverse=False)
        else:  # sort=="count"
            sorted_nodes = sorted(filtered_cum_counts,
                                  key=lambda k: -filtered_cum_counts[k],
                                  reverse=False)

    # Move root to the front to print always first first
    sorted_nodes.insert(0, sorted_nodes.pop(sorted_nodes.index("1")))

    return sorted_nodes


def remove_hierarchy(reports, counts, skip, keep, quiet):
    """
    Remove hiearchical levels and account all matches to root node
    """
    for hierarchy_name in list(reports.keys()):
        # Skiped report
        if hierarchy_name in skip or (keep and hierarchy_name not in keep):
            del reports[hierarchy_name]
            print_log(" - skipped " + str(counts[hierarchy_name]["reads"]) +
                      " reads with " + str(counts[hierarchy_name]["matches"]) +
                      " matches for " + hierarchy_name, quiet)

    return reports
