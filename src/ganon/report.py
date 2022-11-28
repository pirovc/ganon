import time
import os
from math import floor, ceil

from copy import deepcopy
from ganon.util import validate_input_files
from ganon.util import print_log
from ganon.tax_util import get_genome_size, parse_genome_size_tax

from multitax import CustomTx, NcbiTx, GtdbTx


def report(cfg):
    # validate input input files
    rep_files = validate_input_files(cfg.input, cfg.input_extension, cfg.quiet)

    # Parse taxonomy or download new
    tax_args = {"undefined_node": "",
                "undefined_rank": "na",
                "undefined_name": "na",
                "root_rank": "root",
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

        if cfg.report_type in ["abundance", "corr"]:
            try:
                genome_sizes = parse_genome_size_tax(dbp)
            except ValueError:
                print_log(
                    "Failed to get genome sizes from .tax files, run ganon report without -d/--db-prefix")
                return False
    else:
        tx = time.time()
        if cfg.taxonomy_files:
            print_log("Parsing " + cfg.taxonomy + " taxonomy", cfg.quiet)
        else:
            print_log("Downloading and parsing " +
                      cfg.taxonomy + " taxonomy", cfg.quiet)

        if cfg.taxonomy == "ncbi":
            tax = NcbiTx(files=cfg.taxonomy_files, **tax_args)
        elif cfg.taxonomy == "gtdb":
            tax = GtdbTx(files=cfg.taxonomy_files, **tax_args)

        print_log(" - done in " + str("%.2f" %
                                      (time.time() - tx)) + "s.\n", cfg.quiet)

        # In case no tax was provided, generate genome sizes (for the full tree)
        if cfg.report_type in ["abundance", "corr"]:
            genome_sizes = get_genome_size(cfg, tax.leaves(), tax, "./")

    default_ranks = [tax.root_name] + cfg.choices_default_ranks

    # define fixed_ranks or leave it empty for all
    if cfg.ranks and cfg.ranks[0] == "all":
        fixed_ranks = []
    else:
        if not cfg.ranks or cfg.ranks == [""]:
            fixed_ranks = default_ranks
        else:
            fixed_ranks = [tax.root_name] + cfg.ranks

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
            reports = remove_hierarchy(
                reports, counts, cfg.skip_hierarchy, cfg.keep_hierarchy, cfg.quiet)

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
                    r = build_report(
                        {h: reports[h]}, counts, tax, genome_sizes, output_file_h, fixed_ranks, default_ranks, cfg, rep_file)
                    if not r:
                        print_log(" - nothing to report for hierarchy " +
                                  h + " in the " + rep_file, cfg.quiet)
                        continue
                    else:
                        print_log(" - report saved to " +
                                  output_file_h, cfg.quiet)
                        any_rep = True

        else:
            output_file = output_file + ".tre"
            r = build_report(reports, counts, tax, genome_sizes,
                             output_file, fixed_ranks, default_ranks, cfg, rep_file)
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


def build_report(reports, counts, full_tax, genome_sizes, output_file, fixed_ranks, default_ranks, cfg, rep_file):

    # total
    if cfg.report_type == "matches":
        total = counts["total"]["matches"]
    else:
        total = counts["total"]["reads"] + counts["total"]["unclassified"]

    # Merge .rep files into a dict {'unique_reads': 0, 'lca_reads': 0, 'direct_matches': 0}
    if len(reports) == 1:
        merged_rep = list(reports.values())[0]
    else:
        merged_rep = merge_reports(reports)

    # Copy full tax and only keep used subset (for consistency with downloaded taxonomy)
    tax = deepcopy(full_tax)
    tax.filter(list(merged_rep.keys()))

    orphan_nodes = set()
    # Add orphan nodes to taxonomy
    for node in merged_rep.keys():
        if tax.latest(node) == tax.undefined_node:
            tax._nodes[node] = tax.root_node
            tax._ranks[node] = tax.undefined_rank
            orphan_nodes.add(node)
    tax.check_consistency()
    # Pre-build lineages for performance
    tax.build_lineages()

    # Re-distribute lca reads to leaf nodes based on unique matches or shared
    if cfg.report_type in ["abundance", "dist"]:
        redistribute_shared_reads(merged_rep, tax)

    # Count data from merged_rep into a final count per target, depending on report type
    target_counts = count_targets(merged_rep, cfg.report_type)

    # Cummulatively sum leaf counts to its lineage parents
    tree_cum_counts = cummulative_sum_tree(target_counts, tax)

    if cfg.report_type in ["abundance", "corr"]:
        # Correct counts based on estimated genome sizes for default ranks
        # returns tree_cum_counts with corrected FRACTION of counts and use it to calculate the abundances
        # still reports original counts for user
        corr_tree_cum_counts = correct_genome_size(
            target_counts, genome_sizes, tax, default_ranks)
        tree_cum_perc = cummulative_perc_tree(corr_tree_cum_counts, total)
    else:
        # Simple percentage calculation
        tree_cum_perc = cummulative_perc_tree(tree_cum_counts, total)

    # filter with fixed ranks and user parameters (names, taxid)
    # filtered_cum_counts[node] = cum_count
    filtered_cum_counts = filter_report(
        tree_cum_counts, tree_cum_perc, tax, fixed_ranks, default_ranks, orphan_nodes, cfg)

    if not filtered_cum_counts:
        return False

    # sort entries based on report or user-defined
    # sorted_nodes = ["1", "1224", ...]
    sorted_nodes = sort_report(
        filtered_cum_counts, tree_cum_perc, cfg.sort, fixed_ranks, tax, merged_rep)

    # Output file
    tre_file = open(output_file, 'w')

    if cfg.output_format == "bioboxes":
        print("@Version:0.10.0", file=tre_file)
        print("@SampleID:" + rep_file + " " + ",".join(reports.keys()), file=tre_file)
        print("@Ranks:" + "|".join(fixed_ranks[1:]), file=tre_file)
        print("@Taxonomy:" + ",".join(tax.sources), file=tre_file)
        print("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE", file=tre_file)
        for node in sorted_nodes:
            cum_perc = tree_cum_perc[node]*100
            # Do not report root
            if node == tax.root_node:
                continue

            cum_perc = tree_cum_perc[node]*100

            if fixed_ranks:
                r = fixed_ranks.index(tax.rank(node))
                lineage = tax.lineage(node, ranks=fixed_ranks[:r+1])
                name_lineage = tax.name_lineage(node, ranks=fixed_ranks[:r+1])
            else:
                lineage = tax.lineage(node)
                name_lineage = tax.name_lineage(node)

            out_line = [node,
                        tax.rank(node),
                        "|".join(lineage[1:]), #ignore root
                        "|".join(name_lineage[1:]), #ignore root
                        str("%.5f" % cum_perc)]
            print(*out_line, file=tre_file, sep="\t")

    else:
        output_rows = []
        # Reporting reads, first line prints unclassified entries
        if cfg.report_type != "matches":
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
                print(*unclassified_line, file=tre_file,
                      sep="\t" if cfg.output_format == "tsv" else ",")
            else:
                output_rows.append(unclassified_line)

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
                    shared = merged_rep[node]['direct_matches'] - \
                        merged_rep[node]['unique_reads']
                else:
                    shared = merged_rep[node]['lca_reads']

            children = children - unique - shared

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
                print(*out_line, file=tre_file,
                      sep="\t" if cfg.output_format == "tsv" else ",")
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
        print_log(" - " + str(len(orphan_nodes)) + " not found in the taxonomy (orphan nodes). " +
                  "\n   Orphan nodes are reported with 'na' rank with root as a direct parent node. " +
                  "\n   Too show them, use 'na' in --ranks or set --ranks all"
                  "\n   Too ommit them, use --no-orphan", cfg.quiet)
    print_log(" - " + str(len(sorted_nodes)) + " entries reported", cfg.quiet)
    return True


def merge_reports(reports):
    merged_rep = {}
    for hierarchy_name, report in reports.items():
        for target, rep in report.items():
            if target not in merged_rep:
                merged_rep[target] = {'unique_reads': 0,
                                      'lca_reads': 0, 'direct_matches': 0}
            merged_rep[target]['unique_reads'] += rep['unique_reads']
            merged_rep[target]['lca_reads'] += rep['lca_reads']
            merged_rep[target]['direct_matches'] += rep['direct_matches']
    return merged_rep


def count_targets(merged_rep, report_type):
    res = {}
    for target, v in merged_rep.items():
        count = 0
        if report_type == "matches":
            count = v['direct_matches']
        else:
            count = v['unique_reads'] + v['lca_reads']
        if count == 0:
            continue
        res[target] = count

    return res


def redistribute_shared_reads(merged_rep, tax):
    """
    change merged_rep with redistributed reads
    only move counts of lca_reads
    """
    for target in merged_rep.keys():

        # if there are shared reads to redistribute among leaves
        if merged_rep[target]["lca_reads"] > 0:

            # Get leaves or return itself in case of target is already a leaf
            leaves = tax.leaves(target)

            # If returned itself, no redistribution necessary / not found in tax, skip
            if not leaves or leaves == [target]:
                continue

            # Distribute shared reads among leaves with unique reads
            redist_field = "unique_reads"
            total_leaves = 0
            leaves_unique = set()
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

            # If no leaves could be found, skip
            if len(leaves_unique) == 0:
                continue

            total_redist = 0
            for leaf in leaves_unique:
                # redistribue proportionally to the number of unique assignments (or matches) for each leaf
                red = floor(merged_rep[target]["lca_reads"] *
                            (merged_rep[leaf][redist_field]/total_leaves))
                total_redist += red
                if leaf not in merged_rep:
                    merged_rep[leaf] = {'unique_reads': 0,
                                        'lca_reads': 0, 'direct_matches': 0}
                merged_rep[leaf]['lca_reads'] += red

            # If there are left overs to redistribute
            left_overs = merged_rep[target]["lca_reads"] - total_redist
            if left_overs:
                # Distribute left_over for the top leaves
                # (by unique, matches and follow by leaf name to keep consistency in case of tie)
                for leaf in sorted(leaves_unique, key=lambda x: (-merged_rep[x]["unique_reads"], -merged_rep[x]["direct_matches"], x))[:left_overs]:
                    merged_rep[leaf]['lca_reads'] += 1

            # Remove redistributed reads from target
            merged_rep[target]["lca_reads"] = 0


def correct_genome_size(target_counts, genome_sizes, tax, default_ranks):
    """
    Correct counts with genome sizes based only on default ranks.
    Leaf or nodes in between are counted for its closest default parent
    Other ranks are ignored on the correction but re-inserted later
    Correct counts (allowing fractions) independently for each tax. level and regenerate cum tree at the end
    """
    ranked_counts = {}
    lost_targets = {}
    total_rank_ratio = {r: 0 for r in default_ranks}
    total_rank_count = {r: 0 for r in default_ranks}
    for target, count in target_counts.items():

        # Define closest parent based on default ranks and save count
        closest_parent = tax.closest_parent(target, ranks=default_ranks)
        if closest_parent not in ranked_counts:
            ranked_counts[closest_parent] = 0
        ranked_counts[closest_parent] += count

        # Keep track of targets merged into default ranks
        if closest_parent != target:
            lost_targets[target] = closest_parent

        # Sum total counts for each default rank
        gs = genome_sizes[closest_parent] if closest_parent in genome_sizes else genome_sizes[tax.root_node]
        closest_rank = tax.rank(closest_parent)
        total_rank_ratio[closest_rank] += count/gs
        total_rank_count[closest_rank] += count

    # Correct counts by the genome sizes (only default ranks)
    corr_counts = {t: 0 for t in ranked_counts.keys()}
    for node in ranked_counts.keys():
        rank_node = tax.rank(node)
        gs = genome_sizes[node] if node in genome_sizes else genome_sizes[tax.root_node]
        corr_counts[node] = total_rank_count[rank_node] * \
            ((ranked_counts[node]/gs)/total_rank_ratio[rank_node])

    # Corrected counts should match original
    assert sum(target_counts.values()) == round(
        sum(corr_counts.values())), "invalid number of counts after correction"

    # Generate cumulative tree based on corrected counts
    corr_tree = cummulative_sum_tree(corr_counts, tax)

    # Re-insert lost targets in the tree
    # They can be either leaf or ranks in between default ranks (e.g. strain or suborder)
    # Step necessary to report any rank in abundance node (default rank nodes are unchanged)
    for target, closest_parent in lost_targets.items():
        # all ranks in between default ranks will receive proportional counts from parent
        for t in tax.lineage(target, root_node=closest_parent)[1:]:
            # In case of a leaf node, re-add to the tree
            if t not in corr_tree:
                corr_tree[t] = 0
            # Sum proportional amount of reads of default parent node to all nodes in between
            corr_tree[t] += target_counts[target] * \
                (corr_counts[closest_parent]/ranked_counts[closest_parent])

    return corr_tree


def cummulative_sum_tree(target_count, tax):
    """
    Iterate over the taxonomic tree and sum the values of a dict cummulatively {node: value}
    """
    cum_counts = {}
    for target, count in target_count.items():
        # Cummulative sum of all nodes up to root
        for t in tax.lineage(target):
            if t not in cum_counts:
                cum_counts[t] = 0
            cum_counts[t] += count
    return cum_counts


def cummulative_perc_tree(tree_cum_counts, total):
    """
    Calculate percentage based on total counts on the tree
    """
    tree_cum_perc = {}
    for node, cum_count in tree_cum_counts.items():
        tree_cum_perc[node] = cum_count/total

    return tree_cum_perc


def filter_report(tree_cum_counts, tree_cum_perc, tax, fixed_ranks, default_ranks, orphan_nodes, cfg):
    """
    filter with fixed ranks and user parameters (names, taxid)
    """
    filtered_cum_counts = {}

    filter_counts_msg = {}
    filter_counts_msg["orphan"] = {"count": 0, "msg": "orphan entries removed"}
    filter_counts_msg["ranks"] = {"count": 0, "msg": "entries removed not in --ranks [" + ",".join(fixed_ranks[1:]) if fixed_ranks else "" + "]"}
    filter_counts_msg["percentile"] = {"count": 0, "msg": "entries removed with --top-percentile " + str(cfg.top_percentile)}
    filter_counts_msg["min_count"] = {"count": 0, "msg": "entries removed with --min-count " + str(cfg.min_count)}
    filter_counts_msg["max_count"] = {"count": 0, "msg": "entries removed with --max-count " + str(cfg.min_count)}
    filter_counts_msg["taxids"] = {"count": 0, "msg": "entries removed not in --taxids [" + ",".join(cfg.taxids) + "]"}
    filter_counts_msg["names"] = {"count": 0, "msg": "entries removed not in --names [" + ",".join(cfg.names) + "]"}
    filter_counts_msg["names_with"] = {"count": 0, "msg": "entries removed not in --names-with [" + ",".join(cfg.names_with) + "]"}


    rank_cutoff_percentile = {}
    # Detect cut-off for top percentile
    if cfg.top_percentile:
        rank_perc = {r: [] for r in default_ranks}
        # Sorted percentages/abundance by rank (only default)
        for node, perc in sorted(tree_cum_perc.items(), key=lambda x: x[1], reverse=True):
            rank = tax.rank(node)
            if rank in default_ranks:
                rank_perc[rank].append(perc)
        # Define threhsold for percentile

        for rank, perc_list in rank_perc.items():
            top = ceil(cfg.top_percentile * len(perc_list))
            if top < len(perc_list):
                rank_cutoff_percentile[rank] = perc_list[top]

    for node, cum_count in tree_cum_counts.items():
        rank = tax.rank(node)
        # always keep root
        if node == tax.root_node:
            filtered_cum_counts[node] = cum_count
            continue

        # Skip orphan nodes
        if node in orphan_nodes and cfg.no_orphan:
            filter_counts_msg["orphan"]["count"]+=1
            continue

        # skip if not in fixed ranks
        if fixed_ranks and rank not in fixed_ranks:
            filter_counts_msg["ranks"]["count"]+=1
            continue

        # Filter by top percentile (if rank_cutoff_percentile is populated)
        if rank in rank_cutoff_percentile and tree_cum_perc[node] <= rank_cutoff_percentile[rank]:
            filter_counts_msg["percentile"]["count"]+=1
            continue

        # Filter by value
        if cfg.min_count:
            if cfg.min_count > 1 and cum_count < cfg.min_count:
                filter_counts_msg["min_count"]["count"]+=1
                continue
            elif cfg.min_count < 1 and tree_cum_perc[node] < cfg.min_count:
                filter_counts_msg["min_count"]["count"]+=1
                continue

        if cfg.max_count:
            if cfg.max_count > 1 and cum_count > cfg.max_count:
                filter_counts_msg["max_count"]["count"]+=1
                continue
            elif cfg.max_count < 1 and tree_cum_perc[node] > cfg.max_count:
                filter_counts_msg["max_count"]["count"]+=1
                continue

        if cfg.taxids and not any(t in cfg.taxids for t in tax.lineage(node)):
            filter_counts_msg["taxids"]["count"]+=1
            continue

        if cfg.names and not tax.name(node) in cfg.names:
            filter_counts_msg["names"]["count"]+=1
            continue

        if cfg.names_with and not any(n in tax.name(node) for n in cfg.names_with):
            filter_counts_msg["names_with"]["count"]+=1
            continue

        filtered_cum_counts[node] = cum_count


    for f in filter_counts_msg.keys():
        if filter_counts_msg[f]["count"] > 0:
            print_log(" - " + str(filter_counts_msg[f]["count"]) + " " + filter_counts_msg[f]["msg"])

    return filtered_cum_counts


def sort_report(filtered_cum_counts, tree_cum_perc, sort, fixed_ranks, tax, merged_rep):
    # Always keep root at the top
    # Sort report entries, return sorted keys of the dict
    if not sort:  # not user-defined, use defaults
        if not fixed_ranks:
            sorted_nodes = sorted(filtered_cum_counts,
                                  key=lambda k: tax.lineage(k))
        else:
            # Add undefined node to fixed ranks in case of reporting
            sort_fixed_ranks = fixed_ranks + [tax.undefined_rank]
            sorted_nodes = sorted(filtered_cum_counts,
                                  key=lambda k: (sort_fixed_ranks.index(
                                      tax.rank(k)), -tree_cum_perc[k]),
                                  reverse=False)
    else:  # user-defined
        if sort == "lineage":
            sorted_nodes = sorted(filtered_cum_counts,
                                  key=lambda k: tax.lineage(k))
        elif sort == "rank":
            # if sorting by rank showing all ranks
            # order of the ranks is not known, sort them alphabetically
            if not fixed_ranks:
                sorted_nodes = sorted(filtered_cum_counts,
                                      key=lambda k: (
                                          tax.rank(k), -tree_cum_perc[k]),
                                      reverse=False)
            else:
                # Add undefined node to fixed ranks in case of reporting
                sort_fixed_ranks = fixed_ranks + [tax.undefined_rank]
                sorted_nodes = sorted(filtered_cum_counts,
                                      key=lambda k: (sort_fixed_ranks.index(
                                          tax.rank(k)), -tree_cum_perc[k]),
                                      reverse=False)
        elif sort == "unique":
            sorted_nodes = sorted(filtered_cum_counts,
                                  key=lambda k: (-merged_rep[k]['unique']
                                                 if k in merged_rep else 0, -tree_cum_perc[k]),
                                  reverse=False)
        elif sort == "count":
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
