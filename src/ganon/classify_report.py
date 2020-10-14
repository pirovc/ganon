import time
from ganon.tax import Tax
from ganon.util import *

def classify(cfg):
    print_log("Classifying reads (ganon-classify)", cfg.quiet)
    run_ganon_classify = " ".join([cfg.path_exec['classify'],
                                   "--single-reads " +  ",".join(cfg.single_reads) if cfg.single_reads else "",
                                   "--paired-reads " +  ",".join(cfg.paired_reads) if cfg.paired_reads else "",
                                   "--ibf " + ",".join([db_prefix+".ibf" for db_prefix in cfg.db_prefix]),
                                   "--map " + ",".join([db_prefix+".map" for db_prefix in cfg.db_prefix]), 
                                   "--tax " + ",".join([db_prefix+".tax" for db_prefix in cfg.db_prefix]),
                                   "--hierarchy-labels " + ",".join(cfg.hierarchy_labels) if cfg.hierarchy_labels else "",
                                   "--max-error " + ",".join([str(me) for me in cfg.max_error]) if cfg.max_error else "",
                                   "--min-kmers " + ",".join([str(mk) for mk in cfg.min_kmers]) if cfg.min_kmers else "",
                                   "--max-error-unique " + ",".join([str(meu) for meu in cfg.max_error_unique]) if cfg.max_error_unique else "",
                                   "--strata-filter " + ",".join([str(sf) for sf in cfg.strata_filter]) if cfg.strata_filter else "",
                                   "--offset " + str(cfg.offset) if cfg.offset else "",
                                   "--output-prefix " + cfg.output_prefix if cfg.output_prefix else "",
                                   "--output-all" if cfg.output_all else "",
                                   "--output-unclassified" if cfg.output_unclassified else "",
                                   "--output-single" if cfg.output_single else "",
                                   "--threads " + str(cfg.threads) if cfg.threads else "",
                                   "--n-reads " + str(cfg.n_reads) if cfg.n_reads is not None else "",
                                   "--n-batches " + str(cfg.n_batches) if cfg.n_batches is not None else "",
                                   "--verbose" if cfg.verbose else "",
                                   "--quiet" if cfg.quiet else ""])
    stdout, stderr = run(run_ganon_classify)
    if not cfg.output_prefix: print(stdout)
    print_log(stderr, cfg.quiet)

    if cfg.output_prefix:
        tax = Tax([db_prefix+".tax" for db_prefix in cfg.db_prefix])
        total_matches, classified_reads, unclassified_reads, reports = parse_rep(cfg.output_prefix+".rep")
        print_final_report(reports, tax, total_matches, classified_reads, unclassified_reads, cfg.output_prefix+".tre", cfg.ranks, 0, 0, [], "reads")

    return True

def report(cfg):
    tx = time.time()
    print_log("Generating report", cfg.quiet)
    total_matches, classified_reads, unclassified_reads, reports = parse_rep(cfg.rep_file)
    tax = Tax([db_prefix+".tax" for db_prefix in cfg.db_prefix])
    print_final_report(reports, tax, total_matches, classified_reads, unclassified_reads, cfg.output_report, cfg.ranks, cfg.min_count, cfg.min_percentage, cfg.taxids, cfg.report_type)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    return True
    
def parse_rep(rep_file):
    reports = {}
    total_matches = 0
    with open(rep_file, 'r') as rep_file:
        for line in rep_file:
            fields = line.rstrip().split("\t")
            if fields[0] == "#total_classified":
                classified_reads = int(fields[1])
            elif fields[0] == "#total_unclassified":
                unclassified_reads = int(fields[1])
            else:
                hierarchy_name, target, direct_matches, unique_reads, lca_reads, rank, name = fields
                if hierarchy_name not in reports:
                    reports[hierarchy_name] = {}
                direct_matches = int(direct_matches)
                reports[hierarchy_name][target] = {"direct_matches":direct_matches, "unique_reads":int(unique_reads), "lca_reads":int(lca_reads)}
                total_matches+=direct_matches
    return total_matches, classified_reads, unclassified_reads, reports

def print_final_report(reports, tax, total_matches, classified_reads, unclassified_reads, final_report_file, ranks, min_count, min_percentage, taxids, report_type):
    if not reports: return False

    if not ranks:  
        all_ranks = False
        fixed_ranks = ['root','superkingdom','phylum','class','order','family','genus','species','species+','assembly']
    elif ranks[0]=="all":
        all_ranks = True
        fixed_ranks = []
    else:
        all_ranks = False
        fixed_ranks = ['root'] + ranks

    # Count targets in the report by report type, merging multiple db hierarchical levels
    # merged_report[target] = {'count': INT, 'unique': INT}
    merged_report = count_targets(reports, report_type)

    # Iterate over the taxonomic tree and sum the entries
    # final_report[node] = {'cum_count': INT, 'rank': STR}
    final_report = cummulative_count_tree(merged_report, tax, all_ranks, fixed_ranks)

    # build lineage for each entry based on chosen ranks
    # lineage[node] = ["1", "1224", ..., node]
    lineage = build_lineage(final_report, tax, all_ranks, fixed_ranks, taxids)

    frfile = open(final_report_file, 'w') if final_report_file else None
    if report_type=="reads":
        total = classified_reads + unclassified_reads
        print("unclassified" +"\t"+ "-" +"\t"+ "-" +"\t"+ "-" +"\t"+ "-" +"\t"+ "-" +"\t"+ str(unclassified_reads) +"\t"+ str("%.5f" % ((unclassified_reads/total)*100)), file=frfile)
    else:
        total = total_matches

    # Sort entries
    if all_ranks:
        sorted_nodes = sorted(lineage, key=lineage.get)
    else:
        sorted_nodes = sorted(lineage, key=lambda k: (fixed_ranks.index(final_report[k]['rank']), -final_report[k]['cum_count']), reverse=False)
    
    for node in sorted_nodes:
        rank=final_report[node]['rank'] 
        name=tax.nodes[node][2]
        unique = merged_report[node]['unique'] if node in merged_report else 0
        all_count = merged_report[node]['unique'] + merged_report[node]['count'] if node in merged_report else 0
        cum_count=final_report[node]['cum_count']
        if cum_count < min_count: continue
        cum_count_perc=(final_report[node]['cum_count']/total)*100
        if cum_count_perc < min_percentage: continue
        print(rank, node, "|".join(lineage[node]), name, unique, all_count, cum_count, "%.5f" % cum_count_perc, file=frfile, sep="\t")
    
    if final_report_file: frfile.close()

def count_targets(reports, report_type):
    merged_report = {}
    if report_type=="reads":
        for hierarchy_name,report in reports.items():
            for target,rep in report.items():
                # if there were reads assigned to the target (not only shared matches)
                if rep['unique_reads'] + rep['lca_reads']:
                    if target not in merged_report:
                        merged_report[target] = {'unique':0, 'count': 0}
                    merged_report[target]['unique'] += rep['unique_reads']
                    merged_report[target]['count'] += rep['lca_reads']
    else:
        for hierarchy_name,report in reports.items():
            for target,rep in report.items():
                # If there were any matches to the target
                if rep['direct_matches']:
                    if target not in merged_report:
                        merged_report[target] = {'unique':0, 'count': 0}
                    merged_report[target]['unique'] += rep['unique_reads']
                    # count already has unique, so it needs to be subtracted to fit the same model as the report type reads
                    merged_report[target]['count'] += rep['direct_matches']-rep['unique_reads']

    return merged_report

def cummulative_count_tree(merged_report, tax, all_ranks, fixed_ranks):
    final_report = {}
    for leaf in merged_report.keys():
        sum_count = merged_report[leaf]['unique'] + merged_report[leaf]['count']
        if all_ranks: # Use all nodes of the tree
            t = leaf
            r = tax.nodes[t][1]
            while t!="0":
                if t not in final_report: final_report[t] = {'cum_count': 0, 'rank': ""}
                final_report[t]['cum_count']+=sum_count
                final_report[t]['rank']=r
                t = tax.nodes[t][0]
                r = tax.nodes[t][1] if t!="0" else ""
        else: # Use selected nodes
            # get closest node of fixed ranks
            t, r = tax.get_node_rank_fixed(leaf, fixed_ranks)
            while t!="0":
                if t not in final_report: final_report[t] = {'cum_count': 0, 'rank': ""}
                final_report[t]['cum_count']+=sum_count
                final_report[t]['rank']=r
                t, r = tax.get_node_rank_fixed(tax.nodes[t][0], fixed_ranks)

    return final_report

def build_lineage(final_report, tax, all_ranks, fixed_ranks, taxids):
    lineage = {}
    for node in final_report.keys():
        lineage[node]=[]
        if all_ranks:
            t=node
            while t!="0":
                lineage[node].insert(0,t)
                t = tax.nodes[t][0]
        else:
            t, r = tax.get_node_rank_fixed(node, fixed_ranks)
            max_rank_idx = fixed_ranks.index(r) # get index of current rank
            while t!="0":
                # Add empty || if fixed rank is missing
                for i in range(max_rank_idx-fixed_ranks.index(r)):
                    lineage[node].insert(0,"")
                    max_rank_idx-=1
                lineage[node].insert(0,t)
                max_rank_idx-=1
                t, r = tax.get_node_rank_fixed(tax.nodes[t][0], fixed_ranks)

        # if taxids is provided, just keep entries with them (and root)
        if taxids and node!="1":
            if not any(t in taxids for t in lineage[node]):
                del lineage[node]

    return lineage