import time
from collections import defaultdict
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
        tx = time.time()
        print_log("Generating report", cfg.quiet)
        tax = Tax([db_prefix+".tax" for db_prefix in cfg.db_prefix])
        classified_reads, unclassified_reads, reports = parse_rep(cfg.output_prefix+".rep")
        print_final_report(reports, tax, classified_reads, unclassified_reads, cfg.output_prefix+".tre", cfg.ranks, 0, 0, [])
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    return True

def report(cfg):
    classified_reads, unclassified_reads, reports = parse_rep(cfg.rep_file)
    tax = Tax([db_prefix+".tax" for db_prefix in cfg.db_prefix])
    print_final_report(reports, tax, classified_reads, unclassified_reads, cfg.output_report, cfg.ranks, cfg.min_matches, cfg.min_matches_perc, cfg.taxids)

    return True
    
def parse_rep(rep_file):
    reports = defaultdict(lambda: defaultdict(lambda: {'direct_matches': 0, 'unique_reads': 0, 'lca_reads': 0}))
    with open(rep_file, 'r') as rep_file:
        for line in rep_file:
            fields = line.rstrip().split("\t")
            if fields[0] == "#total_classified":
                seq_cla = int(fields[1])
            elif fields[0] == "#total_unclassified":
                seq_unc = int(fields[1])
            else:
                hierarchy_name, target, direct_matches, unique_reads, lca_reads, rank, name = fields
                reports[hierarchy_name][target]["direct_matches"]+=int(direct_matches)
                reports[hierarchy_name][target]["unique_reads"]+=int(unique_reads)
                reports[hierarchy_name][target]["lca_reads"]+=int(lca_reads)
    return seq_cla, seq_unc, reports

def print_final_report(reports, tax, classified_reads, unclassified_reads, final_report_file, ranks, min_matches, min_matches_perc, taxids):
    if not ranks:  
        all_ranks = False
        fixed_ranks = ['root','superkingdom','phylum','class','order','family','genus','species','species+','assembly']
    elif ranks[0]=="all":
        all_ranks = True
        fixed_ranks = []
    else:
        all_ranks = False
        fixed_ranks = ['root'] + ranks

    # sum counts of each report
    merged_rep = defaultdict(int)
    for rep in reports.values():
        for leaf in rep.keys():
            reads_assigned = rep[leaf]['unique_reads']+rep[leaf]['lca_reads']
            if reads_assigned>0:
                merged_rep[leaf]+=reads_assigned

    final_rep = defaultdict(lambda: {'count': 0, 'rank': ""})
    
    # make cummulative sum of the counts on the lineage
    for leaf in merged_rep.keys():
        count = merged_rep[leaf]
        if all_ranks: # use all nodes on the tree
            t = leaf
            r = tax.nodes[t][1]
        else: # use only nodes of the fixed ranks
            t, r = tax.get_node_rank_fixed(leaf, fixed_ranks)

        while t!="0":
            final_rep[t]['count']+=count
            final_rep[t]['rank']=r
            if all_ranks:
                t = tax.nodes[t][0]
                r = tax.nodes[t][1] if t!="0" else ""
            else:
                t, r = tax.get_node_rank_fixed(tax.nodes[t][0], fixed_ranks)

    # build lineage after all entries were defined
    lineage = {}
    for assignment in final_rep.keys():
        lineage[assignment]=[]
        if all_ranks:
            t=assignment
            while t!="0":
                lineage[assignment].insert(0,t)
                t = tax.nodes[t][0]
        else:
            t, r = tax.get_node_rank_fixed(assignment, fixed_ranks)
            max_rank_idx = fixed_ranks.index(r) # get index of current rank
            while t!="0":
                # Add empty || if fixed rank is missing
                for i in range(max_rank_idx-fixed_ranks.index(r)):
                    lineage[assignment].insert(0,"")
                    max_rank_idx-=1

                lineage[assignment].insert(0,t)
                max_rank_idx-=1
                t, r = tax.get_node_rank_fixed(tax.nodes[t][0], fixed_ranks)

        # if taxids is provided, just keep entries with them (and root)
        if taxids and assignment!="1":
            if not any(t in taxids for t in lineage[assignment]):
                del lineage[assignment]

    total_reads = classified_reads + unclassified_reads
    frfile = open(final_report_file, 'w') if final_report_file else None
    print("unclassified" +"\t"+ "-" +"\t"+ "-" +"\t"+ "-" +"\t"+ "-" +"\t"+ "-" +"\t"+ str(unclassified_reads) +"\t"+ str("%.5f" % ((unclassified_reads/total_reads)*100)), file=frfile)
    
    if all_ranks:
        sorted_assignments = sorted(lineage, key=lineage.get)
    else:
        sorted_assignments = sorted(lineage, key=lambda k: (fixed_ranks.index(final_rep[k]['rank']), -final_rep[k]['count']), reverse=False)
    
    for assignment in sorted_assignments:
        rank=final_rep[assignment]['rank'] 
        name=tax.nodes[assignment][2]
        reads_unique = rep[assignment]['unique_reads']
        reads = merged_rep[assignment] 
        matches=final_rep[assignment]['count']
        if matches < min_matches: continue
        matches_perc=(final_rep[assignment]['count']/total_reads)*100
        if matches_perc < min_matches_perc: continue
        print(rank, assignment, "|".join(lineage[assignment]), name, reads_unique, reads, matches, "%.5f" % matches_perc, file=frfile, sep="\t")
    
    if final_report_file: frfile.close()
