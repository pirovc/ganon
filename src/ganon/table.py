import time
from ganon.tax import Tax
from ganon.util import *

def table(cfg):
    #validate input input files
    tre_files = validate_input_files(cfg.tre_files, cfg.input_directory, cfg.input_extension, cfg.quiet)
    
    print_log("Generating table", cfg.quiet)
    # reports[file] = {"taxa": {name: count,...}, 
    #                  "lineage": {name: ["1",...],...}, 
    #                  "label": filename, 
    #                  "total": INT, 
    #                  "unclassified_root": INT, 
    #                  "unclassified_rank": INT, 
    #                  "filtered_rank": INT}
    reports, total_taxa = parse_reports(tre_files, cfg.rank)
    print_log(" - " + str(len(reports)) + " files parsed" , cfg.quiet)
    print_log(" - " + str(total_taxa) + " total taxa selected at " + cfg.rank + " level", cfg.quiet)

    # filter reports
    filtered_total_taxa = filter_reports(reports, cfg)
    if total_taxa-filtered_total_taxa: print_log(" - " + str(total_taxa-filtered_total_taxa) + " taxa filtered out", cfg.quiet)

    # top_sample and top_all are mutually exclusive
    if cfg.top_sample:
        top_sample_total_taxa = select_top_sample(reports, cfg.top_sample)
        print_log(" - keeping " + str(top_sample_total_taxa) + "/" + str(filtered_total_taxa) + " (--top-sample "+ str(cfg.top_sample)+")", cfg.quiet)
        filtered_total_taxa = top_sample_total_taxa
    elif cfg.top_all:
        top_all_total_taxa = select_top_all(reports, cfg.top_all)
        print_log(" - keeping " + str(top_all_total_taxa) + "/" + str(filtered_total_taxa) + " (--top-all "+ str(cfg.top_all)+")", cfg.quiet)
        filtered_total_taxa = top_all_total_taxa

    if cfg.min_occurence or cfg.min_occurence_percentage:
        if cfg.min_occurence_percentage:
            cfg.min_occurence = int(len(reports)*cfg.min_occurence_percentage)
        min_occurence_total_taxa = select_occurence(reports, cfg.min_occurence)
        print_log(" - keeping " + str(min_occurence_total_taxa) + "/" + str(filtered_total_taxa) + " (--min-occurence "+ str(cfg.min_occurence)+")", cfg.quiet)
        filtered_total_taxa = min_occurence_total_taxa

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
        taxa, lineage, total, unclassified_root, unclassified_rank = parse_tre_rank(tre_file, rank)
        total_taxa.update(taxa.keys())
        reports[tre_file]["label"] = tre_file
        reports[tre_file]["taxa"] = taxa
        reports[tre_file]["lineage"] = lineage
        reports[tre_file]["total"] = total
        reports[tre_file]["unclassified_root"] = unclassified_root
        reports[tre_file]["unclassified_rank"] = unclassified_rank
        reports[tre_file]["filtered_rank"] = 0
    return reports, len(total_taxa)

def parse_tre_rank(tre_file, rank):
    taxa = {}
    lineage = {}
    unclassified_root = 0
    classified_root = 0
    classified_rank = 0
    with open(tre_file, "r") as file:
        for line in file:
            r, taxid, l, name, _, _, cum_count, _ = line.rstrip().split("\t")
            cum_count = int(cum_count)
            if r == "unclassified":
                unclassified_root=cum_count
            elif r == "root":
                classified_root=cum_count
            elif r==rank:
                classified_rank += cum_count
                taxa[name] = cum_count
                lineage[name] = l.split("|")

    unclassified_rank = classified_root-classified_rank
    total = unclassified_root + classified_root
    return taxa, lineage, total, unclassified_root, unclassified_rank


def filter_reports(reports, cfg):
    total_taxa = set()
    for file, rep in reports.items():
        for name in list(rep["taxa"].keys()):
            count = rep["taxa"][name]
            filtered = False
            if count < cfg.min_count or count/rep["total"] < cfg.min_percentage:
                filtered = True
            elif cfg.taxids and not any(t in cfg.taxids for t in rep["lineage"][name]): 
                filtered = True
            elif cfg.names and not name in cfg.names:
                filtered = True
            elif cfg.names_with and not any(n in name for n in cfg.names_with):
                filtered = True

            if filtered: 
                rep["filtered_rank"]+=count
                del rep["taxa"][name]
                del rep["lineage"][name]
            else:
                total_taxa.add(name)

    return len(total_taxa)

def select_top_sample(reports, top_sample):
    total_taxa = set()
    for file, rep in reports.items():
        for i, (name, count) in enumerate(sorted(rep["taxa"].items(), key=lambda x: x[1], reverse=True)): # sorted by count
            if i<top_sample: 
                total_taxa.add(name)
                continue
            rep["filtered_rank"]+=count
            del rep["taxa"][name]
            del rep["lineage"][name]
            
    return len(total_taxa)

def select_top_all(reports, top_all):
    total_taxa = set()
    total_counts = get_total_counts(reports)
    top_names = []
    for i,name in enumerate(sorted(total_counts, key=lambda kv: total_counts[kv]["sum_percentage"], reverse=True)):
        if i<top_all:
            top_names.append(name)

    for file, rep in reports.items():
        for name in list(rep["taxa"]):
            if name in top_names:
                total_taxa.add(name)
                continue
            rep["filtered_rank"]+=rep["taxa"][name]
            del rep["taxa"][name]
            del rep["lineage"][name]
            
    return len(total_taxa) 

def select_occurence(reports, min_occurence):
    total_taxa = set()
    total_counts = get_total_counts(reports)
    min_occ_names = []

    for name,val in total_counts.items():
        if val["occurence"]>=min_occurence:
            min_occ_names.append(name)

    for file, rep in reports.items():
        for name in list(rep["taxa"]):
            if name in min_occ_names:
                total_taxa.add(name)
                continue
            rep["filtered_rank"]+=rep["taxa"][name]
            del rep["taxa"][name]
            del rep["lineage"][name]
                
    return len(total_taxa) 

def get_total_counts(reports):
    # return sum of % of all samples
    total_counts = {}
    for d in reports.values():
        for name,count in d["taxa"].items():
            if name not in total_counts: total_counts[name] = {"sum_percentage": 0, "occurence": 0}
            total_counts[name]["sum_percentage"] += count/d["total"]
            total_counts[name]["occurence"] += 1
    return total_counts

def write_tsv(reports, cfg):
    total_counts = get_total_counts(reports)
    out_file = open(cfg.output_file, "w")

    # order by name
    sorted_names = sorted(total_counts.keys())

    if cfg.header=="taxid" or cfg.header=="lineage":
        header = [""] + [name for name in sorted_names]
        lineage = {}
        for file in reports:
            lineage.update(reports[file]["lineage"])
        if cfg.header=="taxid":
            header = [""] + [lineage[name][-1] for name in sorted_names]
        elif cfg.header=="lineage":
            header = [""] + ["|".join(lineage[name]) for name in sorted_names]
    else:
        header = [""] + [name for name in sorted_names]
    
    if cfg.add_unclassified_rank: 
        header.append("unclassified_" + cfg.rank)
    if cfg.add_unclassified: 
        header.append("unclassified")
    if cfg.add_filtered: 
        header.append("filtered")
    print(*header, sep="\t" if cfg.output_format=="tsv" else ",", file=out_file)
    cols=len(header)-1

    lines=0
    for file in sorted(reports):
        res = reports[file]
        tsv_data = [res["label"]]
        for name in sorted_names:
            if name in res["taxa"]:
                v = res["taxa"][name]
                if cfg.output_value=="percentage":
                    v = v/res["total"]
            else:
                v = 0
            tsv_data.append(v)

        if cfg.skip_zeros and (len(tsv_data) > 1 and max(tsv_data[1:]) == 0):
            print_log(" - Skipping line (" + res["label"] + ") with only zeros", cfg.quiet)
        else:
            if cfg.add_unclassified_rank:
                tsv_data.append(res["unclassified_rank"]/res["total"] if cfg.output_value=="percentage" else res["unclassified_rank"])
            if cfg.add_unclassified:
                tsv_data.append(res["unclassified_root"]/res["total"] if cfg.output_value=="percentage" else res["unclassified_root"])
            if cfg.add_filtered:
                tsv_data.append(res["filtered_rank"]/res["total"] if cfg.output_value=="percentage" else res["filtered_rank"])

            print(*tsv_data, sep="\t" if cfg.output_format=="tsv" else ",", file=out_file)
            lines+=1
            
    out_file.close()
    return lines, cols