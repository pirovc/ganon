import time
from ganon.tax import Tax
from ganon.util import *
from collections import OrderedDict

def table(cfg):

    tx = time.time()
    print_log("Generating table", cfg.quiet)
    print_log(" - Parsing " + str(len(cfg.tre_files)) + " files" , cfg.quiet)
    results = OrderedDict()
    for file in sorted(cfg.tre_files): # sorted by file name
        results[file] = dict()
        results[file]["label"] = file
        results[file]["data"], results[file]["total"], results[file]["unclassified_root"], results[file]["unclassified_rank"], results[file]["filtered_rank"] = parse_tre(file, cfg)
    
    total_counts = get_total_counts(results)
    print_log(" - " + str(len(total_counts)) + " taxa found at " + cfg.rank + " level", cfg.quiet)

    # filter results based on top hits to all samples
    if cfg.top_all:
        top_names = []
        for i,name in enumerate(sorted(total_counts, key=lambda kv: total_counts[kv]["sum_percentage"], reverse=True)):
            if i<cfg.top_all:
                top_names.append(name)

        for d in results.values():
            tre, classified_read_count, filtered_read_count = filter_names(d["data"], top_names)
            d["data"] = tre
            d["filtered_rank"]+=filtered_read_count

        print_log(" - keeping " + str(len(top_names)) + "/" + str(len(total_counts)) + " top taxa among all files", cfg.quiet)
    
        # reset total counts
        total_counts = get_total_counts(results)

    if cfg.min_occurence:
        min_occ_names = []
        for name,val in total_counts.items():
            if val["occurence"]>=cfg.min_occurence:
                min_occ_names.append(name)

        for d in results.values():
            tre, classified_read_count, filtered_read_count = filter_names(d["data"], min_occ_names)
            d["data"] = tre
            d["filtered_rank"]+=filtered_read_count

        print_log(" - keeping " + str(len(min_occ_names)) + "/" + str(len(total_counts)) + " taxa occuring in " + str(cfg.min_occurence) + " or more files", cfg.quiet)
    
        # reset total counts
        total_counts = get_total_counts(results)


    lines = write_tsv(results, total_counts.keys(), cfg)
    print_log(" - " + str(len(total_counts.keys())) + "x" + str(lines) + " table saved to " + cfg.output_file, cfg.quiet)
    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    return True

def parse_tre(tre_file, cfg):
    tre = dict()

    unclassified_root = 0
    classified_root = 0
    unclassified_rank = 0
    classified_rank = 0
    filtered_rank = 0

    with open(tre_file, "r") as file:
        for line in file:
            rank, taxid, _, name, _, _, read_count, percentage = line.rstrip().split("\t")
            read_count = int(read_count)
            percentage = float(percentage)/100
            if rank == "unclassified":
                unclassified_root=read_count
            elif rank == "root":
                classified_root=read_count
            elif rank==cfg.rank:
                if read_count<cfg.min_count or percentage<cfg.min_percentage:
                    filtered_rank += read_count
                    continue
                elif cfg.taxids and taxid not in cfg.taxids:
                    filtered_rank += read_count
                    continue
                classified_rank += read_count
                tre[name] = read_count
                    

    # filter only top hits of each sample
    if cfg.top_sample:
        top = dict()
        classified_rank = 0 #restart to count
        for idx, (name, read_count) in enumerate(sorted(tre.items(), key=lambda kv: kv[1], reverse=True)): # sorted by read_count
            if idx<cfg.top_sample:
                top[name] = read_count
                classified_rank += read_count
            else:
                filtered_rank += read_count # add to already filtered
        tre = top

    # keep only selected names
    if cfg.names:
        tre, classified_read_count, filtered_read_count = filter_names(tre, cfg.names)
        classified_rank = classified_read_count # restart to count
        filtered_rank+=filtered_read_count # just add elements removed from classified
    elif cfg.names_with:
        tre, classified_read_count, filtered_read_count = filter_names(tre, cfg.names_with, True)
        classified_rank = classified_read_count # restart to count
        filtered_rank+=filtered_read_count # just add elements removed from classified

    # read without a match on the chosen rank
    unclassified_rank = classified_root-classified_rank-filtered_rank

    total = unclassified_root + classified_root

    return tre, total, unclassified_root, unclassified_rank, filtered_rank

def filter_names(tre, names, name_with=False):
    selected_names = dict()
    classified_read_count = 0
    filtered_read_count = 0
    for name, read_count in tre.items():
        if not name_with:
            found=name in names
        else:
            found=any(n in name for n in names)

        if found:
            selected_names[name] = read_count
            classified_read_count += read_count
        else:
            filtered_read_count += read_count # add to already filtered
    
    return selected_names, classified_read_count, filtered_read_count

def get_total_counts(results):
    # return sum of % of all samples
    total_counts = {}
    for d in results.values():
        for name,count in d["data"].items():
            if name not in total_counts: total_counts[name] = {"sum_percentage": 0, "occurence": 0}
            total_counts[name]["sum_percentage"] += count/d["total"]
            total_counts[name]["occurence"] += 1
    return total_counts

def write_tsv(results, names, cfg):
    sorted_names = sorted(names)
    lines=0
    tsv_file = open(cfg.output_file, "w")
    header = [""] + [name for name in sorted_names]

    if cfg.add_unclassified_rank: header.append("unclassified_" + cfg.rank)
    if cfg.add_unclassified: header.append("unclassified")
    if cfg.add_filtered: header.append("filtered")

    print(*header, sep="\t", file=tsv_file)
    for res in results.values():
        tsv_data = [res["label"]]
        for name in sorted_names:
            if name in res["data"]:
                v = res["data"][name]
                if cfg.output_value=="percentage":
                    v = v/res["total"]
            else:
                v = 0
            tsv_data.append(v)

        if cfg.skip_zeros and (len(tsv_data) > 1 and max(tsv_data[1:]) == 0):
            print_log(" - Skipping line (" + res["label"] + ") with only zeros")
        else:
            if cfg.add_unclassified_rank:
                tsv_data.append(res["unclassified_rank"]/res["total"] if cfg.output_value=="percentage" else res["unclassified_rank"])
            if cfg.add_unclassified:
                tsv_data.append(res["unclassified_root"]/res["total"] if cfg.output_value=="percentage" else res["unclassified_root"])
            if cfg.add_filtered:
                tsv_data.append(res["filtered_rank"]/res["total"] if cfg.output_value=="percentage" else res["filtered_rank"])

            print(*tsv_data, sep="\t", file=tsv_file)
            lines+=1
            
    tsv_file.close()
    return lines