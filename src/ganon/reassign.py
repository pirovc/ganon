from ganon.util import run, print_log, check_file
from ganon.report import report
from ganon.config import Config
from ganon.util import validate_input_files, rm_files

import os
from collections import defaultdict


def reassign(cfg):

    # Look for .rep
    all_files = {}
    rep_file = cfg.input + ".rep"
    rep_file_out = cfg.output_prefix + ".rep"
    rep_file_info = []
    if check_file(rep_file):
        print_log("Report file found " + rep_file, cfg.quiet)
        # look for hiearchies
        with open(rep_file) as rep:
            for line in rep:
                if line[0]!="#":
                    all_files[line.split("\t")[0]] = ""
                else:
                    rep_file_info.append([line.rstrip()])

            split_hiearchy = True
            for h in all_files.keys():
                if check_file(cfg.input + "." + h + ".all"):
                    all_files[h] = cfg.input + "." + h + ".all"
                else:
                    split_hiearchy = False
                    break

            # Single all file
            if not split_hiearchy:
                if check_file(cfg.input + ".all"):
                    all_files[""] = cfg.input + ".all"

            rm_files(rep_file_out)
    else:
        rep_file = ""
        if check_file(cfg.input + ".all"):
            all_files[""] = cfg.input + ".all"

    if not all_files:
        print_log("No .rep or .all file(s) found with prefix --input " + cfg.input, cfg.quiet)
        return False

    print_log("Reassigning reads", cfg.quiet)

    for hierarchy, af in all_files.items():

        print_log(af + (" [" + hierarchy + "]" if hierarchy else ""), cfg.quiet)

        # transoform target string into int
        targets = defaultdict(lambda: len(targets))

        read_matches = {}
        unique_matches = {}
        with open(af, "r") as all_file:
            for line in all_file:
                readid, target, kcount = line.rstrip().split("\t")
                if readid not in read_matches:
                    read_matches[readid] = []
                read_matches[readid].append((targets[target], int(kcount)))

                # Not all targets have unique matches, initialize
                if targets[target] not in unique_matches:
                    unique_matches[targets[target]] = 0

        total_umatches = 0
        for matches in read_matches.values():
            if len(matches) == 1:
                total_umatches += 1
                if matches[0][0] not in unique_matches:
                    unique_matches[matches[0][0]] = 0
                unique_matches[matches[0][0]] += 1

        # Calculate first probabilities based on unique matche
        prob = {}
        for target, unique in unique_matches.items():
            prob[target] = unique / total_umatches

        # EM loop
        for i in range(cfg.max_iter):

            total_reassigned = 0
            reassigned_matches = unique_matches.copy()

            for matches in read_matches.values():
                if len(matches) == 1:
                    continue

                # Set first match as target, in case all targets have no unique matches
                # loop will also randonly get one (last) if prob is equal
                max_target = matches[0][0]
                max_p = 0
                # print(matches)

                for m, _ in matches:
                    if prob[m] > max_p:
                        max_p = prob[m]
                        max_target = m

                reassigned_matches[max_target] += 1
                total_reassigned += 1

            diff = 0
            for target, count in reassigned_matches.items():
                new_prob = count / (total_umatches+total_reassigned)
                diff += abs(prob[target] - new_prob)
                prob[target] = new_prob

            print_log(" - Iteration " + str(i+1) + " (" + str(round(diff,6)) + ")", cfg.quiet)

            # Converged
            if diff <= cfg.threshold:
                break

        # General output file
        if len(all_files) == 1:
            output_file = cfg.output_prefix + ".res"
        else:
            file_pre = os.path.splitext(os.path.basename(af))[0]
            output_file = cfg.output_prefix + "." + hierarchy + ".res"

        #reverse target <-> id
        targets_rev = {val: key for (key, val) in targets.items()}
        reassigned_reads = 0
        with open(output_file, "w") as out_file:
            for readid, matches in read_matches.items():
                if len(matches) == 1:
                    print(readid, targets_rev[matches[0][0]], matches[0]
                          [1], sep="\t", file=out_file)
                else:
                    reassigned_reads += 1
                    max_target = matches[0][0]
                    max_p = 0
                    kcount = 0
                    for m, k in matches:
                        if prob[m] > max_p:
                            max_p = prob[m]
                            max_target = m
                            kcount = k
                    print(readid, targets_rev[max_target], k, sep="\t", file=out_file)

        print_log(" - " + str(reassigned_reads) +
                  " reassigned reads [" + output_file + "]", cfg.quiet)

        # Check if properly working
        # should I zero the lca_reads?
        if rep_file_out:
            with open(rep_file_out, "a") as rep_out:
                with open(rep_file, "r") as rep:
                    for line in rep:
                        if line[0]!="#":
                            hierarchy_name, target, direct_matches, unique_reads, lca_reads, rank, name = line.rstrip().split("\t")
                            if hierarchy_name==hierarchy and targets[target] in reassigned_matches:
                                print(hierarchy_name, target, direct_matches, reassigned_matches[targets[target]], lca_reads, rank, name, sep="\t", file=rep_out)
 
    if rep_file_out:
        with open(rep_file_out, "a") as rep_out:
            for info in rep_file_info:
                print(*info, sep="\t", file=rep_out)
        print_log(" - New report file: " + rep_file_out, cfg.quiet)


    return True
