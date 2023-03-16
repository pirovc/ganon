from ganon.util import run, print_log, check_file
from ganon.report import report
from ganon.config import Config
from ganon.util import validate_input_files

import os
from collections import defaultdict


def reassign(cfg):

    # validate input input files
    all_files = validate_input_files(cfg.input, "all", cfg.quiet)

    print_log("Reassigning reads", cfg.quiet)

    for af in all_files:

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
            output_file = cfg.output_prefix + file_pre + ".res"

        #reverse target <-> id
        targets = {val: key for (key, val) in targets.items()}
        reassigned_reads = 0
        with open(output_file, "w") as out_file:
            for readid, matches in read_matches.items():
                if len(matches) == 1:
                    print(readid, targets[matches[0][0]], matches[0]
                          [1], sep="\t", file=out_file)
                else:
                    reassigned_reads += 1
                    max_target = matches[0]
                    max_p = 0
                    kcount = 0
                    for m, k in matches:
                        if prob[m] > max_p:
                            max_p = prob[m]
                            max_target = m
                            kcount = k
                    print(readid, targets[max_target], k, sep="\t", file=out_file)

        print_log(" - " + str(reassigned_reads) +
                  " reassigned reads [" + output_file + "]", cfg.quiet)

    return True
