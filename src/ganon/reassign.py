from ganon.util import print_log, check_file, find_rep_files

import os
import pathlib
from collections import defaultdict


def reassign(cfg):
    print_log("Reassigning reads", cfg.quiet)
    print_log("", cfg.quiet)

    rep_files = list(find_rep_files(cfg.input_prefix))

    for rep_file in rep_files:
        rep_file_path = pathlib.Path(rep_file)
        rep_file_prefix = str(pathlib.Path(rep_file_path.parent, rep_file_path.stem))

        if cfg.output_prefix:
            if len(rep_files) == 1:
                out_file_prefix = cfg.output_prefix
            else:
                out_file_prefix = cfg.output_prefix + str(rep_file_path.stem)
        else:
            out_file_prefix = rep_file_prefix

        if not cfg.skip_rep:
            rep_file_out = out_file_prefix + ".rep"
        else:
            rep_file_out = ""
        rep_file_info = []

        all_files = {}
        if check_file(rep_file):
            print_log("Ganon report output found: " + str(rep_file), cfg.quiet)
            # look for hierarchies
            with rep_file.open("r") as rep:
                for line in rep:
                    if line[0] != "#":
                        all_files[line.split("\t")[0]] = ""
                    else:
                        rep_file_info.append([line.rstrip()])

            for h in all_files.keys():
                if check_file(rep_file_prefix + "." + h + ".all"):
                    # Check individual .all for multi-level hierarchy
                    all_files[h] = rep_file_prefix + "." + h + ".all"
                elif check_file(rep_file_prefix + ".all"):
                    # Check unique file for for multi-level hierarchy with --output-single
                    all_files = {}
                    all_files[""] = rep_file_prefix + ".all"
                    break
                else:
                    print_log(
                        "No matching files for given .rep ["
                        + rep_file_prefix
                        + "*.all]",
                        cfg.quiet,
                    )
                    return False
        else:
            print_log(
                "No .rep/.all file(s) found with prefix --input-prefix "
                + rep_file_prefix,
                cfg.quiet,
            )
            return False

        new_rep = []
        for hierarchy, af in all_files.items():
            print_log(af + (" [" + hierarchy + "]" if hierarchy else ""), cfg.quiet)

            # Auto-increment dict to save targets (string) into integers, less memory
            targets = defaultdict(lambda: len(targets))

            prob = {}
            read_matches = {}
            initial_weight = {}

            with open(af, "r") as all_file:
                for line in all_file:
                    readid, target, kcount = line.rstrip().split("\t")
                    if readid not in read_matches:
                        read_matches[readid] = []
                    read_matches[readid].append((targets[target], int(kcount)))

                    # Not all targets have unique matches, initialize all targets with zero weights
                    if targets[target] not in initial_weight:
                        initial_weight[targets[target]] = 0

            # Assign unique matches as initial weights
            # Calculate first probabilities based on weights (unique matches)
            total_weight = len(read_matches)
            total_initial_weight = 0
            for matches in read_matches.values():
                if len(matches) == 1:
                    total_initial_weight += 1
                    initial_weight[matches[0][0]] += 1

            if total_initial_weight == 0:
                total_initial_weight = 1

            for target, unique in initial_weight.items():
                prob[target] = unique / total_initial_weight

            # EM loop
            em_ite_cnt = 0
            while True:
                # Make copy of initial weights
                reassigned_matches = initial_weight.copy()

                for matches in read_matches.values():
                    # Skip unique matches
                    if len(matches) > 1:
                        # Get match with highest probability in this round
                        t, _ = get_top_match(matches, prob)
                        reassigned_matches[t] += 1

                diff = 0
                # Calculate new probabilities based on "simulated" re-distributed reads for the next round
                for target, count in reassigned_matches.items():
                    new_prob = count / total_weight
                    diff += abs(prob[target] - new_prob)
                    prob[target] = new_prob

                print_log(
                    " - Iteration "
                    + str(em_ite_cnt + 1)
                    + " ("
                    + str(round(diff, 6))
                    + ")",
                    cfg.quiet,
                )

                # If abs. difference among old and new probabilities already converged (no change)
                if diff <= cfg.threshold:
                    break
                if cfg.max_iter > 0 and em_ite_cnt == cfg.max_iter - 1:
                    break

                em_ite_cnt += 1

            # Skip .one file (just generate .rep)
            if not cfg.skip_one:
                # General output file
                if len(all_files) == 1:
                    one_file_out = out_file_prefix + ".one"
                else:
                    one_file_out = out_file_prefix + "." + hierarchy + ".one"

                # reverse string target <-> integer id
                targets_rev = {val: key for (key, val) in targets.items()}
                reassigned_reads = 0
                with open(one_file_out, "w") as out_file:
                    for readid, matches in read_matches.items():
                        if len(matches) == 1:
                            print(
                                readid,
                                targets_rev[matches[0][0]],
                                matches[0][1],
                                sep="\t",
                                file=out_file,
                            )
                        else:
                            reassigned_reads += 1
                            target, kcount = get_top_match(matches, prob)
                            print(
                                readid,
                                targets_rev[target],
                                kcount,
                                sep="\t",
                                file=out_file,
                            )

                print_log(
                    " - "
                    + str(reassigned_reads)
                    + " reassigned reads to "
                    + one_file_out,
                    cfg.quiet,
                )

            if rep_file_out:
                with open(rep_file, "r") as rep:
                    for line in rep:
                        if line[0] != "#":
                            fields = line.rstrip().split("\t")
                            hierarchy_name = fields[0]
                            target = fields[1]
                            direct_matches = fields[2]
                            rank = fields[5] if len(fields) >= 6 else ""
                            name = fields[6] if len(fields) >= 7 else ""
                            # Only print line of targets
                            if (
                                hierarchy == "" or hierarchy_name == hierarchy
                            ) and targets[target] in reassigned_matches:
                                # LCA matches are zero, since they will be reassigned to other nodes
                                new_rep.append(
                                    [
                                        hierarchy_name,
                                        target,
                                        direct_matches,
                                        reassigned_matches[targets[target]],
                                        0,
                                        rank,
                                        name,
                                    ]
                                )

        if rep_file_out:
            with open(rep_file_out, "w") as rep_out:
                for line in new_rep:
                    print(*line, sep="\t", file=rep_out)
                for info in rep_file_info:
                    print(*info, sep="\t", file=rep_out)
            print_log("New .rep file: " + rep_file_out, cfg.quiet)

        if cfg.remove_all:
            for af in all_files.values():
                os.remove(af)

    return True


def get_top_match(matches, prob):
    # Get top match based on prob
    # In case all targets have equal probabilities (max_p==0 at the end)
    # set first match as target (also the case for no unique matches)
    target = matches[0][0]
    kcount = matches[0][1]

    max_p = 0
    for m, k in matches:
        if prob[m] > max_p:
            max_p = prob[m]
            target = m
            kcount = k

    return target, kcount
