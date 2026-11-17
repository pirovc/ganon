from ganon.util import print_log, check_file

import os
import pathlib
from collections import defaultdict


def reassign(cfg):
    print_log("Reassigning reads", cfg.quiet)
    rep_file = cfg.input_prefix + ".rep"
    rep_file_path = pathlib.Path(rep_file)
    rep_file_prefix = str(pathlib.Path(rep_file_path.parent, rep_file_path.stem))

    if cfg.output_prefix:
        out_file_prefix = cfg.output_prefix
    else:
        out_file_prefix = rep_file_prefix

    if not cfg.skip_rep:
        rep_file_out = out_file_prefix + ".rep"
    else:
        rep_file_out = ""
    rep_file_info = []

    all_files = {}
    if check_file(rep_file):
        print_log("Input: " + str(rep_file), cfg.quiet)
        # look for hierarchies
        with rep_file_path.open("r") as rep:
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
                    "No matching files for given .rep [" + rep_file_prefix + "*.all]",
                    cfg.quiet,
                )
                return False
    else:
        print_log(
            "No .rep/.all file(s) found with prefix --input-prefix " + rep_file_prefix,
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
        unique_counts = {}

        with open(af, "r") as all_file:
            for line in all_file:
                readid, target, kcount = line.rstrip().split("\t")
                if readid not in read_matches:
                    read_matches[readid] = []
                read_matches[readid].append((targets[target], int(kcount)))

        # initialize all targets with zero weights
        unique_counts = {tid: 0 for tid in targets.values()}

        # Initialization
        for matches in read_matches.values():
            if len(matches) == 1:
                unique_counts[matches[0][0]] += 1

        # Initialize prob. target considering unique reads + epsilon
        # epsilon as a pseudocount to not have 0 initial weights
        EPS = 1e-9
        n_targets = len(unique_counts)

        n_reads = len(read_matches)

        n_unique = sum(unique_counts.values())
        for target, unique in unique_counts.items():
            prob[target] = (unique + EPS) / (n_unique + EPS * n_targets)

        # EM loop
        em_ite_cnt = 0
        while True:
            expected_counts = {t: float(u) for t, u in unique_counts.items()}
            for matches in read_matches.values():
                # Unique reads are already fully counted in unique_counts
                if len(matches) == 1:
                    continue

                # E-step: unnormalised responsibilities
                weights = []
                denom = 0.0
                for target, kcount in matches:
                    # Use kmer count as a relative weigth for the match
                    w = prob[target] * kcount
                    weights.append((target, w))
                    denom += w

                # should not happend with EPS
                if denom == 0.0:
                    share = 1.0 / len(matches)
                    for target, _ in matches:
                        expected_counts[target] += share
                else:
                    # Add weighted observation to the counts
                    for target, w in weights:
                        expected_counts[target] += w / denom

            # M-step: p_t = expected_count_t / n_reads
            # Calculate new probability based on expected_counts
            diff = 0.0
            for target, count in expected_counts.items():
                new_prob = count / n_reads
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

        # Final hard E-step: one integer count per read
        reassigned_matches = unique_counts.copy()
        for matches in read_matches.values():
            if len(matches) > 1:
                target, _ = get_top_match(matches, prob)
                reassigned_matches[target] += 1

        # Sum = #unique reads + #multi reads = R
        assert sum(reassigned_matches.values()) == n_reads

        # Skip .one file (just generate .rep)
        if not cfg.skip_one:
            # General output file
            if len(all_files) == 1:
                one_file_out = out_file_prefix + ".one"
            else:
                one_file_out = out_file_prefix + "." + hierarchy + ".one"

            # reverse string target <-> integer id
            targets_rev = {val: key for (key, val) in targets.items()}
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
                        target, kcount = get_top_match(matches, prob)
                        print(
                            readid,
                            targets_rev[target],
                            kcount,
                            sep="\t",
                            file=out_file,
                        )

            print_log(
                "Output (.one): " + one_file_out,
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
                        unique = int(fields[3])
                        rank = fields[5] if len(fields) >= 6 else ""
                        name = fields[6] if len(fields) >= 7 else ""
                        # Only print line of targets
                        if (hierarchy == "" or hierarchy_name == hierarchy) and targets[
                            target
                        ] in reassigned_matches:
                            # LCA matches are zero, since they will be reassigned to other nodes
                            new_rep.append(
                                [
                                    hierarchy_name,
                                    target,
                                    direct_matches,
                                    unique,
                                    reassigned_matches[targets[target]] - unique,
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
        print_log("Output (.rep): " + rep_file_out, cfg.quiet)

    if cfg.remove_all:
        for af in all_files.values():
            os.remove(af)
    print_log("", cfg.quiet)
    return True


def get_top_match(matches, prob):
    target, kcount = matches[0]
    best = prob[target] * kcount
    for t, k in matches[1:]:
        score = prob[t] * k
        if score > best or (score == best and k > kcount):
            target, kcount, best = t, k, score
    return target, kcount
