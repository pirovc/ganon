from ganon.util import print_log, check_file

import os
import pathlib
from collections import defaultdict

HARD_CAP_ITER = 1000


def parse_rep(rep_file_path):
    hierarchies = set()
    rep_file_totals = []
    with rep_file_path.open("r") as rep:
        for line in rep:
            if not line.strip():
                continue
            if line[0] != "#":
                hierarchies.add(line.split("\t")[0])
            else:
                rep_file_totals.append([line.rstrip()])
    return sorted(hierarchies), rep_file_totals


def define_all_files(rep_file_prefix, hierarchies, quiet):
    all_files = {}
    for h in hierarchies:
        if check_file(rep_file_prefix + "." + h + ".all"):
            all_files[h] = rep_file_prefix + "." + h + ".all"
        elif check_file(rep_file_prefix + ".all"):
            all_files = {}
            all_files[""] = rep_file_prefix + ".all"
            if len(hierarchies) > 1:
                print_log(
                    "Reassignment from concatenated single files using ganon classify --output-single and multiple hierarchies is not supported.",
                    quiet,
                )
                return None
            break
        else:
            print_log(
                "Could not find matching file for .rep [" + rep_file_prefix + "*.all]",
                quiet,
            )
            return None
    return all_files


def parse_all_file(all_file, hierarchy):
    # auto increment dict to store (hiearchy, target) as int
    targets = defaultdict(lambda: len(targets))
    read_matches = defaultdict(list)
    with open(all_file, "r") as af:
        for line in af:
            if not line.strip():
                continue
            readid, target, kcount = line.rstrip().split("\t")[:3]
            read_matches[readid].append((targets[(hierarchy, target)], int(kcount)))
    return read_matches, targets


def run_em(cfg, read_matches, targets):

    prob = {}
    # Initialization
    # initialize all targets with zero weights
    # count one with unique matches
    unique_counts = {tid: 0 for tid in targets.values()}
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

            # should not happend with EPS - here as safeguard
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
            " - Iteration " + str(em_ite_cnt + 1) + " (" + str(round(diff, 6)) + ")",
            cfg.quiet,
        )

        # If abs. difference among old and new probabilities already converged (no change)
        if diff <= cfg.threshold:
            break
        if cfg.max_iter > 0 and em_ite_cnt == cfg.max_iter - 1:
            break
        if cfg.max_iter == 0 and em_ite_cnt == HARD_CAP_ITER - 1:
            print_log(
                f"Reach max. default iterations [{HARD_CAP_ITER}]. If you need more iterations, set --max-iter to a value greater than 0."
            )
            break

        em_ite_cnt += 1

    # Final hard E-step: one integer count per read

    top_matches = {}
    reassigned_matches = unique_counts.copy()

    for read, matches in read_matches.items():
        if len(matches) > 1:
            target, kcount = get_top_match(matches, prob)
            reassigned_matches[target] += 1
            if not cfg.skip_one:
                top_matches[read] = (target, kcount)

    if sum(reassigned_matches.values()) != len(read_matches):
        raise RuntimeError("Inconsistent read redistribution")

    return reassigned_matches, top_matches


def write_one(one_file_out, read_matches, targets, top_matches):
    # reverse string target <-> integer id
    targets_rev = {val: key for (key, val) in targets.items()}
    with open(one_file_out, "w") as out_file:
        for readid, matches in read_matches.items():
            if len(matches) == 1:
                print(
                    readid,
                    targets_rev[matches[0][0]][1],
                    matches[0][1],
                    sep="\t",
                    file=out_file,
                )
            else:
                target, kcount = top_matches[readid]
                print(
                    readid,
                    targets_rev[target][1],
                    kcount,
                    sep="\t",
                    file=out_file,
                )


def generate_rep(rep_file_path, hierarchy, reassigned_matches, targets):
    rep_lines = []
    with rep_file_path.open("r") as rep:
        for line in rep:
            if not line.strip():
                continue
            if line[0] != "#":
                fields = line.rstrip().split("\t")
                hierarchy_name = fields[0]
                target = fields[1]
                direct_matches = fields[2]
                unique = int(fields[3])
                rank = fields[5] if len(fields) >= 6 else ""
                name = fields[6] if len(fields) >= 7 else ""

                # Only print line of targets
                h = hierarchy_name if hierarchy else ""
                tid = targets.get((h, target))
                if (
                    tid is not None
                    and (hierarchy == "" or hierarchy_name == hierarchy)
                    and tid in reassigned_matches
                ):
                    rep_lines.append(
                        [
                            hierarchy_name,
                            target,
                            direct_matches,
                            unique,
                            reassigned_matches[tid] - unique,
                            rank,
                            name,
                        ]
                    )
    return rep_lines


def write_rep(output_prefix, new_rep_lines, rep_file_totals):
    with open(output_prefix + ".rep", "w") as rep_out:
        for line in new_rep_lines:
            print(*line, sep="\t", file=rep_out)
        for totals in rep_file_totals:
            print(*totals, sep="\t", file=rep_out)


def reassign(cfg):

    print_log("Reassigning reads", cfg.quiet)

    if cfg.skip_one and cfg.skip_rep:
        print_log("No valid output with --skip-one and --skip-rep", cfg.quiet)
        return False

    rep_file_path = pathlib.Path(cfg.input)
    rep_file_prefix = str(pathlib.Path(rep_file_path.parent, rep_file_path.stem))

    if not check_file(rep_file_path):
        print_log(f"Report file not found [{rep_file_path}]", cfg.quiet)
        return False

    hierarchies, rep_file_totals = parse_rep(rep_file_path)

    if not hierarchies:
        print_log(f"No valid report entries to reassign [{cfg.input}]", cfg.quiet)
        return False

    all_files = define_all_files(rep_file_prefix, hierarchies, cfg.quiet)

    if all_files is None:
        return False

    if cfg.output_prefix:
        output_prefix = cfg.output_prefix
    else:
        output_prefix = rep_file_prefix

    new_rep_lines = []
    for hierarchy, all_file in all_files.items():
        print_log(all_file + (" [" + hierarchy + "]" if hierarchy else ""), cfg.quiet)

        read_matches, targets = parse_all_file(all_file, hierarchy)

        reassigned_matches, top_matches = run_em(cfg, read_matches, targets)

        # Skip .one file (just generate .rep)
        if not cfg.skip_one:
            if len(all_files) == 1:
                one_file_out = output_prefix + ".one"
            else:
                one_file_out = output_prefix + "." + hierarchy + ".one"
            write_one(one_file_out, read_matches, targets, top_matches)
            print_log(
                "Output (.one): " + one_file_out,
                cfg.quiet,
            )

        if not cfg.skip_rep:
            new_rep_lines.extend(
                generate_rep(rep_file_path, hierarchy, reassigned_matches, targets)
            )

    if not cfg.skip_rep:
        write_rep(output_prefix, new_rep_lines, rep_file_totals)
        print_log("Output (.rep): " + output_prefix + ".rep", cfg.quiet)

    if cfg.remove_all:
        for all_file in all_files.values():
            os.remove(all_file)
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
