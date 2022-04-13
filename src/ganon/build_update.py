import time
import pandas as pd

from ganon.util import validate_input_files
from ganon.util import print_log
from ganon.util import set_tmp_folder
from ganon.util import run
from io import StringIO

from multitax import NcbiTx, GtdbTx

def build(cfg):
    return True

def update(cfg):
    return True

def build_custom(cfg):

    # Retrieve and check input files or folders
    input_files = validate_input_files(cfg.input, cfg.input_extension, cfg.quiet)
    if not input_files:
        print_log("ERROR: No valid input files found")
        return False

    # Set --input-target if not manually set
    if not cfg.input_target:
        cfg.input_target = "sequence" if len(input_files)==1 else "file"

    # Set working folder 
    tmp_output_folder = cfg.db_prefix + "_tmp/"
    if not set_tmp_folder(tmp_output_folder, cfg.restart): return False

    # Set-up taxonomy
    if cfg.taxonomy:
        tx = time.time()
        print_log("Parsing " + cfg.taxonomy + " taxonomy", cfg.quiet)
        if cfg.taxonomy=="ncbi":
            tax = NcbiTx(files=cfg.taxonomy_files)
        elif cfg.taxonomy=="gtdb":
            tax = GtdbTx(files=cfg.taxonomy_files)
        print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)
    else:
        tax = None

    # Set-up input info
    info = load_target_info(cfg, input_files)


    # # If taxonomy is given and level is an taxonomic rank, convert nodes
    # if tax:
    #     # get latest tax.latest()
    #     # validate node on given tax (if came from target_info for example and it's outdated or wrong)
    #     if cfg.level!="leaves" and cfg.level in set(tax._ranks):
    #         info = replace_node_rank()

    # else:
    #     # na on col nodes


    # replaced_spec = validate_specialization(info)
    # if replaced_spec:
    #     print_log(str(replaced_spec) + " invalid specialization entries replaced by target\n", cfg.quiet)


    # load input-info
    #or
    # parse input files
    # if --taxonomy: retrieve input info

    # save-state 1

    # filter tax
    # validate specialization
    # write tax

    # write target-info file

    # run ganon-build

    return True

def update_custom(cfg):
    return True


def load_target_info(cfg, input_files):

    target_info_colums=['target','node','specialization']

    tx = time.time()
    if cfg.target_info:
        print_log("Parsing --target-info " + cfg.target_info, cfg.quiet)
        info = pd.read_csv(cfg.target_info, sep='\t', header=None, skiprows=0, index_col=None, dtype='str', names=target_info_colums)
        # Drop cols without values
        info.dropna(how="all", inplace=True)
        # Drop duplicated target (seqid or fileid)
        info.drop_duplicates(subset=['target'], inplace=True)
        # set target as index
        info.set_index('target', inplace=True)
        print_log(" - " + str(info.size) + " unique entries", cfg.quiet)

    else:
        info = pd.DataFrame(columns=target_info_colums)
        print(info)
        if cfg.input_target=="sequence":
            sequence_accessions = parse_sequence_accession(input_files, info)
            print(sequence_accessions)
            # return also dict {acc.ver: file}
            #eutils, accession2taxid...
            # if cfg.level == "assembly": get also assembly
            # do not need to have tax
            info = retrieve_tax_sequence(cfg, tax, targets) 

        # elif cfg.input_target=="file":
        #     targets = input_files
        #     # return file names
        #     info = retrieve_tax_file(cfg, tax, targets) # new - get from assembly_summary
        

    print_log(" - done in " + str("%.2f" % (time.time() - tx)) + "s.\n", cfg.quiet)

    return info

def validate_specialization(info):
    """
    validate specializatio, if given
    each specialization can have only one parent node
    return number of replaced specializations
    """

    # check for invalid specialization entries
    idx_null_spec = info.specialization.isnull()
    # if all entries are null, no specialization was given
    if all(idx_null_spec):
        return 0
    # get unique tuples node-specialization
    node_spec = info[['node', 'specialization']].drop_duplicates()
    # check for duplicated specialization in the tuples
    idx_multi_parent_spec = info.specialization.isin(node_spec.specialization[node_spec.specialization.duplicated(keep=False)].unique())
    # merge indices for invalid entries
    idx_replace = idx_null_spec | idx_multi_parent_spec
    if idx_replace.any():
        # replace invalid specialization entries with target
        info.loc[idx_replace,"specialization"] = info.index[idx_replace]
        return sum(idx_replace)
    return 0

def parse_sequence_accession(input_files, info):
    sequence_accessions = pd.DataFrame(columns=["target","file"])

    for file in input_files:
        # cat | zcat | gawk -> compability with osx
        run_get = "cat {0} {1} | grep -o '^>[^ ]*' | sed 's/>//'".format(file, "| zcat" if file.endswith(".gz") else "")
        stdout, stderr = run(run_get, shell=True)
        #sequence_accessions.update(stdout.split("\n"))
        x = pd.read_csv(StringIO(stdout), header=None, names=['target'])
        x["file"] = file
        sequence_accessions = pd.concat([sequence_accessions, x])
        
    return sequence_accessions