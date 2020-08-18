#!/usr/bin/env python3
import pandas as pd
from io import StringIO
from src.ganon.util import *

class SeqInfo:
    seq_info_colums=['seqid','length','taxid', 'assembly']
    seq_info_types={'seqid': 'str', 'length': 'uint64', 'taxid': 'str', 'assembly': 'str'}
    seqinfo = pd.DataFrame(columns=seq_info_colums)

    def __init__(self, seq_info_file: str=None):
        if seq_info_file:
            self.seqinfo = pd.read_csv(seq_info_file, sep='\t', header=None, skiprows=0, names=self.seq_info_colums, dtype=self.seq_info_types)

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k,v) in vars(self).items()]
        return 'SeqInfo({})'.format(', '.join(args))
    
    def get_csv(self):
        return self.seqinfo.to_csv(sep="\t",header=False, index=False)

    def size(self):
        return self.seqinfo.shape[0]

    def get_seqids(self):
        return self.seqinfo.seqid

    def remove_seqids(self, seqids):
        self.seqinfo = self.seqinfo[~self.seqinfo['seqid'].isin(seqids)]

    def write_seqid_file(self, seqid_file):
        self.seqinfo["seqid"].to_csv(seqid_file, header=False, index=False)

    def parse_seqid(self, input_files):
        for file in input_files:
            # cat | zcat | gawk -> compability with osx
            run_get_header = "cat {0} {1} | gawk 'BEGIN{{FS=\" \"}} /^>/ {{print substr($1,2)}}'".format(file, "| zcat" if file.endswith(".gz") else "")
            stdout, stderr = run(run_get_header, print_stderr=False, shell=True)
            self.seqinfo = self.seqinfo.append(pd.read_csv(StringIO(stdout), header=None, names=['seqid']), ignore_index=True)

    def parse_seqid_length(self, input_files):
        for file in input_files:
            # cat | zcat | gawk -> compability with osx
            run_get_length = "cat {0} {1} | gawk 'BEGIN{{FS=\" \"}} /^>/ {{if (seqlen){{print seqlen}}; printf substr($1,2)\"\\t\";seqlen=0;next;}} {{seqlen+=length($0)}}END{{print seqlen}}'".format(file, "| zcat" if file.endswith(".gz") else "")
            stdout, stderr = run(run_get_length, print_stderr=False, shell=True)
            self.seqinfo = self.seqinfo.append(pd.read_csv(StringIO(stdout), sep="\t", header=None, names=['seqid', 'length']), ignore_index=True)
            # remove entries with length 0
            self.seqinfo = self.seqinfo[self.seqinfo['length']>0]

    def parse_ncbi_eutils(self, seqid_file, path_exec_get_len_taxid, skip_len_taxid=False, get_assembly=False):
        run_get_len_taxid_cmd = '{0} -i {1} {2} {3} -r'.format(
                                    path_exec_get_len_taxid,
                                    seqid_file,
                                    "-a" if get_assembly else "",
                                    "-s" if skip_len_taxid else "")
        stdout, stderr = run(run_get_len_taxid_cmd, print_stderr=True, exit_on_error=False)
        
        if get_assembly and skip_len_taxid:
            # todo - order is always the same?
            self.seqinfo.assembly = pd.read_csv(StringIO(stdout), sep='\t', header=None, skiprows=0, names=['seqid','assembly'], dtype=self.seq_info_types)['assembly']
        else:
            self.seqinfo = pd.read_csv(StringIO(stdout), sep='\t', header=None, skiprows=0, names=self.seq_info_colums, dtype=self.seq_info_types)

    def parse_acc2txid(self, acc2txid_files):
        count_acc2txid = {}
        set_seqids = set(self.seqinfo.seqid)
        for acc2txid in acc2txid_files:
            tmp_seqid_taxids = pd.read_csv(acc2txid, sep='\t', header=None, skiprows=1, usecols=[1,2], names=['seqid','taxid'], converters={'seqid':lambda x: x if x in set_seqids else ""}, dtype={'taxid': 'str'})
            tmp_seqid_taxids = tmp_seqid_taxids[tmp_seqid_taxids['seqid']!=""] #keep only seqids used
            tmp_seqid_taxids = tmp_seqid_taxids[tmp_seqid_taxids['taxid']!="0"] # filter out taxid==0
            # save count to return
            count_acc2txid[acc2txid] = tmp_seqid_taxids.shape[0]
            # append to main file (will match col names)
            self.seqinfo = self.seqinfo.append(tmp_seqid_taxids, sort=False, ignore_index=True)
            del tmp_seqid_taxids
            #if already found all taxid for the seqids (no need to parse all files till the end)
            if sum(count_acc2txid.values()) == self.size(): 
                break 
        return count_acc2txid