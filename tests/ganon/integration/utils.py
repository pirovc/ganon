import shutil, os
import pandas as pd
from pathlib import Path

def setup_dir(results_dir):
    shutil.rmtree(results_dir, ignore_errors=True)
    os.makedirs(results_dir)

def check_files(prefix, extensions):
    for ext in extensions:
        f=prefix+"."+ext if ext else prefix
        if not Path(f).is_file():
            print("File (" + f +") was not created")
            return False
        elif Path(f).stat().st_size==0:
            print("File (" + f +") is empty")
            return False
    return True
    
def parse_seq_info(seq_info_file):
    colums=['seqid', 'length', 'taxid', 'specialization']
    types={'seqid': 'str', 'length': 'uint64', 'taxid': 'str', 'specialization': 'str'}
    return pd.read_table(seq_info_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_map(map_file):
    colums=['target', 'binid']
    types={'target': 'str', 'binid': 'uint64'}
    return pd.read_table(map_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_tax(tax_file):
    colums=['taxid', 'parent', 'rank', 'name']
    types={'taxid': 'str', 'parent': 'str', 'rank': 'str', 'name': 'str'}
    return pd.read_table(tax_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_tre(tre_file):
    colums=['rank', 'target', 'lineage', 'name', 'unique', 'total', 'cumulative', 'cumulative_perc']
    types={'rank': 'str', 'target': 'str', 'lineage': 'str', 'name': 'str', 'unique': 'str', 'total': 'str', 'cumulative': 'str', 'cumulative_perc': 'str'}
    return pd.read_table(tre_file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)
