#!/usr/bin/env python3
import pickle
import gzip


class Gnn:
    def __init__(self,
                 kmer_size: int=None,
                 window_size: int=None,
                 hash_functions: int=None,
                 number_of_bins: int=None,
                 bin_length: int=None,
                 fragment_length: int=None,
                 overlap_length: int=None,
                 rank: str=None,
                 specialization: str=None,
                 bins: list=None,
                 file: str=None):
        if file:
            self.parse(file)
        else:
            self.kmer_size = kmer_size
            self.window_size = window_size
            self.hash_functions = hash_functions
            self.number_of_bins = number_of_bins
            self.rank = rank
            self.specialization = specialization
            self.bin_length = bin_length
            self.fragment_length = fragment_length
            self.overlap_length = overlap_length
            self.bins = bins

    def __repr__(self):
        args = ['{}={}'.format(k, repr(v)) for (k, v) in vars(self).items()]
        return 'Gnn({})'.format(', '.join(args))

    def parse(self, file):
        gnn = pickle.load(gzip.open(file, "rb"))
        self.kmer_size = gnn['kmer_size']
        self.window_size = gnn['window_size']
        self.hash_functions = gnn['hash_functions']
        self.number_of_bins = gnn['number_of_bins']
        self.rank = gnn['rank']
        self.specialization = gnn['specialization']
        self.bin_length = gnn['bin_length']
        self.fragment_length = gnn['fragment_length']
        self.overlap_length = gnn['overlap_length']
        self.bins = gnn['bins']

    def write(self, file):
        gnn = {'kmer_size': self.kmer_size,
               'window_size': self.window_size,
               'hash_functions': self.hash_functions,
               'number_of_bins': self.number_of_bins,
               'rank': self.rank,
               'specialization': self.specialization,
               'bin_length': self.bin_length,
               'fragment_length': self.fragment_length,
               'overlap_length': self.overlap_length,
               'bins': self.bins}
        with gzip.open(file, 'wb') as f:
            pickle.dump(gnn, f)
