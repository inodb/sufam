"""
Functions for oncotator
"""
import pandas as pd
import pyhgvs as hgvs
import sys


class Oncotator(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.df = pd.read_csv(filepath, sep="\t", comment="#")

    def get_hgvs_mutations(self, transcript_id, ignore_introns=True):
        hgvs_muts = []
        if ignore_introns:
            df = self.df[self.df.Variant_Classification != "Intron"]
        else:
            df = self.df
        for r in df.HGVS_coding_DNA_change:
            if r != "Exception_encountered":
                try:
                    h = hgvs.HGVSName(r)
                    # remove .version postfix
                    if h.transcript == transcript_id or ".".join(h.transcript.split(".")[:-1]) == transcript_id:
                        hgvs_muts += [h]
                except hgvs.InvalidHGVSName:
                    sys.stderr.write("Invalid HGVS found: {}\n".format(r))
                    pass
        return hgvs_muts
