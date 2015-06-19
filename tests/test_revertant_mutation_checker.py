from cStringIO import StringIO
import sys
import os
from nose.tools import ok_, assert_equals
import numpy as np
from os.path import join as ospj

FILE_PATH = os.path.realpath(__file__)
TEST_DIR_PATH = os.path.dirname(FILE_PATH)
DATA_PATH = os.path.abspath(ospj(TEST_DIR_PATH, "test_data"))
TMP_DIR_PATH = ospj(TEST_DIR_PATH, "nose_tmp_output")
TMP_BASENAME_DIR = ospj(TMP_DIR_PATH, "validation")
PKG_PATH = ospj(TEST_DIR_PATH, '..')

sys.path.append(PKG_PATH)
from sufam import revertant_mutation_checker
from sufam import utils

from Bio import SeqIO
import pyhgvs as hgvs


class TestRevertantMutation(object):
    def setUp(self):
        """Delete temporary dir if it exists then create it"""
        self.tearDown()
        utils.mkdir_p(TMP_BASENAME_DIR)

    def tearDown(self):
        """remove temp output files"""
        utils.rm_rf(TMP_DIR_PATH)

    def test_apply_hgvs(self):
        # p.Leu1303Phefs == c.3908dupT
        transcripts = \
            SeqIO.to_dict(SeqIO.parse("tests/test_data/BRCA1_transcripts.fa",
            "fasta"))
        brca1_mut = hgvs.HGVSName("ENST00000357654:c.3908dupT")
        normal_p = transcripts["ENST00000357654"].seq.translate()
        assert_equals("L", normal_p[1302])
        mut_c = revertant_mutation_checker.apply_hgvs(transcripts["ENST00000357654"], brca1_mut)
        assert_equals("TT", mut_c[3907:3909])
        mut_p = mut_c.translate()
        assert_equals("F", mut_p[1302])

    def test_revertant_mutation_checker(self):
        out = StringIO()
        revertant_mutation_checker.check_revertant_mutations(
            ospj(DATA_PATH, "to_be_reverted_mutations.txt"),
            ospj(DATA_PATH, "oncotator_del_maf.txt"),
            ospj(DATA_PATH, "BRCA1_transcripts.fa")
        )
        assert_equals("1\n0\n1\n0\n0\n", out.getvalue())
