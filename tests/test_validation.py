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
import sufam.__main__
import sufam.mpileup_parser as mpileup_parser
from sufam import utils


class TestValidation(object):
    def setUp(self):
        """Delete temporary dir if it exists then create it"""
        self.tearDown()
        utils.mkdir_p(TMP_BASENAME_DIR)

    def tearDown(self):
        """remove temp output files"""
        utils.rm_rf(TMP_DIR_PATH)

    def test_validate_mutations(self):
        out = StringIO()
        sufam.__main__.validate_mutations(ospj(DATA_PATH, "mutations.vcf"),
                                            ospj(DATA_PATH, "OCT10T.bam"),
                                            ospj(DATA_PATH, "human_g1k_v37.fa"),
                                            "test",
                                            "matrix",
                                            out
                                            )
        assert_equals("0\n1\n", out.getvalue())

    def test_validate_mutations_indel(self):
        out = StringIO()
        sufam.__main__.validate_mutations(ospj(DATA_PATH, "mutations_indel.vcf"),
                                            ospj(DATA_PATH, "OCT9T.bam"),
                                            ospj(DATA_PATH, "human_g1k_v37.fa"),
                                            "test",
                                            "matrix",
                                            out
                                            )
        assert_equals("1\n0\n1\n0\n0\n", out.getvalue())

    def test_mpileup_parser_two_digit_indel(self):
        two_digit_indel = "X\t150349557\tC\t24\t.$,-12caccactggcca.-12CACCACTGGCCA.,.,,,,.-12CACCACTGGCCA,..,,-12caccactggcca..-12CACCACTGGCCA,,,..,\t;FCDDDDDDD/FDCC/C/E<FBDC\n"
        assert_equals("X\t150349557\tC\t24\t0\t19\t0\t0\t0\tCACCACTGGCCA,CACCACTGGCCA,CACCACTGGCCA,CACCACTGGCCA,CACCACTGGCCA\t\t",
            mpileup_parser.parse(two_digit_indel))
