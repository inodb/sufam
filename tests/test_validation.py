from cStringIO import StringIO
import sys
import os
from nose.tools import ok_, assert_equals, assert_almost_equals
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
        assert_equals("X\t150349557\tC\t24\t0\t24\t0\t0\t0\tCACCACTGGCCA,CACCACTGGCCA,CACCACTGGCCA,CACCACTGGCCA,CACCACTGGCCA\t",
            mpileup_parser.parse(two_digit_indel))

    def test_mpileup_parser_indel_followed_by_snv(self):
        indel_followed_by_snv = "X\t150349557\tC\t28\t.$,-12caccactggccat-12CACCACTGGCCAT.,.,,,,.-12CACCACTGGCCAT,..,,-12caccactggccaG..-12CACCACTGGCCAA,,,..,\t;FCDDDDDDD/FDCC/C/E<FBDC\n"
        assert_equals("X\t150349557\tC\t28\t1\t23\t1\t3\t0\tCACCACTGGCCA,CACCACTGGCCA,CACCACTGGCCA,CACCACTGGCCA,CACCACTGGCCA\t",
            mpileup_parser.parse(indel_followed_by_snv))

    def test_mpileup_test1(self):
        test = open(ospj(DATA_PATH, "mpileup_test1.tsv")).read()
        bpdf = sufam.__main__.get_baseparser_extended_df("test", [mpileup_parser.parse(test)], "AG", "A")
        assert_equals(test.count(",") + test.count("."), int(bpdf['cov'].iloc[0]))
        assert_equals(int(bpdf['cov'].iloc[0]), int(bpdf.A.iloc[0]))
        assert_almost_equals(0.7225, float(bpdf.val_maf.iloc[0]), places=3)
        assert_almost_equals(0.7225, float(bpdf.most_common_indel_maf.iloc[0]), places=3)
        assert_equals("-", bpdf.most_common_indel_type.iloc[0])

    def test_mpileup_test2(self):
        test = open(ospj(DATA_PATH, "mpileup_test2.tsv")).read()
        bpdf = sufam.__main__.get_baseparser_extended_df("test", [mpileup_parser.parse(test)], "G", "GAA")
        assert_equals(int(bpdf['cov'].iloc[0]), int(bpdf.G.iloc[0]))
        assert_equals(test.count(",") + test.count("."), int(bpdf['cov'].iloc[0]))
        assert_almost_equals(0.4324, float(bpdf.val_maf.iloc[0]), places=3)
        assert_almost_equals(0.4324, float(bpdf.most_common_indel_maf.iloc[0]), places=3)
        assert_equals("+", bpdf.most_common_indel_type.iloc[0])

    def test_mpileup_test3(self):
        test = open(ospj(DATA_PATH, "mpileup_test3.tsv")).read()
        bpdf = sufam.__main__.get_baseparser_extended_df("test", [mpileup_parser.parse(test)], "G", "A")
        assert_equals(int(bpdf['cov'].iloc[0]), int(bpdf.G.iloc[0]) + int(bpdf.A.iloc[0]) + int(bpdf["T"].iloc[0]))
        assert_equals(1, int(bpdf["T"].iloc[0]))
        assert_equals("AA", bpdf.most_common_indel.iloc[0])
        assert_equals("+", bpdf.most_common_indel_type.iloc[0])
        assert_almost_equals(0.0139, float(bpdf.val_maf.iloc[0]), places=3)
        assert_almost_equals(0.0139, float(bpdf.most_common_al_maf.iloc[0]), places=3)
