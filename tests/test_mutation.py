import sys
import os
from nose.tools import ok_, assert_equals, assert_not_equal
import numpy as np
from os.path import join as ospj

FILE_PATH = os.path.realpath(__file__)
TEST_DIR_PATH = os.path.dirname(FILE_PATH)
DATA_PATH = os.path.abspath(ospj(TEST_DIR_PATH, "test_data"))
TMP_DIR_PATH = ospj(TEST_DIR_PATH, "nose_tmp_output")
TMP_BASENAME_DIR = ospj(TMP_DIR_PATH, "validation")
PKG_PATH = ospj(TEST_DIR_PATH, '..')

sys.path.append(PKG_PATH)
from sufam import mutation


class TestMutation(object):
    def test_equality_mutation(self):
        mutations = mutation.parse_vcf(ospj(DATA_PATH, "mutation_tests.vcf"))
        ok_(mutations[0] == mutations[1])
        ok_(mutations[3] == mutations[4])
        ok_(mutations[2] != mutations[3])

    def test_correct_mutation(self):
        mutations = mutation.parse_vcf(ospj(DATA_PATH, "mutation_tests.vcf"))
        assert_equals(mutations[0].type, ".")
        assert_equals(mutations[0].change, "G")
        assert_equals(mutations[2].type, "-")
        assert_equals(mutations[2].change, "A")
        assert_equals(mutations[3].type, "+")
        assert_equals(mutations[3].change, "A")

    def test_normal_filter(self):
        nmuts = mutation.parse_vcf(ospj(DATA_PATH, "normal_mutation_tests.vcf"))
        for m in nmuts:
            m.count = 21
            m.cov = 100
        nmuts = mutation.MutationsAtSinglePosition.from_mutation_list(nmuts)

        muts = mutation.parse_vcf(ospj(DATA_PATH, "mutation_tests.vcf"))
        for m in muts:
            m.count = 21
            m.cov = 100
        muts = mutation.MutationsAtSinglePosition.from_mutation_list(muts)

        filt_muts = muts.filter_against_normal(nmuts)
        assert_equals(len(filt_muts), 1)
