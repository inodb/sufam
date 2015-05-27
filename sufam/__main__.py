#!/usr/bin/env python
"""
So U Found A Mutation? (SUFAM)

Found a mutation in one or more samples? Now you want to check if they are in
another sample. Unfortunately mutect, varscan or whatever other variant caller
is not calling them. Use SUFAM. The super sensitive validation caller that
calls everything on a given position. All you need is a vcf with the mutations
that you are interested in and the sam/bam file of the sample where you want to
find the same inconsipicuous mutation.

Author: inodb
"""
import argparse
from collections import Counter
import glob
import subprocess
import sys

import numpy as np
import pandas as pd

from sufam import mpileup_parser


def get_pile_up_baseparser(bam, chrom, pos1, pos2, reffa):
    """Get vcf line of given chrom and pos from mpileup and baseparser"""
    posmin = min(pos1, pos2)
    posmax = max(pos1, pos2)
    cmd = "samtools view -bh {bam} {chrom}:{pos1}-{pos2} " \
            "| samtools mpileup -R -q 1 -f {reffa} -".format(bam=bam,chrom=chrom,pos1=posmin,pos2=posmax, reffa=reffa)
    if pos1 == pos2:
        cmd += " | grep -P '^{chrom}\\t{pos}\\t'".format(chrom=chrom,pos=pos1)
    else:
        cmd += " | tail -n +2"
    sys.stderr.write("Running:\n{}\n".format(cmd))
    rv = subprocess.Popen(cmd,
            shell=True, stdout=subprocess.PIPE).communicate()[0]
    return [mpileup_parser.parse(line) for line in rv.split("\n")[:-1]]


def _most_common_al(x):
    if x["ref"]:
        bl = ["A","C","G","T"]
        bl.remove(x["ref"])
        mc = Counter({k:int(v) for k, v in dict(x.ix[bl]).iteritems()}).most_common(1)[0]
        return pd.Series({"most_common_al":str(mc[0]),
                        "most_common_al_count":str(mc[1]),
                        "most_common_al_maf":str(mc[1]/float(x["cov"])) if float(x["cov"]) > 0 else "0"})
    else:
        return pd.Series({"most_common_al":None,
                        "most_common_al_count":None,
                        "most_common_al_maf":None})


def _most_common_indel(x):
    dels = Counter(x["-"].split(",")).most_common(1)[0] if x["-"] else None
    ins = Counter(x["+"].split(",")).most_common(1)[0] if x["+"] else None
    if ins and dels:
        mc = dels if dels[1] >= ins[1] else ins
        indel_type = "-" if dels[1] >= ins[1] else "+"
    elif ins:
        mc = ins
        indel_type = "+"
    elif dels:
        mc = dels
        indel_type = "-"
    else:
        return pd.Series({
                        "most_common_indel_type":None,
                        "most_common_indel":None,
                        "most_common_indel_count":None,
                        "most_common_indel_maf":None})
    return pd.Series({
                    "most_common_indel_type":indel_type,
                    "most_common_indel":str(mc[0]),
                    "most_common_indel_count":str(mc[1]),
                    "most_common_indel_maf":str(mc[1]/float(x["cov"]) if float(x["cov"]) > 0 else "0")})


def get_baseparser_extended_df(bam, sample, chrom, pos1, pos2, reffa):
    """Turn baseParser results into a dataframe"""
    columns = "chrom\tpos\tref\tcov\tA\tC\tG\tT\t*\t-\t+\tX".split()
    bp_lines = get_pile_up_baseparser(bam, chrom, pos1, pos2, reffa)

    # change baseparser output to get most common maf per indel
    bpdf = pd.DataFrame([[bam.split("/")[1][:-4]] + l.rstrip('\n').split("\t") for l in bp_lines if len(l) > 0],
                        columns=["sample"] + columns, dtype=np.object)

    if len(bpdf) == 0:
        return None

    bpdf = pd.concat([bpdf, bpdf.apply(_most_common_indel, axis=1)], axis=1)
    bpdf = pd.concat([bpdf, bpdf.apply(_most_common_al, axis=1)], axis=1)
    bpdf["most_common_count"] = bpdf.apply(lambda x: max([x.most_common_al_count, x.most_common_indel_count]), axis=1)
    bpdf["most_common_maf"] = bpdf.apply(lambda x: max([x.most_common_al_maf, x.most_common_indel_maf]), axis=1)

    return bpdf


def filter_out_mutations_in_normal(tumordf, normaldf, most_common_maf_min=0.2,
                                   most_common_count_maf_threshold=20,
                                   most_common_count_min=1):
    """Remove mutations that are in normal"""
    df = tumordf.merge(normaldf, on=["chrom","pos"], suffixes=("_T","_N"))

    # filters
    common_al = (df.most_common_al_count_T == df.most_common_count_T) & (df.most_common_al_T == df.most_common_al_N)
    common_indel = (df.most_common_indel_count_T == df.most_common_count_T) & (df.most_common_indel_T == df.most_common_indel_N)
    normal_criteria = ((df.most_common_count_N >= most_common_count_maf_threshold) & (df.most_common_maf_N > most_common_maf_min)) | \
                      ((df.most_common_count_N < most_common_count_maf_threshold) & (df.most_common_count_N > most_common_count_min))
    df = df[~(common_al | common_indel) & normal_criteria]

    # restore column names of tumor
    for c in df.columns:
        if c.endswith("_N"):
            del df[c]
    df.columns = [c[:-2] if c.endswith("_T") else c for c in df.columns]

    return df


def select_only_revertant_mutations(bpdf, snv=None, ins=None, dlt=None):
    """
    Selects only mutations that revert the given mutations in a single event.
    """
    if sum([bool(snv), bool(ins), bool(dlt)]) != 1:
        raise(Exception("Should be either snv, ins or del".format(snv)))

    if snv:
        if snv not in ["A","C","G","T"]:
            raise(Exception("snv {} should be A, C, G or T".format(snv)))
            return bpdf[(bpdf.most_common_al == snv) & (bpdf.most_common_al_count == bpdf.most_common_count)]
    elif bool(ins):
        return \
            bpdf[((bpdf.most_common_indel.apply(lambda x: len(x) + len(ins) % 3 if x else None) == 0 ) & (bpdf.most_common_indel_type == "+") & (bpdf.most_common_count == bpdf.most_common_indel_count)) |
                    ((bpdf.most_common_indel.apply(lambda x: len(ins) - len(x) % 3 if x else None) == 0 ) & (bpdf.most_common_indel_type == "-") & (bpdf.most_common_count == bpdf.most_common_indel_count))]
    elif bool(dlt):
        return \
            bpdf[((bpdf.most_common_indel.apply(lambda x: len(x) - len(dlt) % 3 if x else None) == 0) & (bpdf.most_common_indel_type == "+") & (bpdf.most_common_count == bpdf.most_common_indel_count)) |
                    ((bpdf.most_common_indel.apply(lambda x: -len(dlt) - len(x) % 3 if x else None) == 0 ) & (bpdf.most_common_indel_type == "-") & (bpdf.most_common_count == bpdf.most_common_indel_count))]
    else:
        # should never happen
        raise(Exception("No mutation given?"))


def validate_mutations(vcffile, bam, reffa, sample):
    """Check if mutations in vcf are in bam"""
    header = []
    row = []
    for line in open(vcffile):
        if line.startswith("#CHROM"):
            header = line[1:].rstrip('\n').split("\t")
        if line.startswith("#"):
            continue
        record = dict(zip(header, line.rstrip('\n').split("\t")))
        bpdf = get_baseparser_extended_df(bam, sample, record["CHROM"], record["POS"], record["POS"], reffa)
        if len(record["REF"]) == len(record["ALT"]):
            if len(bpdf[(bpdf.most_common_al == record["ALT"]) & (bpdf.most_common_al_maf > 0)]) > 0:
                row += [True]
            else:
                row += [False]
        else:
            # deletion
            if len(record["REF"]) > len(record["ALT"]):
                if len(bpdf[(bpdf.most_common_indel == record["REF"][1:]) & (bpdf.most_common_indel_maf > 0)]) > 0:
                    row += [True]
                else:
                    row += [False]
            # insertion
            else:
                if len(bpdf[(bpdf.most_common_indel == record["ALT"][1:]) & (bpdf.most_common_indel_maf > 0)]) > 0:
                    row += [True]
                else:
                    row += [False]
    return "\t".join([str(int(x)) for x in row])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reffa", type=str, help="Reference genome (fasta)")
    parser.add_argument("vcf", type=str, help="VCF with mutations to be validated")
    parser.add_argument("bam", type=str, help="BAM to find mutations in")
    parser.add_argument("--sample_name", type=str, default=None, help="Set name "
                        "of sample, used in output.")
    args = parser.parse_args()
    if args.sample_name is None:
        args.sample_name = args.bam
    print validate_mutations(args.vcf, args.bam, args.reffa, args.sample_name)


if __name__ == "__main__":
    main()
