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
import subprocess
import sys
import warnings

import numpy as np
import pandas as pd

from sufam import mpileup_parser


def get_pile_up_baseparser(bam, chrom, pos1, pos2, reffa):
    """Get vcf line of given chrom and pos from mpileup and baseparser"""
    posmin = min(pos1, pos2)
    posmax = max(pos1, pos2)
    cmd = "samtools view -bh {bam} {chrom}:{pos1}-{pos2} " \
        "| samtools mpileup -R -q 1 -f {reffa} -".format(bam=bam, chrom=chrom, pos1=posmin, pos2=posmax, reffa=reffa)
    if pos1 == pos2:
        cmd += " | grep -P '^{chrom}\\t{pos}\\t'".format(chrom=chrom, pos=pos1)
    else:
        cmd += " | tail -n +2"
    sys.stderr.write("Running:\n{}\n".format(cmd))
    child = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    stdout, stderr = child.communicate()
    if child.returncode != 0:
        if len(stdout) == 0 and stderr is None:
            warnings.warn("Command:\n{cmd}\n did not exit with zero exit code. "
                          "Possibly no coverage for sample.".format(cmd=cmd))
        else:
            raise(Exception("Command:\n{cmd}\n did not exit with zero exit code. "
                            "Check command.".format(cmd=cmd)))
    else:
        return [mpileup_parser.parse(line) for line in stdout.split("\n")[:-1]]


def _most_common_al(x):
    if x["ref"]:
        bl = ["A", "C", "G", "T"]
        bl.remove(x["ref"])
        mc = Counter({k: int(v) for k, v in dict(x.ix[bl]).iteritems()}).most_common(1)[0]
        return pd.Series({"most_common_al": str(mc[0]),
                          "most_common_al_count": str(mc[1]),
                          "most_common_al_maf": str(mc[1]/float(x["cov"])) if float(x["cov"]) > 0 else "0"})
    else:
        return pd.Series({"most_common_al": None,
                          "most_common_al_count": None,
                          "most_common_al_maf": None})


def _val_al(x):
    if x["ref"]:
        if len(x["val_ref"]) == 1 and len(x["val_alt"]) == 1:  # SNV
            al_count = int(x.ix[x["val_alt"]])
            al_type = "snv"
        elif len(x["val_alt"]) > len(x["val_ref"]):  # insertion
            query = x["val_alt"][len(x["val_ref"]):]
            al_count = Counter(x["+"].split(","))[query] if x["+"] is not None else 0
            al_type = "insertion"
        else:  # deletion
            query = x["val_ref"][len(x["val_alt"]):]
            al_count = Counter(x["-"].split(","))[query] if x["-"] is not None else 0
            al_type = "deletion"
        al_maf = al_count / float(x["cov"]) if float(x["cov"]) > 0 else None
        return pd.Series({"val_al_type": al_type,
                          "val_al_count": al_count,
                          "val_maf": al_maf})
    else:
        return pd.Series({"val_al_type": None,
                          "val_al_count": None,
                          "val_maf": None})


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
            "most_common_indel_type": None,
            "most_common_indel": None,
            "most_common_indel_count": None,
            "most_common_indel_maf": None})
    return pd.Series({
        "most_common_indel_type": indel_type,
        "most_common_indel": str(mc[0]),
        "most_common_indel_count": str(mc[1]),
        "most_common_indel_maf": str(mc[1]/float(x["cov"]) if float(x["cov"]) > 0 else "0")})


def get_baseparser_extended_df(bam, sample, chrom, pos1, pos2, reffa, ref, alt):
    """Turn baseParser results into a dataframe"""
    columns = "chrom\tpos\tref\tcov\tA\tC\tG\tT\t*\t-\t+".split()
    bp_lines = get_pile_up_baseparser(bam, chrom, pos1, pos2, reffa)
    if bp_lines is None:
        return None

    # change baseparser output to get most common maf per indel
    bpdf = pd.DataFrame([[bam.split("/")[1][:-4]] + l.rstrip('\n').split("\t") for l in bp_lines if len(l) > 0],
                        columns=["sample"] + columns, dtype=np.object)
    bpdf[bpdf == ""] = None

    if len(bpdf) == 0:
        return None

    # add columns for validation allele
    bpdf = pd.concat([bpdf, pd.DataFrame({"val_ref": pd.Series(ref), "val_alt": pd.Series(alt)})], axis=1)
    bpdf = pd.concat([bpdf, bpdf.apply(_val_al, axis=1)], axis=1)

    bpdf = pd.concat([bpdf, bpdf.apply(_most_common_indel, axis=1)], axis=1)
    bpdf = pd.concat([bpdf, bpdf.apply(_most_common_al, axis=1)], axis=1)
    bpdf["most_common_count"] = bpdf.apply(lambda x: max([x.most_common_al_count, x.most_common_indel_count]), axis=1)
    bpdf["most_common_maf"] = bpdf.apply(lambda x: max([x.most_common_al_maf, x.most_common_indel_maf]), axis=1)

    return bpdf


def filter_out_mutations_in_normal(tumordf, normaldf, most_common_maf_min=0.2,
                                   most_common_count_maf_threshold=20,
                                   most_common_count_min=1):
    """Remove mutations that are in normal"""
    df = tumordf.merge(normaldf, on=["chrom", "pos"], suffixes=("_T", "_N"))

    # filters
    common_al = (df.most_common_al_count_T == df.most_common_count_T) & (df.most_common_al_T == df.most_common_al_N)
    common_indel = (df.most_common_indel_count_T == df.most_common_count_T) & \
        (df.most_common_indel_T == df.imost_common_indel_N)
    normal_criteria = ((df.most_common_count_N >= most_common_count_maf_threshold) &
                       (df.most_common_maf_N > most_common_maf_min)) | \
        ((df.most_common_count_N < most_common_count_maf_threshold) &
         (df.most_common_count_N > most_common_count_min))
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
        if snv not in ["A", "C", "G", "T"]:
            raise(Exception("snv {} should be A, C, G or T".format(snv)))
        return bpdf[(bpdf.most_common_al == snv) & (bpdf.most_common_al_count == bpdf.most_common_count)]
    elif bool(ins):
        return \
            bpdf[((bpdf.most_common_indel.apply(lambda x: len(x) + len(ins) % 3 if x else None) == 0) &
                  (bpdf.most_common_indel_type == "+") & (bpdf.most_common_count == bpdf.most_common_indel_count)) |
                 ((bpdf.most_common_indel.apply(lambda x: len(ins) - len(x) % 3 if x else None) == 0) &
                  (bpdf.most_common_indel_type == "-") & (bpdf.most_common_count == bpdf.most_common_indel_count))]
    elif bool(dlt):
        return \
            bpdf[((bpdf.most_common_indel.apply(lambda x: len(x) - len(dlt) % 3 if x else None) == 0) &
                  (bpdf.most_common_indel_type == "+") & (bpdf.most_common_count == bpdf.most_common_indel_count)) |
                 ((bpdf.most_common_indel.apply(lambda x: -len(dlt) - len(x) % 3 if x else None) == 0) &
                  (bpdf.most_common_indel_type == "-") & (bpdf.most_common_count == bpdf.most_common_indel_count))]
    else:
        # should never happen
        raise(Exception("No mutation given?"))


def _write_bp(outfile, bp, header, output_format):
    if output_format == "sufam":
        outfile.write("\t".join(bp.get(header, None).astype(str)) + "\n")
    elif output_format == "matrix":
        outfile.write("1" if bp.val_al_count > 0 else "0")
    else:
        raise(Exception("Unrecognized output format"))


def validate_mutations(vcffile, bam, reffa, sample, output_format, outfile):
    """Check if mutations in vcf are in bam"""
    header = []
    output_header = "sample chrom pos ref cov A C G T * - + " \
        "val_ref val_alt val_al_type val_al_count val_maf "\
        "most_common_indel most_common_indel_count most_common_indel_maf most_common_indel_type most_common_al " \
        "most_common_al_count most_common_al_maf most_common_count most_common_maf".split()

    if output_format == "sufam":
        outfile.write("\t".join(output_header))
        outfile.write("\n")
        for line in open(vcffile):
            if line.startswith("#CHROM"):
                header = line[1:].rstrip('\n').split("\t")
            if line.startswith("#"):
                continue
        if len(header) == 0:
            raise(Exception("No header found in vcf file, #CHROM not found"))
        record = dict(zip(header, line.rstrip('\n').split("\t")))
        record_type = "snv"
        if len(record["REF"]) > len(record["ALT"]):
            record_type = "deletion"
        elif len(record["ALT"]) > len(record["REF"]):
            record_type = "insertion"
            no_cov = pd.Series({
                "sample": sample,
                "chrom": record["CHROM"], "pos": record["POS"],
                "ref": record["REF"],
                "cov": 0, "A": 0, "C": 0, "G": 0, "T": 0,
                "val_ref": record["REF"], "val_alt": record["ALT"],
                "val_al_type": record_type, "val_al_count": 0, "val_maf": 0})
            bpdf = get_baseparser_extended_df(bam, sample, record["CHROM"], record["POS"], record["POS"], reffa,
                                              record["REF"], record["ALT"])
            if bpdf is None:
                bp = no_cov
            else:
                bp = bpdf.ix[0, :]
                _write_bp(outfile, bp, output_header, output_format)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reffa", type=str, help="Reference genome (fasta)")
    parser.add_argument("vcf", type=str, help="VCF with mutations to be validated")
    parser.add_argument("bam", type=str, help="BAM to find mutations in")
    parser.add_argument("--sample_name", type=str, default=None, help="Set name "
                        "of sample, used in output.")
    parser.add_argument("--format", type=str, choices=["matrix", "sufam"], default="sufam", help="Set output format")
    args = parser.parse_args()
    if args.sample_name is None:
        args.sample_name = args.bam
        validate_mutations(args.vcf, args.bam, args.reffa, args.sample_name, args.format, sys.stdout)


if __name__ == "__main__":
    main()
