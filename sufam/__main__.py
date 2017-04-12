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
from collections import Counter, namedtuple
import sys
import warnings
import six

import numpy as np
import pandas as pd
import vcf

import sufam
from sufam import mpileup_parser


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
            al_maf = al_count / float(x["cov"]) if float(x["cov"]) > 0 else None
        elif len(x["val_alt"]) > len(x["val_ref"]):  # insertion
            query = x["val_alt"][len(x["val_ref"]):]
            al_count = Counter(x["+"].split(","))[query] if x["+"] is not None else 0
            al_type = "insertion"
            al_maf = al_count / float(x["cov"])
        else:  # deletion
            query = x["val_ref"][len(x["val_alt"]):]
            al_count = Counter(x["-"].split(","))[query] if x["-"] is not None else 0
            al_type = "deletion"
            al_maf = al_count / float(x["cov"])
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


def get_baseparser_extended_df(sample, bp_lines, ref, alt):
    """Turn baseParser results into a dataframe"""
    columns = "chrom\tpos\tref\tcov\tA\tC\tG\tT\t*\t-\t+".split()
    if bp_lines is None:
        return None

    # change baseparser output to get most common maf per indel
    bpdf = pd.DataFrame([[sample] + l.rstrip('\n').split("\t") for l in bp_lines if len(l) > 0],
                        columns=["sample"] + columns, dtype=np.object)
    bpdf[bpdf == ""] = None

    # remove zero coverage rows
    bpdf = bpdf[bpdf["cov"].astype(int) > 0]

    if len(bpdf) == 0:
        return None

    if ref and alt:
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
        outfile.write("\t".join(bp.where(pd.notnull(bp), np.nan).get(header, None).astype(str)) + "\n")
    elif output_format == "matrix":
        outfile.write("1\n" if bp.val_al_count > 0 else "0\n")
    else:
        raise(Exception("Unrecognized output format"))


def _write_bp_vcf(outfile, bps, vcf_writer, record):
    def determine_genotype(bp):
        ref = int(int(bp['cov']) - bp.val_al_count > 0)
        alt = int(bp.val_al_count > 0)
        return '{}/{}'.format(ref, alt)

    _CallDataFormat = namedtuple('CallDataFormat', 'GT AD DP'.split())

    samp_fmt = vcf.parser.Reader
    calls = []
    for bp in bps:
        call = vcf.model._Call(None,
                            None,
                            _CallDataFormat(GT=determine_genotype(bp),
                                            AD=[int(bp['cov']) - bp.val_al_count, bp.val_al_count],
                                            DP=bp['cov']))
        calls += [call]
    record.samples = calls
    vcf_writer.write_record(record)


def validate_mutations(vcffile, bams, reffa, samples, output_format, outfile,
                       mpileup_parameters=mpileup_parser.MPILEUP_DEFAULT_PARAMS):
    """Check if mutations in vcf are in bam"""
    output_header = "sample chrom pos ref cov A C G T * - + " \
        "val_ref val_alt val_al_type val_al_count val_maf "\
        "most_common_indel most_common_indel_count most_common_indel_maf most_common_indel_type most_common_al " \
        "most_common_al_count most_common_al_maf most_common_count most_common_maf".split()

    # for backwards compatibility
    # if bam or samples is a string, convert to list instead
    if isinstance(samples, six.string_types):
        samples = [samples]
    if isinstance(bams, six.string_types):
        bams = [bams]

    if output_format == 'vcf':
        vcf_reader = vcf.Reader(open(vcffile))
        vcf_reader.samples = samples
        for f in 'GT AD DP'.split():
            vcf_reader.formats[f] = vcf.parser._Format(id=f,
                                                       num='R' if f == 'AD' else 1,
                                                       type='String' if f == 'GT' else 'Integer',
                                                       desc='')
        vcf_writer = vcf.Writer(outfile, vcf_reader)
    else:
        vcf_reader = open(vcffile)


    if output_format == "sufam":
        outfile.write("\t".join(output_header))
        outfile.write("\n")
    for record in vcf_reader:
        if output_format != 'vcf':
            line = record
            if line.startswith("#CHROM"):
                header = line[1:].rstrip('\n').split("\t")
                # create spoof pyvcf record if vcf_reader is not used
                _Record = namedtuple('Record', header)
            if line.startswith("#"):
                continue
            if len(header) == 0:
                raise(Exception("No header found in vcf file #CHROM not found"))
            # zip all column values, except alt (needs to be list in pyvcf)
            record_args = dict(zip(header, line.rstrip('\n').split("\t")))
            record_args['ALT'] = [record_args['ALT']]
            record = _Record(**record_args)

        # determine type of mutation
        record_type = "snv"
        if len(record.ALT) > 1:
            warnings.warn("Multiple ALT in one record is not implemented - using first")
        if len(record.REF) > len(record.ALT[0]):
            record_type = "deletion"
        elif len(record.ALT[0]) > len(record.REF):
            record_type = "insertion"

        # no coverage results
        no_cov = pd.Series({
            "chrom": str(record.CHROM), "pos": str(record.POS),
            "ref": str(record.REF),
            "cov": 0, "A": 0, "C": 0, "G": 0, "T": 0,
            "val_ref": str(record.REF), "val_alt": str(record.ALT[0]),
            "val_al_type": record_type, "val_al_count": 0, "val_maf": 0})

        # collect mpileup baseparser results per bam
        bps = []
        for i, bam in enumerate(bams):
            sample = samples[i]
            no_cov['sample'] = sample
            bp_lines = mpileup_parser.run_and_parse(bam, str(record.CHROM), str(record.POS), str(record.POS), reffa, mpileup_parameters)
            bpdf = get_baseparser_extended_df(sample, bp_lines, str(record.REF), str(record.ALT[0]))
            if bpdf is None:
                bp = no_cov
            else:
                bp = bpdf.ix[0, :]
            bps += [bp]

        # output call
        if output_format == "vcf":
            _write_bp_vcf(outfile, bps, vcf_writer, record)
        else:
            # only one bam file supported for outputs other than vcf
            _write_bp(outfile, bps[0], output_header, output_format)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reffa", type=str, help="Reference genome (fasta)")
    parser.add_argument("vcf", type=str, help="VCF with mutations to be validated")
    parser.add_argument("bam", type=str, nargs='+', help="BAMs to find mutations in (only --format vcf supports > 1)")
    parser.add_argument("--sample_name", type=str, nargs='+', default=None, help="Set name "
                        "of sample, used in output [name of bam].")
    parser.add_argument("--format", type=str, choices=["matrix", "sufam", "vcf"], default="sufam",
                        help="Set output format [sufam]")
    parser.add_argument("--mpileup-parameters", type=str,  default=mpileup_parser.MPILEUP_DEFAULT_PARAMS,
                        help="Set options for mpileup [{}]".format(mpileup_parser.MPILEUP_DEFAULT_PARAMS))
    parser.add_argument("--version", action='version', version=sufam.__version__)
    args = parser.parse_args()
    if args.sample_name is None:
        args.sample_name = args.bam
    if len(args.bam) > 1 and args.format != 'vcf':
        raise(Exception('Multiple bam files is only supported for --format vcf'))
    if len(args.sample_name) != len(args.bam):
        raise(Exception('# of --sample_name arguments should be equal to # of bams'))
    validate_mutations(args.vcf, args.bam, args.reffa, args.sample_name,
                       args.format, sys.stdout, mpileup_parameters=args.mpileup_parameters)


if __name__ == "__main__":
    main()
