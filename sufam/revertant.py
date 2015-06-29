"""
Functions for revertant mutations
"""
def filter_out_mutations_in_normal(tumordf, normaldf):
    df = tumordf.merge(normaldf, on=["chrom","pos"], suffixes=("_T","_N"))

    # filters
    common_al = (df.most_common_al_count_T == df.most_common_count_T) & (df.most_common_al_T == df.most_common_al_N)
    common_indel = (df.most_common_indel_count_T == df.most_common_count_T) & (df.most_common_indel_T == df.most_common_indel_N)
    normal_criteria = ((df.most_common_count_N >= 20) & df.most_common_maf_N > 0.2) | ((df.most_common_count_N < 20) & df.most_common_count_N > 1)
    df = df[~(common_al | common_indel) & normal_criteria]

    # restore column names of tumor
    for c in df.columns:
        if c.endswith("_N"):
            del df[c]
    df.columns = [c[:-2] if c.endswith("_T") else c for c in df.columns]

    return df


def select_only_revertant_mutations(bpdf, pos, snv=None, ins=None, dlt=None):
    """
    Selects only mutations that revert the given mutations in a single event.
    """
    if sum([bool(snv), bool(ins), bool(dlt)]) != 1:
        raise(Exception("Should be either snv, ins or del".format(snv)))

    if bool(snv):
        if snv not in ["A","C","G","T"]:
            raise(Exception("snv {} should be A, C, G or T".format(snv)))
        rv = bpdf[(bpdf.pos.astype(float) == pos) & (bpdf.most_common_al == snv) & (bpdf.most_common_al_count == bpdf.most_common_count)]
    elif bool(ins):
        rv = bpdf[((bpdf.most_common_indel.apply(lambda x: len(x) + len(ins) % 3 if x else None) == 0 ) & (bpdf.most_common_indel_type == "+") & (bpdf.most_common_count == bpdf.most_common_indel_count)) |
                ((bpdf.most_common_indel.apply(lambda x: len(ins) - len(x) % 3 if x else None) == 0 ) & (bpdf.most_common_indel_type == "-") & (bpdf.most_common_count == bpdf.most_common_indel_count))]
    elif bool(dlt):
        rv = bpdf[((bpdf.most_common_indel.apply(lambda x: len(x) - len(dlt) % 3 if x else None) == 0) & (bpdf.most_common_indel_type == "+") & (bpdf.most_common_count == bpdf.most_common_indel_count)) |
                ((bpdf.most_common_indel.apply(lambda x: -len(dlt) - len(x) % 3 if x else None) == 0 ) & (bpdf.most_common_indel_type == "-") & (bpdf.most_common_count == bpdf.most_common_indel_count))]
    else:
        # should never happen
        raise(Exception("No mutation given?"))

    # find mutations that delete the mutation partially or completely
    del_mut_df = bpdf[(bpdf.most_common_indel_type == "-") & (bpdf.pos.astype(float) <= pos) & (bpdf.pos.astype(float) + bpdf.most_common_indel.apply(lambda x: len(x) if x else None) >= pos)]
    return pd.concat([rv, del_mut_df], axis=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", type=str, help="VCF with mutations to find revertant mutations for")
    parser.add_argument("bam", type=str, help="BAM to find mutations in")


if __name__ == "__main__":
    main()
