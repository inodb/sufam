"""
Check if a mutation is reverted
"""
import argparse
import sys
from sufam import oncotator
import logging

import pyhgvs as hgvs
from Bio import SeqIO


def apply_hgvs(seq, h):
    """Apply 1-based HGVS mutation to 0-based biopython Sequence"""
    if h.kind == "c":
        start = h.cdna_start.coord - 1
        end = h.cdna_end.coord - 1
        if h.mutation_type == ">":
            assert(seq[start] == h.ref_allele)
            seq = seq[start] = h.alt_allele
        elif h.mutation_type == "del":
            assert(seq[start:end+1] == h.ref_allele)
            seq = seq[:start] + seq[end:]
        elif h.mutation_type == "ins":
            assert(seq[start] == h.ref_allele)
            # TODO: is this correct?
            seq = seq[:start] + h.alt_allele + seq[end:]
        elif h.mutation_type == "dup":
            assert(seq[start:end+1] == h.ref_allele)
            seq = seq[:start] + h.alt_allele + seq[end:]
        else:
            raise(Exception("Unexpected mutation_type {}".format(h.mutation_type)))
        return seq
    else:
        raise(Exception("Only cDNA mutations have been implemented"))


def is_revertant(record, hgvs_mut, hgvs_rev_mut):
    normal_p = record.seq.translate(to_stop=True)
    mut_p = apply_hgvs(record.seq, hgvs_mut).translate(to_stop=True)
    revmut_p = apply_hgvs(apply_hgvs(record.seq, hgvs_mut), hgvs_rev_mut).translate(to_stop=True)
    logging.info("---")
    logging.info("mutation: {}".format(hgvs_mut))
    logging.info("putative revertant mutation: {}".format(hgvs_rev_mut))
    logging.info("transcript: {}".format(record.id))
    logging.info("normal length: {}".format(len(normal_p)))
    logging.info("mutation length: {}".format(len(mut_p)))
    logging.info("putative revertant length: {}".format(len(revmut_p)))
    logging.info("---")
    return len(normal_p) == len(revmut_p)


def check_revertant_mutations(vcf, oncotator_file, fasta):
    ot = oncotator.Oncotator(oncotator_file)
    transcripts = \
        SeqIO.to_dict(SeqIO.parse("tests/test_data/BRCA1_transcripts.fa",
        "fasta"))

    for m in [l.rstrip('\n') for l in open(vcf).readlines()]:
        hgvs_mut = hgvs.HGVSName(m)
        for hgvs_rev_mut in ot.get_hgvs_mutations(hgvs_mut.transcript):
            if is_revertant(transcripts[hgvs_mut.transcript], hgvs_mut, hgvs_rev_mut):
                logging.info("REVERTANT MUTATION FOUND ZOMG")


def main():
    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO,
        # format='%(asctime)s:%(levelname)s:%(name)s:%(message)s',
    )
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", type=str, help="VCF with mutation to be reverted")
    parser.add_argument("oncotator_file", type=str, help="MAF from Oncotator to find revertant mutations in")
    parser.add_argument("fasta", type=str, help="Fasta file with transcripts")
    args = parser.parse_args()
    check_revertant_mutations(args.vcf, args.oncotator_file, args.fasta)


if __name__ == "__main__":
    main()
