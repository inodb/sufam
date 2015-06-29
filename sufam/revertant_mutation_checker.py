"""
Check if a mutation is reverted
"""
import argparse
import sys
from sufam import oncotator
import logging

import pyhgvs as hgvs
from Bio import SeqIO


def alter_coords_hgvs_sequential(h1, h2):
    """Change HGVS coords of h2 after applying h1"""
    if h1.kind == "c" and h2.kind == "c":
        if h1.mutation_type == ">":
            h3 = h2
        elif h1.mutation_type == "del":
            if h1.cdna_start.coord > h2.cdna_end.coord:
                h3 = h2
            elif h1.cdna_end.coord < h2.cdna_start.coord:
                h3 = h2
                h3.cdna_start = hgvs.CDNACoord(coord=h3.cdna_start.coord-len(h1.ref_allele))
                h3.cdna_end = hgvs.CDNACoord(coord=h3.cdna_end.coord-len(h1.ref_allele))
            else:
                raise(Exception("Overlapping del not implemented"))
        elif h1.mutation_type == "ins":
            if h1.cdna_start.coord > h2.cdna_end.coord:
                h3 = h2
            elif h1.cdna_end.coord < h2.cdna_start.coord:
                h3 = h2
                h3.cdna_start = hgvs.CDNACoord(coord=h3.cdna_start.coord+len(h1.alt_allele))
                h3.cdna_end = hgvs.CDNACoord(coord=h3.cdna_end.coord+len(h1.alt_allele))
            else:
                raise(Exception("Overlapping ins not implemented"))
        elif h1.mutation_type == "dup":
            if h1.cdna_start.coord > h2.cdna_end.coord:
                h3 = h2
            elif h1.cdna_end.coord < h2.cdna_start.coord:
                h3 = h2
                h3.cdna_start = hgvs.CDNACoord(coord=h3.cdna_start.coord+len(h1.alt_allele))
                h3.cdna_end = hgvs.CDNACoord(coord=h3.cdna_end.coord+len(h1.alt_allele))
            else:
                raise(Exception("Overlapping dup not implemented"))
        elif h1.mutation_type == "delins":
            if h1.cdna_start.coord > h2.cdna_end.coord:
                h3 = h2
            elif h1.cdna_end.coord < h2.cdna_start.coord:
                h3 = h2
                h3.cdna_start = hgvs.CDNACoord(coord=h3.cdna_start.coord-len(h1.ref_allele)+len(h1.alt_allele))
                h3.cdna_end = hgvs.CDNACoord(coord=h3.cdna_end.coord-len(h1.ref_allele)+len(h1.alt_allele))
            else:
                raise(Exception("Overlapping delins not implemented"))
        else:
            raise(Exception("Unexpected mutation_type {}".format(h1.mutation_type)))
        return h3
    else:
        raise(Exception("Only cDNA mutations have been implemented"))


def apply_hgvs(seq, h):
    """Apply 1-based HGVS mutation to 0-based biopython Sequence"""
    if h.kind == "c":
        start = h.cdna_start.coord - 1
        end = h.cdna_end.coord - 1
        if h.mutation_type == ">":
            assert(seq[start] == h.ref_allele)
            new_seq = seq[:start] + h.alt_allele + seq[start+1:]
        elif h.mutation_type == "del":
            assert(seq[start:end+1] == h.ref_allele)
            new_seq = seq[:start] + seq[end+1:]
            assert(len(seq) == len(new_seq) + len(h.ref_allele))
        elif h.mutation_type == "ins":
            new_seq = seq[:start+1] + h.alt_allele + seq[end:]
            assert(len(seq) + len(h.alt_allele) == len(new_seq))
        elif h.mutation_type == "dup":
            assert(seq[start:end+1] == h.ref_allele)
            new_seq = seq[:end+1] + h.alt_allele + seq[end+1:]
        elif h.mutation_type == "delins":
            assert(seq[start:end+1] == h.ref_allele)
            new_seq = seq[:start] + h.alt_allele + seq[end+1:]
            assert(len(seq) - len(h.ref_allele) + len(h.alt_allele) == len(new_seq))
        else:
            raise(Exception("Unexpected mutation_type {}".format(h.mutation_type)))
        return new_seq
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
    transcripts = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))

    for m in [l.rstrip('\n') for l in open(vcf).readlines()]:
        hgvs_mut = hgvs.HGVSName(m)
        for hgvs_rev_mut in ot.get_hgvs_mutations(hgvs_mut.transcript):
            if is_revertant(transcripts[hgvs_mut.transcript], hgvs_mut, alter_coords_hgvs_sequential(hgvs_mut, hgvs_rev_mut)):
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
    parser.add_argument("oncotator_file", type=str, help="MAF to find mutations in")
    parser.add_argument("fasta", type=str, help="Fasta file with transcripts")
    args = parser.parse_args()
    check_revertant_mutations(args.vcf, args.oncotator_file, args.fasta)


if __name__ == "__main__":
    main()
