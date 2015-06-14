"""
Check if a mutation is reverted
"""
import argparse
from sufam import mutation
from sufam import oncotator
import pyhgvs as hgvs


def is_revertant_mutation(record, hgvs_mut, hgvs_rev_mut):


def check_revertant_mutations(vcf, oncotator_file, fasta):
    ot = oncotator.Oncotator(oncotator_file)
    trascripts = \
        SeqIO.to_dict(SeqIO.parse("tests/test_data/BRCA1_transcripts.fa",
        "fasta"))

    for m in [l.rstrip('\n') for l in open(vcf).readlines()]:
        hgvs_mut = hgvs.HGVSName(m)
        for hgvs_rev_mut in ot.get_hgvs_mutations(h.transcript):
            is_revertant(transcripts[h.transcript], hgvs_mut, hgvs_rev_mut)

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", type=str, help="VCF with mutation to be reverted")
    parser.add_argument("oncotator_file", type=str, help="BAM to find mutations in")
    parser.add_argument("fasta", type=str, help="Fasta file with transcripts")
    args = parser.parse_args()
    check_revertant_mutations(args.vcf, args.oncotator_file, args.fasta)


if __name__ == "__main__:
    main()
