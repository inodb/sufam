MUTATION_TYPES = ["-", "+", "."]

class Mutation(object):
    def __init__(self, chrom, pos, type, change):
        if type in MUTATION_TYPES:
            self.type = type
        else:
            raise(Exception("Unexpected mutation type should be one of "
                            "{}".format(MUTATION_TYPES)))
        self.chrom = chrom
        self.pos = pos
        self.change = change

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


def get_mutation(chrom, pos, ref, alt):
    if len(ref) == len(alt):
        return Mutation(chrom, pos, ".", alt)
    elif len(ref) > len(alt):
        assert(len(alt) == 1)
        assert(ref[0] == alt)
        return Mutation(chrom, pos, "-", ref[1:])
    elif len(ref) < len(alt):
        assert(len(ref) == 1)
        assert(alt[0] == ref)
        return Mutation(chrom, pos, "+", alt[1:])


def parse_vcf(vcf_file):
    header = []
    mutations = []
    for line in open(vcf_file):
        if line.startswith("#CHROM"):
            header = line[1:].rstrip('\n').split("\t")
        if line.startswith("#"):
            continue
        if len(header) == 0:
            raise(Exception("No header found in vcf file, #CHROM not found"))
        record = dict(zip(header, line.rstrip('\n').split("\t")))
        mutations += [get_mutation(record["CHROM"], record["POS"], record["REF"], record["ALT"])]

    return mutations
