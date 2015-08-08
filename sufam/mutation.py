MUTATION_TYPES = ["-", "+", "."]


class MutationsAtSinglePosition(object):
    def __init__(self, chrom, pos, cov, ref=None):
        self.snvs = {}
        self.deletions = {}
        self.insertions = {}
        self.chrom = chrom
        self.pos = pos
        self.cov = cov
        if ref:
            self.ref = ref

    def _sanity_check(self, mutation):
        return self._check_same_loc(mutation) and \
            mutation.count and mutation.cov

    def _check_same_loc(self, mutation):
        return self.pos == mutation.pos and self.chrom == mutation.chrom

    def add_snv(self, mutation):
        self.snvs[mutation.change] = mutation

    def add_insertion(self, mutation):
        assert(self._sanity_check)
        self.insertions[mutation.change] = mutation

    def add_deletion(self, mutation):
        assert(self._sanity_check)
        self.deletions[mutation.change] = mutation

    def add_mutation(self, mutation):
        if mutation.type == ".":
            self.add_snv(mutation)
        elif mutation.type == "+":
            self.add_insertion(mutation)
        elif mutation.type == "-":
            self.add_deletion(mutation)
        else:
            raise(Exception("Unknown mutation type: {}".format(mutation)))

    @staticmethod
    def from_mutation_list(l):
        muts = MutationsAtSinglePosition(l[0].chrom, l[0].pos, l[0].cov)

        for m in l:
            muts.add_mutation(m)
        return muts

    def filter_against_normal(self, normal_mutations, maf_min=0.2,
                              maf_count_threshold=20, count_min=1):
        """Filters mutations that are in the given normal"""
        assert(normal_mutations.chrom == self.chrom)
        assert(normal_mutations.pos == self.pos)

        def passes_normal_criteria(mut):
            return (mut.count >= maf_count_threshold and mut.maf > maf_min) or \
                (mut.count < maf_count_threshold and mut.count > count_min)

        nms = normal_mutations
        muts = MutationsAtSinglePosition(self.chrom, self.pos, self.cov)

        for snv in self.snvs:
            if not (snv in nms.snvs and passes_normal_criteria(nms.snvs[snv])):
                muts.add_snv(self.snvs[snv])

        for dlt in self.deletions:
            if not (dlt in nms.deletions and passes_normal_criteria(nms.deletions[dlt])):
                muts.add_deletion(self.deletions[dlt])

        for ins in self.insertions:
            if not (ins in nms.insertions and passes_normal_criteria(nms.insertions[ins])):
                muts.add_insertion(self.insertions[ins])

        return muts

    def __len__(self):
        return len(self.snvs) + len(self.deletions) + len(self.insertions)

    def __iter__(self):
        for v in self.snvs.itervalues():
            yield v
        for v in self.insertions.itervalues():
            yield v
        for v in self.deletions.itervalues():
            yield v

    def __str__(self):
        return str(self.__dict__)


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

    @property
    def maf(self):
        """Requires one to set count and coverage"""
        return self.count / float(self.cov)

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        return self.chrom == other.chrom and self.pos == other.pos and \
            self.change == other.change and self.type == other.type


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
