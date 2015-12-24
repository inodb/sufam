#! /usr/bin/env python
"""
Author: inodb and limr
Loosely based on niknafs original
"""
import sys
import subprocess
import warnings
import re

from sufam.mutation import Mutation, MutationsAtSinglePosition
from collections import Counter

MPILEUP_DEFAULT_PARAMS = '--ignore-RG --min-MQ 1 --max-depth 250000 --max-idepth 250000'


class ParseString(object):

    def __init__(self, ref, string):
        self.ref = ref.upper()
        self.string = string.upper()
        self.types = {'A': 0, 'G': 0, 'C': 0, 'T': 0, '-': [], '*': 0, '+': []}
        self.process()

    def process(self):
        # rm two characters when encountering '^' as it indicates
        # a read start mark and the read mapping quality
        self.string = re.sub(re.escape("^") + ".", '', self.string)
        # process and remove all indels from the pileup
        while 1:
            m = re.search(r"([+-])(\d+)([ATCG]+)", self.string)
            if m is None:
                break
            typ = m.group(1)
            le = int(m.group(2))
            # get indel of given length
            seq = m.group(3)[0:le]
            self.types[typ].append(seq)
            # remove only indel (eg not C after +2ATC)
            self.oldstring = self.string
            self.string = self.string[:m.start()] + \
                self.string[m.start() + 1 + len(m.group(2)) + le:]
        self.types[self.ref] = len(re.findall(r"[.,]", self.string))
        for x in re.findall(r"[ATCG*]", self.string):
            self.types[x] += 1
        return

    def get_mutations(self, chrom, pos, cov, ref):
        pos = int(pos)
        muts = MutationsAtSinglePosition(chrom, pos, cov, ref[0])

        for base in {"A", "C", "G", "T"} - set([self.ref[0]]):
            m = Mutation(chrom, pos, '.', base)
            m.count = self.types[base]
            m.cov = cov
            m.ref = self.ref[0]
            muts.add_snv(m)

        insertions = Counter(self.types["+"])
        for ins in insertions:
            m = Mutation(chrom, pos, '+', ins)
            m.count = insertions[ins]
            m.cov = cov
            m.ref = self.ref[0]
            muts.add_insertion(m)

        deletions = Counter(self.types["-"])
        for deletion in deletions:
            m = Mutation(chrom, pos, '-', deletion)
            m.count = deletions[deletion]
            m.cov = cov
            m.ref = self.ref[0]
            muts.add_deletion(m)

        return muts

    def __repr__(self):
        types = self.types
        return '\t'.join(map(str, [types['A'], types['C'], types['G'], types['T'],
                         types['*']]) + map(','.join, [types['-'], types['+']]))


def get_mutations(line):
    toks = line.strip('\n').split('\t')
    ref = toks[2].upper()
    cov = toks[3]
    return ParseString(ref, toks[4]).get_mutations(toks[0], int(toks[1]), cov, ref)


def parse(line):
    toks = line.strip('\n').split('\t')
    ref = toks[2].upper()
    cov = toks[3]
    return '\t'.join([toks[0], toks[1], ref, cov]) + '\t' + str(ParseString(ref, toks[4]))


def run(bam, chrom, pos1, pos2, reffa, parameters):
    """Run mpileup on given chrom and pos"""
    posmin = min(pos1, pos2)
    posmax = max(pos1, pos2)
    cmd = "samtools view -bh {bam} {chrom}:{pos1}-{pos2} " \
        "| samtools mpileup {parameters} -f {reffa} -".format(bam=bam, chrom=chrom,
                                                              pos1=posmin, pos2=posmax,
                                                              reffa=reffa, parameters=parameters)
    if pos1 == pos2:
        cmd += " | awk '$2 == {pos}'".format(pos=pos1)
    else:
        cmd += " | tail -n +2 | awk '$2 >= {posmin} && $2 <= {posmax}'".format(posmin=posmin, posmax=posmax)
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
        return stdout


def run_and_parse(bam, chrom, pos1, pos2, reffa, mpileup_parameters=MPILEUP_DEFAULT_PARAMS):
    return [parse(line) for line in run(bam, chrom, pos1, pos2, reffa, mpileup_parameters).split("\n")[:-1]]


def run_and_get_mutations(bam, chrom, pos1, pos2, reffa, mpileup_parameters=MPILEUP_DEFAULT_PARAMS):
    return [get_mutations(line) for line in run(bam, chrom, pos1, pos2, reffa, mpileup_parameters).split("\n")[:-1]]


def main():
    print >>sys.stdout, "chrom\tpos\tref\tcov\tA\tC\tG\tT\t*\t-\t+"
    for line in sys.stdin:
        print >>sys.stdout, parse(line)

if __name__ == '__main__':
    main()
