#! /usr/bin/env python
"""
Author: inodb and limr
Loosely based on niknafs original
"""
import sys
import re


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

    def __repr__(self):
        types = self.types
        return '\t'.join(map(str, [types['A'], types['C'], types['G'], types['T'],
                         types['*']]) + map(','.join, [types['-'], types['+']]))


def parse(line):
    toks = line.strip('\n').split('\t')
    ref = toks[2].upper()
    cov = toks[3]
    return '\t'.join([toks[0], toks[1], ref, cov]) + '\t' + str(ParseString(ref, toks[4]))


def main():
    print >>sys.stdout, "chrom\tpos\tref\tcov\tA\tC\tG\tT\t*\t-\t+"
    for line in sys.stdin:
        print >>sys.stdout, parse(line)

if __name__ == '__main__':
    main()
