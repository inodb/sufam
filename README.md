# So U Found A Mutation? (SUFAM)

Found a mutation in one or more samples? Now you want to check if they are in
another sample. Unfortunately mutect, varscan or whatever other variant caller
is not calling them. Use SUFAM. The super sensitive validation caller that
calls everything on a given position. All you need is a vcf with the mutations
that you are interested in and the sam/bam file of the sample where you want to
find the same inconsipicuous mutation.

## Installation

```
git clone http://github.com/inodb/sufam
cd sufam
python setup.py install
```

## Run

<!--
(echo '```'; echo '$ sufam --help'; sufam --help | head -1; sufam --help | awk 'BEGIN {flip=0} {if (!flip) { if ($0 ~ "Author") {flip=1}} else {print $0}}'; echo '``'; ) >> README.md
-->

```
$ sufam --help
usage: sufam [-h] [--sample_name SAMPLE_NAME] reffa vcf bam

positional arguments:
  reffa                 Reference genome (fasta)
  vcf                   VCF with mutations to be validated
  bam                   BAM to find mutations in

optional arguments:
  -h, --help            show this help message and exit
  --sample_name SAMPLE_NAME
                        Set name of sample, used in output.
```
