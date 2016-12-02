So U Found A Mutation? (SUFAM)
==============================
.. image:: https://zenodo.org/badge/14279/inodb/sufam.svg
   :target: https://zenodo.org/badge/latestdoi/14279/inodb/sufam
.. image:: https://badge.fury.io/py/sufam.svg
    :target: http://badge.fury.io/py/sufam
.. image:: http://anaconda.org/inodb/sufam/badges/version.svg
    :target: http://anaconda.org/inodb/sufam
.. image:: https://travis-ci.org/inodb/sufam.svg?branch=master
    :target: https://travis-ci.org/inodb/sufam

Found a mutation in one or more samples? Now you want to check if they are in
another sample. Unfortunately mutect, varscan or whatever other variant caller
is not calling them. Use SUFAM. The super sensitive validation caller that
calls everything on a given position. All you need is a vcf with the mutations
that you are interested in and the sam/bam file of the sample where you want to
find the same inconsipicuous mutation.

Installation
------------
::

    pip install sufam

Run
---
::

	usage: sufam [-h] [--sample_name SAMPLE_NAME [SAMPLE_NAME ...]]
				[--format {matrix,sufam,vcf}]
				[--mpileup-parameters MPILEUP_PARAMETERS] [--version]
				reffa vcf bam [bam ...]

	So U Found A Mutation? (SUFAM)

	Found a mutation in one or more samples? Now you want to check if they are in
	another sample. Unfortunately mutect, varscan or whatever other variant caller
	is not calling them. Use SUFAM. The super sensitive validation caller that
	calls everything on a given position. All you need is a vcf with the mutations
	that you are interested in and the sam/bam file of the sample where you want to
	find the same inconsipicuous mutation.

	Author: inodb

	positional arguments:
	reffa                 Reference genome (fasta)
	vcf                   VCF with mutations to be validated
	bam                   BAMs to find mutations in (only --format vcf supports
							> 1)

	optional arguments:
	-h, --help            show this help message and exit
	--sample_name SAMPLE_NAME [SAMPLE_NAME ...]
							Set name of sample, used in output [name of bam].
	--format {matrix,sufam,vcf}
							Set output format [sufam]
	--mpileup-parameters MPILEUP_PARAMETERS
							Set options for mpileup [--ignore-RG --min-MQ 1 --max-
							depth 250000 --max-idepth 250000]
	--version             show program's version number and exit



Example
~~~~~~~
VCF file like::

	#CHROM POS     ID      REF     ALT     QUAL    FILTER  INFO
	17     7574012 COSM11286,COSM214290    C       G       100     PASS    .
	17     7574012 COSM11286,COSM214290    C       A       100     PASS    .


Check if given mutations are in a bam file::

    sufam human_g1k_v37_chr17.fa mutations.vcf subset1.bam 2> example/sufam.log > example/sufam.tsv

Output:

- `example/sufam.log <example/sufam.log>`_
- `example/sufam.tsv <example/sufam.tsv>`_

Or if you want vcf output where sufam output is added as GT/AD/DP format fields in input vcf::

    sufam --format vcf human_g1k_v37_chr17.fa mutations.vcf subset1.bam > example/sufam.vcf

Output:

- `example/sufam.vcf <example/sufam.vcf>`_

Or for multiple bams in a single vcf::

	sufam --sample_name subset1 subset2 subset3 \
	      --format vcf \
		  human_g1k_v37_chr17.fa \
		  mutations.vcf \
		  subset1.bam subset2.bam \
		  subset3.bam \
	> example/sufam_multibam.vcf

Output:

- `example/sufam_multibam.vcf <example/sufam_multibam.vcf>`_
 
Developers
----------
Tests
~~~~~
In root dir run::

    nosetests

For individual tests::

    nosetests -s tests/test_validation.py:TestValidation.test_validate_mutations_indel
