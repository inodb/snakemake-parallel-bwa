The purpose of this code is to demonstrate using snakemake
(https://bitbucket.org/johanneskoester/snakemake/overview) to make pipelines
and schedule each rule as a separate job with proper dependencies and sbatch.
See also this issue:

https://bitbucket.org/johanneskoester/snakemake/issue/28/clustering-jobs-with-snakemake

The Snakefile aligns reads using bwa. The file is split up into a user defined
number of reads. This way the alignment of each chunk of reads can be scheduled
separately.

# Requirements
* Snakemake (https://bitbucket.org/johanneskoester/snakemake/overview)
* bwa (http://bio-bwa.sourceforge.net/bwa.shtml)
* samtools (http://samtools.sourceforge.net/)
* genomeCoverageBed (http://code.google.com/p/bedtools/)
* Picard (http://sourceforge.net/apps/mediawiki/picard/index.php?title=Main_Page)

# Install
* Change the path to the location of ``sbatch_job`` in  ``Snakefile-sbatch.py``
* Change the location of Picard tools in ``Snakefile``.

# Test
```
# make sure you are in the root dir of this repo
cd test
snakemake -j 99 --debug --immediate-submit --cluster './../Snakefile-sbatch.py {dependencies}' all
```
