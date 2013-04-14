#!/usr/bin/env python
import sys
import subprocess
import os
import math


class SnakeJob:
    """Snakemake can generate bash scripts that can be sumbitted by a
    scheduler.  This class reads the bash script and stores the number of the
    rule, name of bash file and the supplied input files."""
    def __init__(self, snakebashfile):
        fh = open(snakebashfile, 'r')
        self.scriptname = snakebashfile
        fh.readline()
        self.rule = fh.readline().split()[1]
        self.ifiles = fh.readline().split()[1:]
        self.ofiles = fh.readline().split()[1:]
        fh.close()


class UndefinedJobRule(Exception):
    """Exception in case an sbatch job has no defined resource usage in the
    code."""
    def __init__(self, msg):
        self.msg = msg


class SnakeJobSbatch(SnakeJob):
    def schedule(self):
        """Schedules a snakemake job with sbatch and determines resource usage
        based on input files."""
        if self.rule == 'split1' or self.rule == 'split2' or self.rule == 'bwt':
            #hours = math.ceil(os.path.getsize(self.ifiles[0]) / float(500 * 1024 ** 2))
            hours = 1
            sbatch_cmd = 'sbatch -A b2010008 -p core -t %i:00:00 --output=snakemake-%%j.out -J %s ~/bin/sbatch_job bash %s' % (hours, self.rule, self.scriptname)
        elif self.rule == 'aln':
            hours = 1
            sbatch_cmd = 'sbatch -A b2010008 -p core -t %i:00:00 --output=snakemake-%%j.out -J %s ~/bin/sbatch_job bash %s' % (hours, self.rule, self.scriptname)
        elif self.rule == 'sampe':
            hours = 1
            sbatch_cmd = 'sbatch -A b2010008 -p core -t %i:00:00 --output=snakemake-%%j.out -J %s ~/bin/sbatch_job bash %s' % (hours, self.rule, self.scriptname)
        elif self.rule == 'samtobam':
            hours = 1
            sbatch_cmd = 'sbatch -A b2010008 -p core -t %i:00:00 --output=snakemake-%%j.out -J %s ~/bin/sbatch_job bash %s' % (hours, self.rule, self.scriptname)
        elif self.rule == 'merge':
            hours = 1
            sbatch_cmd = 'sbatch -A b2010008 -p core -t %i:00:00 --output=snakemake-%%j.out -J %s ~/bin/sbatch_job bash %s' % (hours, self.rule, self.scriptname)
        elif self.rule in ['split_all', 'clean', 'all', 'merge_all', 'test']:
            # no scheduling just run
            sbatch_cmd = 'bash %s' % (self.scriptname)
        else:
            raise UndefinedJobRule('Undefined resource usage %s' % (self.rule))
            return 2

        print(sbatch_cmd)
        subprocess.Popen(sbatch_cmd, shell=True)


if __name__ == '__main__':
    sj = SnakeJobSbatch(sys.argv[1])
    try:
        sj.schedule()
    except UndefinedJobRule as err:
        print(err.msg, file=sys.stderr)
        sys.exit(2)
