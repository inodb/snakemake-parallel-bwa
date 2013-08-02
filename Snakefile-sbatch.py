#!/usr/bin/env python
import sys
import subprocess
import os
import math
import ipdb
import errno


def make_dir(directory):
    """Make directory unless existing. Ignore error in the latter case."""
    try:
        os.makedirs(directory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


class SnakeJob:
    """Snakemake can generate bash scripts that can be sumbitted by a
    scheduler.  This class reads the bash script and stores the number of the
    rule, name of bash file and the supplied input files."""
    def __init__(self, snakebashfile, dependencies=None):
        fh = open(snakebashfile, 'r')
        self.scriptname = snakebashfile
        fh.readline()
        self.rule = fh.readline().split()[1]
        self.ifiles = fh.readline().split()[1:]
        self.ofiles = fh.readline().split()[1:]
        if dependencies == None:
            self.dependencies = None
        else:
            # expects snakemake like string
            self.dependencies = dependencies.split(',')
            assert len(self.dependencies) >= 1
        fh.close()


class UndefinedJobRule(Exception):
    """Exception in case an sbatch job has no defined resource usage in the
    code."""
    def __init__(self, msg):
        self.msg = msg


class SnakeJobSbatch(SnakeJob):
    def __init__(self, snakebashfile, dependencies=None):
        SnakeJob.__init__(self, snakebashfile, dependencies)
        if self.dependencies == None:
            self.dep_str = ''
        else:
            self.dep_str = '-d ' + ','.join(["afterok:%s" % d for d in self.dependencies])

    def schedule(self):
        """Schedules a snakemake job with sbatch and determines resource usage
        based on input files."""
        # create the output directory, so slurm output can go there
        #make_dir(os.path.dirname(os.path.abspath(self.ofiles[0])))

        if self.rule == 'split1' or self.rule == 'split2' or self.rule == 'bwt':
            #hours = math.ceil(os.path.getsize(self.ifiles[0]) / float(500 * 1024 ** 2))
            minutes = 15
            sbatch_cmd = 'sbatch %s -A b2010008 -p core --qos=short -t 00:%i:00 --output=snakemake-%%j.out -J %s ~/bin/sbatch_job bash %s' % (self.dep_str, minutes, self.rule, self.scriptname)
        elif self.rule in ['aln', 'sort']:
            minutes = 15
            sbatch_cmd = 'sbatch %s -A b2010008 -p node --qos=short -t 00:%i:00 --output=snakemake-%%j.out -J %s ~/bin/sbatch_job bash %s' % (self.dep_str, minutes, self.rule, self.scriptname)
        elif self.rule == 'sampe':
            minutes = 15
            sbatch_cmd = 'sbatch %s -A b2010008 -p core --qos=short -t 00:%i:00 --output=snakemake-%%j.out -J %s ~/bin/sbatch_job bash %s' % (self.dep_str, minutes, self.rule, self.scriptname)
        elif self.rule == 'samtobam':
            minutes = 15
            sbatch_cmd = 'sbatch %s -A b2010008 -p core --qos=short -t 00:%i:00 --output=snakemake-%%j.out -J %s ~/bin/sbatch_job bash %s' % (self.dep_str, minutes, self.rule, self.scriptname)
        elif self.rule == 'merge':
            minutes = 15
            sbatch_cmd = 'sbatch %s -A b2010008 -p core --qos=short -t 00:%i:00 --output=snakemake-%%j.out -J %s ~/bin/sbatch_job bash %s' % (self.dep_str, minutes, self.rule, self.scriptname)
        elif self.rule in ['removeduplicates', 'cleansam', 'index', 'coverage', 'mean_coverage_per_contig']:
            minutes = 15
            sbatch_cmd = 'sbatch %s -A b2010008 -p core --qos=short -t 00:%i:00 --output=snakemake-%%j.out -J %s ~/bin/sbatch_job bash %s' % (self.dep_str, minutes, self.rule, self.scriptname)
        elif self.rule in ['split_all', 'clean', 'all', 'merge_all', 'test']:
            # no scheduling just run
            sbatch_cmd = 'bash %s' % (self.scriptname)
        else:
            raise UndefinedJobRule('Undefined resource usage %s' % (self.rule))
            return 2

        print(sbatch_cmd)
        popenrv = subprocess.Popen(sbatch_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).communicate()
        print("%i" % int(popenrv[0].split()[-1]))


if __name__ == '__main__':

    if len(sys.argv) == 2:
        sj = SnakeJobSbatch(sys.argv[1])
    elif len(sys.argv) == 3:
        sj = SnakeJobSbatch(sys.argv[1], sys.argv[2])
    else:
        raise Exception('Expected snakemake arguments, should max be 3')
    try:
        sj.schedule()
    except UndefinedJobRule as err:
        print(err.msg, file=sys.stderr)
        sys.exit(2)
