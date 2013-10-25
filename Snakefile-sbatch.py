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
        if dependencies == None or len(dependencies) < 1:
            self.dependencies = None
        else:
            # expects snakemake like list of numbers
            self.dependencies = dependencies
            assert len(self.dependencies) >= 1
        fh.close()


class UndefinedJobRule(Exception):
    """Exception in case an sbatch job has no defined resource usage in the
    code."""
    def __init__(self, msg):
        self.msg = msg


class SnakeJobSbatch(SnakeJob):
    # Change this to the path of the sbatch_job wrapper script
    sbatch_job_path = "~/bin/sbatch_job"

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

        run_locally = False

        if self.rule == 'split1' or self.rule == 'split2' or self.rule == 'bwt':
            #hours = math.ceil(os.path.getsize(self.ifiles[0]) / float(500 * 1024 ** 2))
            minutes = 10
            sbatch_cmd = 'sbatch %s -A b2010008 -p core -t 00:%i:00 --output=snakemake-%%j.out -J %s %s bash %s' % (self.dep_str, minutes, self.rule, self.sbatch_job_path, self.scriptname)
        elif self.rule in ['aln', 'sort']:
            minutes = 10
            sbatch_cmd = 'sbatch %s -A b2010008 -p node -t 00:%i:00 --output=snakemake-%%j.out -J %s %s bash %s' % (self.dep_str, minutes, self.rule, self.sbatch_job_path, self.scriptname)
        elif self.rule == 'sampe':
            minutes = 10
            sbatch_cmd = 'sbatch %s -A b2010008 -p core -t 00:%i:00 --output=snakemake-%%j.out -J %s %s bash %s' % (self.dep_str, minutes, self.rule, self.sbatch_job_path, self.scriptname)
        elif self.rule == 'samtobam':
            minutes = 10
            sbatch_cmd = 'sbatch %s -A b2010008 -p core -t 00:%i:00 --output=snakemake-%%j.out -J %s %s bash %s' % (self.dep_str, minutes, self.rule, self.sbatch_job_path, self.scriptname)
        elif self.rule == 'merge':
            minutes = 10
            sbatch_cmd = 'sbatch %s -A b2010008 -p core -t 00:%i:00 --output=snakemake-%%j.out -J %s %s bash %s' % (self.dep_str, minutes, self.rule, self.sbatch_job_path, self.scriptname)
        elif self.rule in ['removeduplicates', 'cleansam', 'index', 'coverage', 'mean_coverage_per_contig']:
            minutes = 10
            sbatch_cmd = 'sbatch %s -A b2010008 -p core -t 00:%i:00 --output=snakemake-%%j.out -J %s %s bash %s' % (self.dep_str, minutes, self.rule, self.sbatch_job_path, self.scriptname)
        elif self.rule in ['split_all', 'clean', 'all', 'merge_all', 'test']:
            minutes = 1
            sbatch_cmd = 'sbatch %s -A b2010008 -p core -t 00:%i:00 --output=snakemake-%%j.out -J %s %s bash %s' % (self.dep_str, minutes, self.rule, self.sbatch_job_path, self.scriptname)
            #sbatch_cmd = 'bash %s' % (self.scriptname)
            #run_locally = True
        else:
            raise UndefinedJobRule('Undefined resource usage %s' % (self.rule))
            return 2

        #print(sbatch_cmd)
        popenrv = subprocess.Popen(sbatch_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True).communicate()
        #print(popenrv[0])
        if not run_locally:
            try:
                print("%i" % int(popenrv[0].split()[-1]))
            except ValueError:
                print("Not a submitted job: %s" % popenrv[0])
                sys.exit(2)


if __name__ == '__main__':
    sj = SnakeJobSbatch(sys.argv[-1], sys.argv[1:-1])
    try:
        sj.schedule()
    except UndefinedJobRule as err:
        print(err.msg, file=sys.stderr)
        sys.exit(2)
