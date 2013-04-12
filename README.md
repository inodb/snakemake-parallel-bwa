This is a Snakefile to do alignments using bwa. The file is split up into a
user defined number of reads. This way the alignment on each chunk of reads can
be scheduled separately. Currently only supports sbatch, but it should be
straightforward to add your own.
