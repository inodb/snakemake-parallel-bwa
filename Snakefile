# Split parameters
SPLIT_OUT="split-out"
LINES_PER_FILE=os.environ.get("LINES_PER_FILE", "400000000")

# aln parameters
MAP_IM_OUT="map-intermediate-out"

# merge parameters
MERGE_OUT="map"

# picard location
PICARD_DIR="/bubo/home/h16/inod/src/picard-tools-1.89"

# Given there's a pair1.fastq and pair2.fastq, split up all reads
rule all:
    message: """(1) Splitting up pair1.fastq and pair2.fastq for parallelization (use --cluster)
                (2) map them to ref.fa
                (3) sort bam and get coverage info
    """
    input: MERGE_OUT + "/pair_ref-s-md-s.bam", MERGE_OUT + "/pair_ref-s-md-s.bam.bai",
           MERGE_OUT + "/pair_ref-s-md-s.coverage"

rule merge_all:
    input: MERGE_OUT + "/pair_ref.bam"

rule split_all:
    input: dynamic(SPLIT_OUT + "/pair1.fastq.split.{splitid}"), dynamic(SPLIT_OUT + "/pair2.fastq.split.{splitid}")

rule split1:
    input: "pair1.fastq"
    params: lines_per_file=str(LINES_PER_FILE), out_dir=SPLIT_OUT
    output: dynamic(SPLIT_OUT + "/pair1.fastq.split.{splitid}")
    shell: """
    mkdir -p {params.out_dir}
    split -d -l {params.lines_per_file} {input[0]} {params.out_dir}/pair1.fastq.split.
    """

rule split2:
    input: "pair2.fastq"
    params: lines_per_file=str(LINES_PER_FILE), out_dir=SPLIT_OUT
    output: dynamic(SPLIT_OUT + "/pair2.fastq.split.{splitid}")
    shell: """
    mkdir -p {params.out_dir}
    split -d -l {params.lines_per_file} {input[0]} {params.out_dir}/pair2.fastq.split.
    """

rule bwt:
    input: "ref.fa"
    output: "ref.fa.bwt"
    log: "ref.fa.bwt.sm.log"
    shell: """
	bwa index {input}
    """

rule aln:
    input: pairn=SPLIT_OUT + "/pair{pair_nr}.fastq.split.{id}", refbwt="ref.fa.bwt", ref="ref.fa"
    threads: 8
    params: out_dir=MAP_IM_OUT
    output: MAP_IM_OUT + "/pair{pair_nr}.fastq.split.{id}.sai"
    log: MAP_IM_OUT + "/pair{pair_nr}.fastq.split.{id}.sai.sm.log"
    shell: """
    mkdir -p {params.out_dir}
	bwa aln -t {threads} {input.ref} -{wildcards.pair_nr} {input.pairn} > {output}
    """

rule sampe:
    input: sai1=MAP_IM_OUT + "/pair1.fastq.split.{id}.sai", sai2=MAP_IM_OUT + "/pair2.fastq.split.{id}.sai", pair1=SPLIT_OUT + "/pair1.fastq.split.{id}", pair2=SPLIT_OUT + "/pair2.fastq.split.{id}", refbwt="ref.fa.bwt", ref="ref.fa"
    params: out_dir=MAP_IM_OUT
    output: MAP_IM_OUT + "/pair.fastq.split.{id}.sam"
    log: MAP_IM_OUT + "/pair.fastq.split.{id}.sam.sm.log"
    shell: """
    mkdir -p {params.out_dir}
	bwa sampe {input.ref} {input.sai1} {input.sai2} {input.pair1} {input.pair2} > {output}
    """

rule samtobam:
    input: "{name}.sam"
    output: "{name}.bam"
    log: "{name}.bam.sm.log"
    shell: """
    samtools view -bS {input} > {output}
    """

rule merge:
    input: dynamic(MAP_IM_OUT + "/pair.fastq.split.{splitid}.bam")
    params: out_dir=MERGE_OUT
    output: MERGE_OUT + "/pair_ref.bam"
    log: MERGE_OUT + "/pair_ref.bam.sm.log"
    shell: """
    mkdir -p {params.out_dir}
    samtools merge {output} {input}
    """

rule sort:
    input: "{name}.bam"
    output: "{name}-s.bam"
    log: "{name}-s.bam.sm.log"
    params: srtsam_jar=PICARD_DIR + "/SortSam.jar"
    shell: """
    java -jar {params.srtsam_jar} \
        INPUT={input} \
        OUTPUT={output} \
        SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=SILENT
    """

# Removes reads mapping to the exact same location, probably duplicates
rule removeduplicates:
    input: "{name}-s.bam"
    output: "{name}-s-md.bam", "{name}-s-md.metrics"
    log: "{name}-s-md.bam.sm.log"
    params: mrkdup_jar=PICARD_DIR + "/MarkDuplicates.jar"
    shell: """
    java -jar {params.mrkdup_jar} \
        INPUT={input} \
        OUTPUT={output[0]} \
        METRICS_FILE={output[1]} \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        REMOVE_DUPLICATES=TRUE \
        VALIDATION_STRINGENCY=SILENT
    """

# Gets rid of the MAPQ should be 0 warning
rule cleansam:
    input: "{name}.bam"
    output: "{name}-c.bam"
    log: "{name}-c.bam.sm.log"
    params: clnsam_jar=PICARD_DIR + "/CleanSam.jar"
    shell: """
    java -jar {params.clnsam_jar} \
        INPUT={input} \
        OUTPUT={output} \
        QUIET=TRUE \
    """

# Create index file for fast bam access
rule index:
    input: "{name}.bam"
    output: "{name}.bam.bai"
    log: "{name}.bam.bai.sm.log"
    shell: """
    samtools index {input}
    """

# Use BEDTools to get coverage of bam file
rule coverage:
    input: "{name}.bam"
    output: "{name}.coverage"
    shell: """
    genomeCoverageBed -ibam {input} > {output}
    """

# Calculate mean coverage per contig given BEDTools output
# TODO: doesn't work, {} should be escaped somehow
rule mean_coverage_per_contig:
    input: "{name}.coverage"
    output: "{name}.mean.coverage.per.contig"
    shell: """
    awk 'BEGIN {pc=""}
    {
        c=$1;
        if (c == pc) {
            cov=cov+$2*$5;
        } else {
            print pc,cov;
            cov=$2*$5;
            pc=c
        }
    } END {print pc,cov}' {input} | tail -n +2 > {output}
    """

rule clean:
    run:
        shell("rm -rf " + " ".join([SPLIT_OUT, MAP_IM_OUT, MERGE_OUT]))
        shell("rm -rf " + " ".join(expand("ref.fa.{ext}",ext=["bwt","pac","fai","sa","amb","ann"])))
