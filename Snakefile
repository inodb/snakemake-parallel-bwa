# Split parameters
SPLIT_OUT="split-out"
LINES_PER_FILE=400000

# aln parameters
MAP_IM_OUT="map-intermediate-out"

# merge parameters
MERGE_OUT="map"

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
    shell: """
	bwa index {input}
    """

rule aln:
    input: pairn=SPLIT_OUT + "/pair{pair_nr}.fastq.split.{id}", refbwt="ref.fa.bwt", ref="ref.fa"
    threads: 1
    params: out_dir=MAP_IM_OUT
    output: MAP_IM_OUT + "/pair{pair_nr}.fastq.split.{id}.sai"
    shell: """
    mkdir -p {params.out_dir}
	bwa aln -t {threads} {input.ref} -{wildcards.pair_nr} {input.pairn} > {output}
    """

rule sampe:
    input: sai1=MAP_IM_OUT + "/pair1.fastq.split.{id}.sai", sai2=MAP_IM_OUT + "/pair2.fastq.split.{id}.sai", pair1=SPLIT_OUT + "/pair1.fastq.split.{id}", pair2=SPLIT_OUT + "/pair2.fastq.split.{id}", refbwt="ref.fa.bwt", ref="ref.fa"
    params: out_dir=MAP_IM_OUT
    output: MAP_IM_OUT + "/pair.fastq.split.{id}.sam"
    shell: """
    mkdir -p {params.out_dir}
	bwa sampe {input.ref} {input.sai1} {input.sai2} {input.pair1} {input.pair2} > {output}
    """

rule fai:
    input: "ref.fa"
    output: "ref.fa.fai"
    shell: """
    samtools faidx {input}
    """

rule samtobam:
    input: "{name}.sam", "ref.fa.fai"
    output: "{name}.bam"
    shell: """
    samtools view -bt {input[0]} {input[1]} > {output}
    """

rule merge:
    input: dynamic(MAP_IM_OUT + "/pair.fastq.split.{splitid}.bam")
    params: out_dir=MERGE_OUT
    output: MERGE_OUT + "/pair_ref.bam"
    shell: """
    mkdir -p {params.out_dir}
    samtools merge {output} {input}
    """

rule clean:
    shell: "rm -rf " + SPLIT_OUT
