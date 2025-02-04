import os

configfile: "RNAseq_config_featureCounts.yaml"

# 提取文件名的基部分（去除路径和扩展名）
ref_basename = os.path.splitext(os.path.basename(config["ref"]))[0]
fastq_suffix = config.get("fastq_suffix", ".fq.gz")

rule all:
    input:
        expand("mapping/{sample}.sorted.bam", sample=config["sample"]),
        "counts/counts.txt"

rule HISAT2_index:
    input:
        reference_genome=config["ref"]
    output:
        "genome_index/{ref_basename}.1.ht2",
        "genome_index/{ref_basename}.2.ht2",
        "genome_index/{ref_basename}.3.ht2",
        "genome_index/{ref_basename}.4.ht2",
        "genome_index/{ref_basename}.5.ht2",
        "genome_index/{ref_basename}.6.ht2",
        "genome_index/{ref_basename}.7.ht2",
        "genome_index/{ref_basename}.8.ht2"
    log:
        "logs/index/hisat2_index_{ref_basename}.log"
    shell:
        """
        hisat2-build \
        {input.reference_genome} \
        genome_index/{wildcards.ref_basename} \
        2> {log}
        """

rule QualityControlfastp:
    input:
        f"raw_data/{{sample}}_clean_1{fastq_suffix}",
        f"raw_data/{{sample}}_clean_2{fastq_suffix}"
    output:
        "clean_data/{sample}_1_clean.fq.gz",
        "clean_data/{sample}_2_clean.fq.gz",
        "logs/fastp/fastp_report/{sample}.fastp.html"
    threads: 2
    params:
        qualified_quality_phred=config["qualified_quality_phred"],
        unqualified_percent_limit=config["unqualified_percent_limit"],
        trim_front=config["trim_front"]
    log:
        "logs/fastp/{sample}.log"
    shell:
        """
        fastp \
        --thread {threads} \
        -i {input[0]} \
        -I {input[1]} \
        -o {output[0]} \
        -O {output[1]} \
        -h {output[2]} \
        -j /dev/null \
        -q {params.qualified_quality_phred} \
        -u {params.unqualified_percent_limit} \
        -f {params.trim_front} \
        2> {log}
        """

rule HISAT2_map:
    input:
        hisat_index=expand("genome_index/{ref_basename}.{ext}", ref_basename=ref_basename, ext=["1.ht2", "2.ht2", "3.ht2", "4.ht2", "5.ht2", "6.ht2", "7.ht2", "8.ht2"]),
        R1="clean_data/{sample}_1_clean.fq.gz",
        R2="clean_data/{sample}_2_clean.fq.gz"
    output:
        "mapping/{sample}.sorted.bam"
    threads: 8
    params:
        hisat_index=f"genome_index/{ref_basename}"
    log:
        "logs/hisat2/hisat2_map_{sample}.log"
    shell:
        """
        hisat2 \
        -p {threads} \
        -x {params.hisat_index} \
        -1 {input.R1} \
        -2 {input.R2} \
        | samtools sort -@ {threads} -o {output} \
        2> {log}
        """

rule Subread_featureCounts:
    input:
        bam=expand("mapping/{sample}.sorted.bam", sample=config["sample"]),
        annotation=config["annotation"]
    output:
        "counts/counts.txt"
    threads: 8
    params:
        feature=config["feature_type"],
        attribution=config["attribution"]
    log:
        "logs/subread/featureCounts.log"
    shell:
        """
        featureCounts \
        -T {threads} \
        -p \
        --countReadPairs \
        -B \
        -C \
        -t {params.feature} \
        -g {params.attribution} \
        -a {input.annotation} \
        -o {output} \
        {input.bam} \
        2> {log}
        """
