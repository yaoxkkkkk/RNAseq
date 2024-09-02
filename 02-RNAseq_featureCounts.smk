import os

# 提取文件名的基部分（去除路径和扩展名）
ref_basename = os.path.splitext(os.path.basename(config["ref"]))[0]
fastq_suffix = config.get("fastq_suffix", ".fq.gz")

qualified_quality_phred = config.get("qualified_quality_phred", 20)
unqualified_percent_limit = config.get("unqualified_percent_limit", 40)
trim_front = config.get("trim_front", 10)

rule all:
    input:
        expand("clean_data/{sample}_1_clean.fq.gz", sample=config["sample"]),
        expand("clean_data/{sample}_2_clean.fq.gz", sample=config["sample"]),
        expand("mapping/{sample}.sorted.bam", sample=config["sample"]),
        "counts/counts.txt"

rule QualityControlfastp:
    input:
        f"raw_data/{{sample}}_R1{fastq_suffix}",
        f"raw_data/{{sample}}_R2{fastq_suffix}"
    output:
        "clean_data/{sample}_1_clean.fq.gz",
        "clean_data/{sample}_2_clean.fq.gz",
        "logs/fastp/fastp_report/{sample}.fastp.html"
    log:
        "logs/fastp/{sample}.log"
    threads: 2
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
        -q {qualified_quality_phred} \
        -u {unqualified_percent_limit} \
        -f {trim_front} \
        &> {log}
        """

rule Hisat2_map:
    input:
        "clean_data/{sample}_1_clean.fq.gz",
        "clean_data/{sample}_2_clean.fq.gz"
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
        -1 {input[0]} \
        -2 {input[1]} \
        --very-sensitive \
        --dta \
        | samtools sort -@ {threads} -o {output} \
        &> {log}
        """

rule Subread_featureCounts:
    input:
        bam=expand("mapping/{sample}.sorted.bam", sample=config["sample"]),
        annotation=config["annotation"]
    output:
        "counts/counts.txt"
    threads: 8
    params:
        feature="exon",
        attribute="gene_id"
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
        -g {params.attribute} \
        -a {input.annotation} \
        -o {output} \
        {input.bam}
        """
