import os

# 提取文件名的基部分（去除路径和扩展名）
ref_basename = os.path.splitext(os.path.basename(config["ref"]))[0]
fastq_suffix = config.get("fastq_suffix", ".fq.gz")

qualified_quality_phred = config.get("qualified_quality_phred", 20)
unqualified_percent_limit = config.get("unqualified_percent_limit", 40)
trim_front = config.get("trim_front", 10)

rule all:
    input:
        expand("mapping/{sample}.sorted.markdup.bam", sample=config["sample"]),
        

rule QualityControlfastp:
    input:
        f"raw_data/{{sample}}_1{fastq_suffix}",
        f"raw_data/{{sample}}_2{fastq_suffix}"
    output:
        "clean_data/{sample}_1_clean.fq.gz",
        "clean_data/{sample}_2_clean.fq.gz",
        "clean_data/{sample}.fastp.html"
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
        -q {qualified_quality_phred} \
        -u {unqualified_percent_limit} \
        -f {trim_front} \
        &> {log}
        """

rule Hisat2_map:
    input:
        "clean_data/{sample}.1_clean.fq.gz",
        "clean_data/{sample}.2_clean.fq.gz",
        "genome_index/{ref_basename}.1.ht2",
        "genome_index/{ref_basename}.2.ht2",
        "genome_index/{ref_basename}.3.ht2",
        "genome_index/{ref_basename}.4.ht2",
        "genome_index/{ref_basename}.5.ht2",
        "genome_index/{ref_basename}.6.ht2",
        "genome_index/{ref_basename}.7.ht2",
        "genome_index/{ref_basename}.8.ht2"
    output:
        "mapping/{sample}.sorted.bam"
    threads: 8
    params:
        hisat_index="genome_index/{ref_basename}"
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
        | samtools view -Sb \
        | samtools sort > \
        {output} \
        &> {log}
        """

rule Stringtie_assembly:
    input:
        bam = "mapping/{sample}.sorted.bam",
        annotation = config["annotation"]
    output:
        "gtf/{sample}.gtf"
    threads: 8
    log:
        "logs/stringtie/stringtie_assembly_{sample}.log"
    shell:
        """
        stringtie \
        -p {threads} \
        -G {input.annotation} \
        -o {output} \
        {input.bam} \
        &> {log}
        """

rule Stringtie_GTF_merge:
    input:
        bam = "mapping/{sample}.sorted.bam",
        annotation = config["annotation"]
    output:
        "gtf/{sample}.gtf"
    threads: 8
    log:
        "logs/stringtie/stringtie_assembly_{sample}.log"
    shell:
        """
        stringtie \
        -p {threads} \
        -G {input.annotation} \
        -o {output} \
        {input.bam} \
        &> {log}
        """