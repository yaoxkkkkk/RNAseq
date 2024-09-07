import os

# 提取文件名的基部分（去除路径和扩展名）
ref_basename = os.path.splitext(os.path.basename(config["ref"]))[0]
fastq_suffix = config.get("fastq_suffix", ".fq.gz")

qualified_quality_phred = config.get("qualified_quality_phred", 20)
unqualified_percent_limit = config.get("unqualified_percent_limit", 40)
trim_front = config.get("trim_front", 10)

rule all:
    input:
        expand("genome_index/{ref_basename}.{ext}", ref_basename=ref_basename, ext=["1.ht2", "2.ht2", "3.ht2", "4.ht2", "5.ht2", "6.ht2", "7.ht2", "8.ht2"]),
        expand("clean_data/{sample}_{pair}_clean.fq.gz", sample=config["sample"], pair=["1", "2"]),
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
        &> {log}
        """

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

rule HISAT2_map:
    input:
        expand("genome_index/{ref_basename}.{ext}", ref_basename=ref_basename, ext=["1.ht2", "2.ht2", "3.ht2", "4.ht2", "5.ht2", "6.ht2", "7.ht2", "8.ht2"]),
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
