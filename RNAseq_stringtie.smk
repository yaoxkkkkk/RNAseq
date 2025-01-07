import os

configfile: "RNAseq_config.yaml"

# 提取文件名的基部分（去除路径和扩展名）
ref_basename = os.path.splitext(os.path.basename(config["ref"]))[0]
fastq_suffix = config.get("fastq_suffix", ".fq.gz")

rule all:
    input:
        expand("mapping/{sample}.sorted.bam", sample=config["sample"]),
        "count/tab/{sample}.tab",
        "count/{sample}.count.gtf",
        "gtf/stringtie.transdecoder.clean.cds.fasta",
        "gtf/stringtie.transdecoder.clean.pep.fasta"


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
        --dta \
        | samtools sort -@ {threads} -o {output} \
        2> {log}
        """

rule StringtieAssembly:
    input:
        bam="mapping/{sample}.sorted.bam",
        annotation=config["annotation"]
    output:
        gtf="gtf/{sample}.gtf"
    threads: 8
    log:
        "logs/stringtie/stringtie_assembly_{sample}.log"
    shell:
        """
        stringtie \
        -p {threads} \
        -G {input.annotation} \
        -o {output.gtf} \
        {input.bam} \
        2> {log}
        """

rule ExtractGTFlist:
    input:
        expand("gtf/{sample}.gtf", sample=config["sample"])
    output:
        "gtf/gtf.list"
    params:
        gtf_dir="gtf/"
    shell:
        """
        find {params.gtf_dir} -name '*.gtf' > {output}
        """
        
rule StringtieGTFmerge:
    input:
        gtflist="gtf/gtf.list",
        annotation=config["annotation"]
    output:
        "gtf/stringtie_merged.gtf"
    threads: 8
    log:
        "logs/stringtie/stringtie_gtf_merge.log"
    shell:
        """
        stringtie \
        --merge \
        -p {threads} \
        -G {input.annotation} \
        -o {output} \
        {input.gtflist} \
        2> {log}
        """
        
rule Modify_Merged_GTF:
    input:
        mstrgscript=config["mstrg_prep_revise_path"],
        mergedgtf="gtf/stringtie_merged.gtf"
    output:
        "gtf/stringtie_merged.mod.gtf"
    log:
        "logs/stringtie/merged_gtf_modify.log"
    shell:
        """
        perl {input.mstrgscript} \
        {input.mergedgtf} \
        | sed 's/ref_gene_id.*//' \
        1> {output} \
        2> {log}
        """
        
rule GTF2fasta:
    input:
        modifiedgtf="gtf/stringtie_merged.mod.gtf",
        genomefile=config["ref"]
    output:
        "gtf/stringtie.exon.fasta"
    conda:
        config["conda_env"]["transdecoder_conda"]
    log:
        "logs/stringtie/gtf_genome_to_cdna_fasta.log"
    shell:
        """
        gtf_genome_to_cdna_fasta.pl \
        {input.modifiedgtf} \
        {input.genomefile} \
        1> {output} \
        2> {log}
        """
        
rule GTF2AlignGFF3:
    input:
        "gtf/stringtie_merged.mod.gtf"
    output:
        "gtf/stringtie.alignment.gff3"
    conda:
        config["conda_env"]["transdecoder_conda"]
    log:
        "logs/stringtie/gtf_to_alignment_gff3.log"
    shell:
        """
        gtf_to_alignment_gff3.pl \
        {input} \
        1> {output} \
        2> {log}
        """
        
rule TransDecoder_LongOrfs:
    input:
        "gtf/stringtie.exon.fasta"
    output:
        temp(directory("stringtie.exon.fasta.transdecoder_dir/"))
    conda:
        config["conda_env"]["transdecoder_conda"]
    log:
        "logs/stringtie/TransDecoder.LongOrfs.log"
    shell:
        """
        TransDecoder.LongOrfs \
        -t {input} \
        2> {log}
        """
        
rule TransDecoder_Predict:
    input:
        "gtf/stringtie.exon.fasta",
        "stringtie.exon.fasta.transdecoder_dir/"
    output:
        temp(expand("stringtie.exon.fasta.transdecoder.{ext}", ext=["bed", "cds", "gff3", "pep"]))
    conda:
        config["conda_env"]["transdecoder_conda"]
    log:
        "logs/stringtie/TransDecoder.Predict.log"
    shell:
        """
        TransDecoder.Predict \
        -t {input[0]} \
        2> {log}
        """
    
rule cdna_alignment_orf_to_genome_orf:
    input:
        "stringtie.exon.fasta.transdecoder.gff3",
        "gtf/stringtie.alignment.gff3",
        "gtf/stringtie.exon.fasta"
    output:
        "gtf/stringtie.transdecoder.gff3"
    conda:
        config["conda_env"]["transdecoder_conda"]
    log:
        "logs/stringtie/cdna_alignment_orf_to_genome_orf.log"
    shell:
        """
        cdna_alignment_orf_to_genome_orf.pl \
        {input[0]} \
        {input[1]} \
        {input[2]} \
        1> {output} \
        2> {log}
        """
        
rule filter_incomplete_gene_coding_models:
    input:
        gff3="gtf/stringtie.transdecoder.gff3",
        genome_file=config["ref"]
    output:
        "gtf/stringtie.transdecoder.clean.gff3"
    container:
        config["singularity"]["AGAT_image"]
    log:
        "logs/stringtie/filter_incomplete_gene_coding_models.log"
    shell:
        """
        agat_sp_filter_incomplete_gene_coding_models.pl \
        --gff {input.gff3} \
        --fasta {input.genome_file} \
        -o {output} \
        2> {log}
        """
        
rule remove_attributes:
    input:
        "gtf/stringtie.transdecoder.clean.gff3"
    output:
        "gtf/stringtie.transdecoder.clean.filter.gff3"
    container:
        config["singularity"]["AGAT_image"]
    log:
        "logs/stringtie/remove_attributes.log"
    shell:
        """
        agat_sp_manage_attributes.pl \
        -gff {input} \
        --att all_attributes \
        -o {output} \
        2> {log}
        """
        
rule extract_cds:
    input:
        gff3="gtf/stringtie.transdecoder.clean.filter.gff3",
        genome_file=config["ref"]
    output:
        "gtf/stringtie.transdecoder.clean.cds.fasta"
    container:
        config["singularity"]["AGAT_image"]
    log:
        "logs/stringtie/extract_cds.log"
    shell:
        """
        agat_sp_extract_sequences.pl \
        --gff {input.gff3} \
        --fasta {input.genome_file} \
        -t cds \
        -o {output} \
        2> {log}
        """
        
rule extract_pep:
    input:
        gff3="gtf/stringtie.transdecoder.clean.filter.gff3",
        genome_file=config["ref"]
    output:
        "gtf/stringtie.transdecoder.clean.pep.fasta"
    container:
        config["singularity"]["AGAT_image"]
    log:
        "logs/stringtie/extract_pep.log"
    shell:
        """
        agat_sp_extract_sequences.pl \
        --gff {input.gff3} \
        --fasta {input.genome_file} \
        -t cds \
        -p \
        -o {output} \
        2> {log}
        """
        
rule gene_quantification:
    input:
        gff3="gtf/stringtie.transdecoder.clean.filter.gff3",
        bam="mapping/{sample}.sorted.bam"
    output:
        countGTF="count/{sample}.count.gtf",
        countTAB="count/tab/{sample}.tab"
    threads: 8
    log:
        "logs/stringtie/{sample}.gene_quantification.log"
    shell:
        """
        stringtie \
        -e \
        -p {threads} \
        -G {input.gff3} \
        -o {output.countGTF} \
        -A {output.countTAB} \
        {input.bam} \
        2> {log}
        """
