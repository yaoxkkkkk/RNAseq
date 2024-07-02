import os

# Extract basename of fasta file (part between path and extension name)
ref_basename = os.path.splitext(os.path.basename(config["ref"]))[0]

rule all:
    input:
        expand("genome_index/{ref_basename}.{ext}", ref_basename=ref_basename, ext=["1.ht2", "2.ht2", "3.ht2", "4.ht2", "5.ht2", "6.ht2", "7.ht2", "8.ht2"])

rule hisat2_index:
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