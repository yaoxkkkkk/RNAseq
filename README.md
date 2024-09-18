# RNAseq-snakemake-pipeline

## Dependent Software

- fastp
- HISAT2
- samtools
- featureCounts
- stringtie

## What the pipeline does

- Genome file index creation (HISAT2)
- RNA-seq reads quality control
- RNA-seq reads map to reference genome
- Transcripts abundance quantification

## What to input

- Reference genome sequence file
- Reference genome annotation file (in `.gtf` format)
- RNA-seq fastq files

## What to output

Gene expression count matrix

## Usage

### 1. Prepare your working directory

```shell
├── script
│   └── snake_pipeline
├── raw_data
├── genome_index
└── logs
```

Please storage your resequence data in `raw_data/` folder and genome file in `genome_index/` folder. Script files, pipeline files and configuration files can be stored in the way you are used to.

### 2. Prepare the config file

The config file needs to be at the same folder of snakefile.

2.1 Move the genome file to `genome_index/` folder and add the genome fasta file absolute path like:

```shell
# Absolute path to the genome fasta file
ref: "/workingdir/genome_index/genome.fasta" 
```

2.2 Sometimes the fastq files may be ended with `.fastq.gz` or `.fq.gz`, specify the suffix of the fastq files if it's necessary.

```shell
# Fastq file suffix
fastq_suffix: " " # Default value is ".fq.gz"
```

2.3 Fill in the name of the samples. The samples name need to be filled with specific format like:

```shell
# Sample list, samples' name should start with letters.
sample:
    - sample1
    - sample2
    - sample3
    - sample4
    - ...
    - samplen
```

You can use following command to add sample list to the config file if you have a sample list txt file (for example `sample.list`):

```shell
# sample.list
sample1
sample2
sample3
sample4

# Add samples to the config file:
sed 's/^/    - /' sample.list >> ${working_dir}/RNAseq_config.yaml
```

### 3. Submit the pipeline to HPC cluster

For example:

```shell
snakemake \
	--snakefile ${working_dir}/00-script/snake_pipeline/${snakemake_file} \
	-d ${working_dir} \
    	--configfile ${working_dir}/RNAseq_config.yaml \
	--cores ${cores_num} \
	--rerun-incomplete \
	--latency-wait 360 \
	--keep-going
```
