# Variant Calling Workflow

This Snakemake workflow is designed for variant calling using a series of rules to perform quality control, read trimming, mapping, and variant calling. The workflow is configured using a YAML file (`config/config.yaml`) and a samples file (`config/samples.tsv`).

## Prerequisites

Before running the workflow, make sure you have the following prerequisites installed:

- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Fastp](https://github.com/OpenGene/fastp)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [Samtools](http://www.htslib.org/)
- [Bcftools](http://samtools.github.io/bcftools/)

Conda environments for these tools are defined in the `envs` directory.

## Configuration

- Edit the configuration file (`config/config.yaml`) to set the paths to the genome and other parameters.

- Define your samples in the samples file (`config/samples.tsv`) with columns `sample_name` and `fastq`.

## Workflow Execution

To execute the workflow, run the following command:

```bash
snakemake --use-conda
```

## Workflow Structure (Continued)

### Rules Description

- **QualityControl_1**: This rule performs initial quality control using FastQC on the raw reads.

- **trim**: Trims reads using Fastp to generate trimmed fastq files.

- **QualityControl_2**: Performs quality control on trimmed reads using FastQC.

- **bowtie_index**: Builds the Bowtie2 index for the reference genome.

- **map_reads**: Maps trimmed reads to the reference genome using Bowtie2.

- **sam_to_bam**: Converts the mapped reads to BAM format using Samtools.

- **index_bam**: Indexes the BAM files using Samtools.

- **call_variants**: Calls variants using Bcftools.

### Output Directory Structure

The workflow generates the following directory structure:

- **QC_BeforeTrim**: Initial quality control results before trimming.

- **QC_AfterTrim**: Quality control results after trimming.

- **00_mapped_reads**: Intermediate directory containing unmapped SAM files.

- **01_called_variants**: Final directory containing called variants in VCF format.

### Running Individual Rules

If you want to run a specific rule or set of rules, you can use the following command:

```bash
snakemake <rule_name> --use-conda

