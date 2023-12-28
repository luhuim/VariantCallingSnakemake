import pandas as pd                                                                                                                    
import re  

configfile: 'config/config.yaml'                                                                                                                                                                                                                                              
sample_df = (pd.read_csv(config['samples'], sep='\t',dtype={'sample_name':str, 'fastq':str}).set_index('sample_name', drop=False))  

rule all:                                                                                                                                      
	input:                                                                                                                                         
		expand('results/01_called_variants/{sample}.vcf',sample=sample_df.sample_name) 

rule QualityControl_1:                                                                                                                         
	input:                                                                                                                                         
		"./resources/{sample}Reads.fastq"                                                         
	output:                                                                                                                                        
		before= directory("resources/QC_BeforeTrim/{sample}"),                                                                                                  
	conda:                                                                                                                                         
		"envs/fastqc.yaml"                                                                                                             
	shell:    
		"""             
		# mkdir {output.before} ;\                                                                                                                     
		fastqc {input} -o {output.before}		
		"""

rule trim:
	input:
		rules.QualityControl_1.input
	output:
		"./resources/{sample}Reads_trimmed.fastq"
	conda:
		"envs/fastqc.yaml"
	shell:
		"fastp -i {input} -o {output}"


rule QualityControl_2:                                                                                                                         
	input:                                                                                                                                         
		rules.trim.output                                                          
	output:                                                                                                                                        
		after= directory("resources/QC_AfterTrim/{sample}"),                                                                                            
	conda:                                                                                                                                         
		"envs/fastqc.yaml"                                                                                                             
	shell:    
		"""             
		# mkdir {output.after};\
		fastqc {input} -o {output.after}
		"""

rule bowtie_index:
	input:
		genome = config['genome']
	output:
		multiext(config['genome'],
		".1.bt2", ".2.bt2", ".3.bt2",
		".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
	conda: 'envs/bowtie2.yaml'
	shell:
		"bowtie2-build {input.genome} {input.genome} "

rule map_reads:
	input:
		idx = rules.bowtie_index.output,
		#reads = lambda wildcards: sample_df.loc[wildcards.sample, 'fastq']
		#reads = rules.QualityControl.output.trim #Here is the problem
		reads = lambda wildcards: rules.trim.output
	output:
		temp('results/00_mapped_reads/{sample}.unsorted.sam')
	params:
		idx = config['genome']
	conda: 'envs/bowtie2.yaml'
	shell:
		"bowtie2 -x {params.idx} -U {input.reads} -S {output} "

rule sam_to_bam:
	input:
		rules.map_reads.output
	output:
		'results/00_mapped_reads/{sample}.bam'
	threads: 2
	conda: 
		'envs/htslib.yaml'
	shell:
		"samtools sort -@ {threads} -o {output} {input} "
rule index_bam:
	input:
		rules.sam_to_bam.output
	output:
		'results/00_mapped_reads/{sample}.bam.bai'
	conda: 
		'envs/htslib.yaml'
	shell:
		'samtools index {input} '

rule call_variants:
	input:
		rules.index_bam.output,
		aligned_reads = rules.sam_to_bam.output,
		genome = config['genome']

	output:
		'results/01_called_variants/{sample}.vcf'
	conda:
		'envs/htslib.yaml'
	shell:
		'bcftools mpileup -Ou '
		'-f {input.genome} '
		'{input.aligned_reads} '
		'| bcftools call -m -v '
		'> {output} '
    
