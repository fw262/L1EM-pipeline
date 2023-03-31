# use STAR to generate alignment BAM files
import pdb
configfile:'config.yaml'

import pandas as pd
import os
import re
import glob

text_file=open(config['Samples'])
SRRIDs=text_file.read().splitlines()
print(SRRIDs)
#SRRIDs=SRRIDs[1:100]

CORES=config['CORES']
L1EM=config['L1EM']
GENOMEREF=config['GENOMEREF']
ROOT_DIR=config['ROOT_DIR']

rule all:
        #input: expand('{PIPELINE_MAJOR}/{sample}/Aligned.out.bam', PIPELINE_MAJOR=config['PIPELINE_MAJOR'],sample=SRRIDs)
        input: expand('{PIPELINE_MAJOR}/{sample}/L1EM/full_counts.txt', PIPELINE_MAJOR=config['PIPELINE_MAJOR'],sample=SRRIDs)
        #input:expand('{ROOT_DIR}/{PIPELINE_MAJOR}/{sample}/L1EM/full_counts.txt', ROOT_DIR=config['ROOT_DIR'],PIPELINE_MAJOR=config['PIPELINE_MAJOR'],sample=SRRIDs)

rule convertFQ:
	input:	read1 = config['DATADIR']+'/{sample}_1.fasta',
		read2 = config['DATADIR']+'/{sample}_2.fasta'
	output:	read1 = config['DATADIR']+'/{sample}_1.fastq',
                read2 = config['DATADIR']+'/{sample}_2.fastq',
	shell:
		"""
		seqtk seq -F "F" {input.read1} > {output.read1}
		seqtk seq -F "F" {input.read2} > {output.read2}
		"""


rule STAR:
	input:	read1 = config['DATADIR']+'/{sample}_R1_001.fastq.gz',
		read2 = config['DATADIR']+'/{sample}_R2_001.fastq.gz',
		index = config['STAR_IND']
	output:	ou1 = '{path}/{sample}/Aligned.sortedByCoord.out.bam',
		ou2 = '{path}/{sample}/Aligned.out.bam'
	params: outputPref = '{path}/{sample}/'
	threads: CORES
	shell:	"""
		STAR --runThreadN {CORES} \
		--readFilesCommand zcat \
		--genomeDir {input.index} \
		--readFilesIn {input.read1} {input.read2} \
		--genomeLoad NoSharedMemory \
		--outSAMtype BAM Unsorted SortedByCoordinate \
		--outSAMunmapped Within KeepPairs \
		--outFileNamePrefix {params.outputPref} \
		--quantMode GeneCounts
		"""

rule indexBam:
        input:  sorted = '{path}/{sample}/Aligned.sortedByCoord.out.bam'
        output: ou1 = '{path}/{sample}/Aligned.sortedByCoord.out.bam.bai'
        threads: CORES
        shell:  """
                samtools index {input.sorted}
		"""

rule L1EM:
        input:  bam = '{path}/{sample}/Aligned.sortedByCoord.out.bam', ind = '{path}/{sample}/Aligned.sortedByCoord.out.bam.bai'
        params:	L1EMPath = L1EM, genomeRef = GENOMEREF, dir = ROOT_DIR+'/{path}/{sample}/L1EM', rootDir=ROOT_DIR
	output: '{path}/{sample}/L1EM/full_counts.txt'
        threads: CORES
        shell:  """
		mkdir -p {params.dir}
		cd {params.dir}
		rm -rf *
		inputFull={params.rootDir}/{input.bam}
		bash -e {params.L1EMPath}/run_L1EM.sh \
		$inputFull \
		{params.L1EMPath} \
		{params.genomeRef}
		cd ..
		ln -s {params.dir}/full_counts.txt ./
                """
