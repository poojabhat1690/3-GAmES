# 3'GAmES (3´-terminal Gene Annotation from mRNA 3´End Sequencing datasets)
3' GAmES is a pipeline to refine and extend mRNA 3' end annotations using 3' end sequencing datasets. 

## Installation 

### Clone from github

git clone https://github.com/poojabhat1690/3-GAmES.git

Please make sure the singularity version you have is > 3.0. 

	cd 3-GAmES/bin
	run the script "getDependencies.sh"

## Description
the main steps involved in the pipeline are outlined :
[The main steps of the pipeline are outlined here](primingSites/flowchart.pdf).
A description of the different steps outlined can be found in the wiki.



## Requirements
After cloning 3'GAmES and pulling the dependencies
1. Organize the input data: 3'GAmES requires an input folder that contains. 

		a. A folder named 'quantseq', which contains fastq files in formal *.fq, *.fastq or *.gz
		b. A folder namesd 'rnaseq' [optional], which contains mapped RNAseq files, and the corresponding indices. The data in this folder is needed for extending the analysis to intergenic regions. 
The pipeline can be run without RNAseq data - in which case, please provide only a 'quantseq' in the input directory. In this mode, intergenic ends will not be identified.

	 
2. create ENSEMBL and refSeq annotation files (descripton provided below)

### Annotations required
To start, 3' GAmES requires  the dependencies  as singularity images and an annotation set from refSeq and ENSEMBL. 

5 different annotation files are required. Please retain the file names as below.

	1. refSeq_mrna_utrsPresent.bed - refSeq 3' UTR annotations  
	2. proteinCoding_annotatedUTRs.bed - ENSEMBL 3' UTR annotations 
	3. exonInfo_proteinCodingGenes.bed - ENSEMBL exon annotations 
	4. intronInfo_proteinCodingGenes.bed - ENSEMBL intron annotations 
	5. transcriptStartsAndEnds_all.txt - ensembl transcript annotation  	

All the annotations should be tab separated and should contain the following information (without column headers): chromosome, start, end, geneName, score, strand, transcript id.

Please also include the transcript biotype (protein coding, non-coding etc) as the 8th column of  transcriptStartsAndEnds_all.txt


All the above annotations can be obtained from the UCSC table browser:

	For example:
		
	for 3' UTR annotations from (refSeq) UCSC genome browser                                                                                                                                    
             1. Form UCSC table browser, select:
                a. clade = mammal 
                b. genome = mouse
                c. assembly = Dec. 2011 (GRCm38/mm10)
                d. group = genes and gene prediction
                e. track = refSeq genes
                f. table = refGene
              2. select bed format and 3' UTR. 

## Quick start -  Running 3'GAmES

The script required to run the whole pipeline is run_3GAmES.sh. You will find this in the 3GAmES/bin/

run_3GAmES.sh -a [adapter] -i [input directory] -o [output directory] -g [genome file] -t [threshold for priming sites]
-u [ucscDir] -e [ensemblDir] -m [mode rnaseq p/s/S] -c [condition] -p [path]
 
 
 	-a 3' adapter sequences that has to be removed using cutadapt
 
 	-i input directory containing two folders named - quantseq, rnaseq
                   quantseq: contains *.fastq, *.fq.gz, *fq  files. 
                   rnaseq: mapped, sorted and indexed bam files. 
 
 	-o Output directory
 
 	-t threshold of the number of reads to define a priming site.
 
 	-u ucsc directory containing annotations downloaded from ucsc table browser. 
 
 	-e ensembl directory containing ensembl annotations ontained from biomart. 
 
 	-m mode of counting for RNAseq coverage, derived from bedtools multicov (s: counting on the same strand, 
            p: paired end reads, S: counting on the opposite strand
 
 	-c condition of sample (example: timepoint or organism)
	
	-p path to the 3-GAmES folder

	-s step after which pipeline will stop [options:preprocessing, primingsites, intergenicends, all]

example : run_3GAmES.sh -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -i /scratch/3-GAmES/xcondition_input/ -o /scratch/xcondition_output/ -g /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa -t 2  -e /scratch/3-GAmES/annotations_all/ -m p -c trialSample -p /scratch/3-GAmES -s all

## Output description
The final output and all intermediate files are organized in the following folders:

#### finalEnds:
 1. highConfidenceEnds.bed :  contains the high confidence mRNA 3' ends 
 2. countingWindows.bed :  250nt counting windows (overlapping counting windows merged), used to count quantSeq reads.
 3. countingWindows_transcriptionalOutput.bed : genomic loci to filter multimappers using SLAMdunk. These include all counting windows + all 3' UTRs + all extended counting windows. 
 4. highConfidenceIntergenicEnds.bed: list of the intergenic counting windows created using presence of continuous RNAseq signal.

#### PASplots:
Contains nucleotide profiles for priming sites overlapping with annotations, sparated by presence or absence of the poly A signal (PAS) and separated by downstream genomic A content. 



## Contact
Pooja Bhat (pooja.bhat@imba.oeaw.ac.at)        


