# 3'GAmES (3´-terminal Gene Annotation from mRNA 3´End Sequencing datasets)
3' GAmES is a pipeline to refine and extend mRNA 3' end annotations using 3' end sequencing datasets. 

## Installation 

### Clone from github

git clone https://github.com/poojabhat1690/3-GAmES.git

Please make sure the singularity version you have is > 3.0. 

	cd 3GAmES/bin
	run the script "getDependencies.sh"

## Description
the main steps involved in the pipline are outlined :
[The main steps of the pipeline are outlined here](primingSites/flowchart.pdf).
A description of the different steps outlined can be found in the wiki.



## Requirements
After cloning 3'GAmES and pulling the dependencies, the last thing required are ENSEMBL and refSeq annotations... 

### Annotations required
To start, 3' GAmES requrires the dependecies as singularity images and an annotation set from refSeq and ENSEMBL. 

5 different annotation files are required. Please reatain the file names as below.

	1. refSeq_mrna_utrsPresent.bed - refSeq 3' UTR annotations  
	2. proteinCoding_annotatedUTRs.bed - ENSEMBL 3' UTR annotations 
	3. exonInfo_proteinCodingGenes.bed - ENSEMBL exon annotations 
	4. intronInfo_proteinCodingGenes.bed - ENSEMBL intron annotations 
	5. transcriptStartsAndEnds_all.txt - ensembl transcript annotation  	

All the annotations should be tab separated and should contain the following information (without column headers): chromosome, start, end , geneName, score, strand, transcript id
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
	
	-p link to the 3GAmES folder



## Output description
The final ouput and all intermediate files are organized in the follwing folders:

#### final90percent :
 1. ends_greater90percent_intergenic_n100 :  contains the high condience mRNA 3' ends 
 2. allAnnotations.bed :  250nt counting windows (overlapping counting windows merged), used to count quantSeq reads.
 3. countingWindows_transcriptionalOutput.bed : genomic loci to filter multimappers using SLAMdunk. These include all counting windows + all 3' UTRs + all extended counting windows. 
 4. onlyIntergenic_90percent_n100: list of the intergenic counting windows created using presence of continuous RNAseq signal.

#### PASplots :
 contains nucleotide profiles for priming sites overlapping with annotations, sparated by presence or absence of the poly A signal (PAS) and separated by downstream genomic A content. 



## Contact
Pooja Bhat (pooja.bhat@imba.oeaw.ac.at)        


