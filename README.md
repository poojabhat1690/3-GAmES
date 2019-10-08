# 3'GAmES (3´-terminal Gene Annotation from mRNA 3´End Sequencing datasets)
3' GAmES is a pipeline to refine and extend mRNA 3' end annotations using 3' end sequencing datasets. 

## Author
Pooja Bhat (pooja.bhat@imba.oeaw.ac.at)

## Installation (change this). 

#### Clone from github

git clone https://github.com/poojabhat1690/refine-pipeline.git
cd pipeline/pre-processing/

All dependencies required for 3' GAmES are provided as 3 separate singularity modules. 
 1. all dependencies 
 2. R and R packages required
 3. slamdunk
Please make sure the singularity version you have is > 3.0. 

## Quickstart
 the main steps involved in the pipline are outlined : 
[Contribution guidelines for this project](primingSites/flowchart.pdf)
### Annotations required
To start, 3' GAmES requrires the dependecies as singularity images and an annotation set from refSeq and ENSEMBL. 

6 different annotation files are required. Please reatain the file names as below.
	1. refSeq_mrna_utrsPresent.bed - refSeq 3' UTR annotations  
	2. proteinCoding_annotatedUTRs.bed - ENSEMBL 3' UTR annotations 
	3. exonInfo_proteinCodingGenes.bed - ENSEMBL exon annotations 
	4. intronInfo_proteinCodingGenes.bed - ENSEMBL intron annotations 
	
The annotations should be tab separated and have the following columns : chromosome, start, end , geneName, score, strand, transcript id. 
One way to obtain the annotations is from the UCSC table browser.

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

### Running 3'GAmES

The script required to run the whole pipeline is beforeMapping.new.sh

beforeMapping.new.sh -a [adapter] -i [input directory] -o [output directory] -g [genome file] -t [threshold for priming sites]
-u [ucscDir] -e [ensemblDir] -m [mode rnaseq p/s/S] -c [condition]
 
 
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
                
 ## Prerequisites: 



#### Annotations from refSeq, downloaded from the UCSC table browser manually. 

for 3' UTR annotations from UCSC genome browser 
             1. Form UCSC table browser, select:
                a. clade = mammal 
                b. genome = mouse
                c. assembly = Dec. 2011 (GRCm38/mm10)
                d. group = genes and gene prediction
                e. track = refSeq genes
                f. table = refGene
              2. select bed format and 3' UTR. 
      

## Output description
The final ouput and all intermediate files are organized in the follwing folders:

#### final90percent :
 1. ends_greater90percent_intergenic_n100 :  contains the high condience mRNA 3' ends 
 2. allAnnotations.bed :  250nt counting windows (overlapping counting windows merged), used to count quantSeq reads.
 3. countingWindows_transcriptionalOutput.bed : genomic loci to filter multimappers using SLAMdunk. These include all counting windows + all 3' UTRs + all extended counting windows. 4. onlyIntergenic_90percent_n100: list of the intergenic counting windows created using presence of continuous RNAseq signal.

#### PASplots :
 contains nucleotide profiles for priming sites overlapping with annotations, sparated by presence or absence of the poly A signal (PAS) and separated by downstream genomic A content. 



        


