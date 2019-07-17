#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
##SBATCH --cpus-per-task=7
#SBATCH --mem-per-cpu=20G
#SBATCH --time=0-24:00:00     # 2 minutes
#SBATCH --job-name=runSLAMdunk
####SBATCH --array=1-4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pooja.bhat@imba.oeaw.ac.at
#SBATCH  --qos=long


PIPELINE=/groups/ameres/Pooja/Projects/refiningPipeline/pipeline/





### reading in the variables passed in the command line... 

while getopts 'a: i: o: g: t: u: e: m: c:' OPTION; do
	  case "$OPTION" in
		      a)
			      avalue="$OPTARG"
			  echo "the adapter trimmed has the sequence $OPTARG"
				          ;;
	
					      i)
						      ivalue="$OPTARG"
						        echo "the input directory containing QuantSeq data is $OPTARG"
							          ;;
										
							o) 
								ovalue="$OPTARG"
								echo "the output dorectory is $OPTARG"
									;;

							
							 g) 
								 genome="$OPTARG"
								 echo "the genome is specified in $OPTARG"
								 ;;

							t) 
								threshold="$OPTARG"
								echo "the threshold to consider priming sites is $OPTARG"
								;;

							u) 
								ucscdir="$OPTARG"
								echo "the ucsc directory specified is $OPTARG"
								;;
							e)
								ensembldir="$OPTARG"
								echo "the ensembl directory specified is $OPTARG"
								;;
							m)
								RNAseqMode="$OPTARG"
								echo "the rnaseq mode is $OPTARG"
								;;
							c)
								Condition="$OPTARG"
								echo "the condition is $OPTARG"
								;;

					  			  ?)
										 echo "script usage: $(basename $0) [-a adapter] [-i input directory] [-o output directory] [-g genome file] [-t threshold for priming sites] [-u ucscDir] [-e ensemblDir] [-m mode rnadeq p/s/S] [-c condition]" >&2
														        exit 1
															      ;;
															        esac
															done
															shift "$(($OPTIND -1))"




if [ "x" == "x$avalue" ]; then
	  echo "-a [adapter] is required"
	    exit
 fi



 if [ "x" == "x$ivalue" ]; then
	echo "-i [input directory] is required"
	exit
fi



 if [ "x" == "x$ovalue" ]; then
	         echo "-o [output directory] is required"
		         exit
fi
		     









#### i value is the input directory ... this should consist of two directories... 

	### quantseq - containing all the zipped files of the quantSeq data
	### rnaseq - containing all the rnaseq data for the condition, mapeed and indexed. 




	################ adding some R packages to a private library (creatied according to the following):https://thecoatlessprofessor.com/programming/working-with-r-on-a-cluster/

	### the R library has been modified to ~/Rlibs
	
	### the following packages are downloaded

	# Rscript -e "install.packages('devtools', '~/Rlibs', 'http://ftp.ussg.iu.edu/CRAN/')"
	#


	
QUANT_ALIGN=$ovalue/polyAmapping_allTimepoints
QUANT_MAP=$ovalue/polyAmapping_allTimepoints/n_100_global_a0/
QUANT_PASPLOTS=$ovalue/PASplot
QUANT_INTERGENIC=$ovalue/intergenicPeaks

mkdir -p "$QUANT_ALIGN"
mkdir -p "$QUANT_MAP"
mkdir -p "$QUANT_PASPLOTS"
mkdir -p "$QUANT_INTERGENIC"
mkdir -p $ovalue/ExtendingINtergenicRegions
mkdir -p $ovalue/coverage

ls "$ivalue"/quantseq/*.{gz,fastq,fq} | perl -pe "s#$ivalue/quantseq/##" > $QUANT_ALIGN/sampleInfo.txt


INPUT="$ivalue"/quantseq/

PARAMETER="$QUANT_ALIGN/sampleInfo.txt"

#set dirs

OUTDIR=$QUANT_ALIGN
LOG=$QUANT_ALIGN/logs/
mkdir -p $LOG




module load cutadapt/1.9.1-foss-2017a-python-2.7.13

 #arrayfile=`ls $INPUT | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
#cutadapt -a $avalue -o "$OUTDIR"/"$arrayfile"_trimmed.fastq  --trim-n $INPUT/"$arrayfile" 2>"$LOG"/stderr_"$index".txt 1>>"$LOG"/stdo_"$index".txt





#init files
touch "$OUTDIR"/logFile.txt ## initialise a log file
touch "$OUTDIR"/readStats.txt ## initialise a file to store your read stats at different processing steps
touch "$OUTDIR"/processingSteps.txt ### just another file to list the processing steps 

touch "$LOG"/stderr_"$index".txt
touch "$LOG"/stdo_"$index".txt


date >"$OUTDIR"/logFile.txt
 date >>"$LOG"/stderr_"$index".txt
 date >>"$LOG"/stdo_"$index".txt

#################
#ADAPTER trimming
#################

module load python/2.7.13-foss-2017a
module load cutadapt/1.9.1-foss-2017a-python-2.7.13
module load fastx-toolkit/0.0.14-foss-2017a
module load fastqc/0.11.5-java-1.8.0_121
module load anaconda2/5.1.0
module load  bedtools/2.25.0-foss-2017a
module load r/3.4.1-foss-2017a-x11-20170314


cutadapt_version=$(cutadapt --version)
echo running cutadapt version "$cutadapt_version" >> "$LOG"/stdo_"$index".txt

## the step below trims the adapter : you should replace this based on your library prep method










	




























########################## Once this is done. I need to find the priming site ##################################



#################### identification of priming sites... 

module load  bedtools/2.25.0-foss-2017a




























 

######################################
###3 getting the sequences +/- 60 nts around the priming sites... 
######################################








#### getting sequences in the 120 nucleotide window




#####################
#### creating nucleotie profile plots...







############# #######
#### addin intergenic peak information 
######################



	if [ -d "$ivalue"/rnaseq/ ]; then
	 
		echo "using the RNAseq samples in the folder "$ivalue"/rnaseq/ "  
	fi


	if [ ! -d "$ivalue"/rnaseq/ ]; then
	
		echo "the rnaseq directory does not exist. Please check if it is named rnaseq."
	fi


	rmd="$PIPELINE/intergenicPeaks/getLongestUTR.R"
	





	## adding intergenic peaks

	module load  bedtools/2.25.0-foss-2017a

	 Rscript --vanilla -e "BIn='$ivalue'; BOut='$ovalue'; ucscDir='$ucscdir'; ensemblDir='$ensembldir'; mode='$RNAseqMode' ; source('$PIPELINE/intergenicPeaks/addingIntergenicPeaks.R')"





######################## final filtering steps

	rmd="$PIPELINE/90PercentFiltering_merging_countingWindows/assignToUTRs.R"





	rmd="$PIPELINE/90PercentFiltering_merging_countingWindows/90PercentFiltering.R"
















######### transferring the data to the input folgrt



#find $ovalue* | egrep ".bed$" | egrep -v "final90percent" | xargs -P 1 gzip

#mkdir /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/$Condition/
#rsync -rva  --exclude "*.bam" --exclude "*.fastq" --exclude "*.fq" $ovalue/* /groups/ameres/Pooja/Projects/zebrafishAnnotation/sequencingRun_december2017/analysis/annotation/stageSpecific/outputs/$Condition/  























