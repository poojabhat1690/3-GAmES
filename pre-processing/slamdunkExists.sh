#!/bin/bash


#### this pipeline used SLAMdunk for mapping quantSeq datasets - SLAMdunk is available as a singulatity module on docker.   
######## checking if the module exists in the folder, if not pulling it from docker. 
module load singularity/2.5.2

PIPELINE=/groups/ameres/Pooja/Projects/refiningPipeline/pipeline/

cd "$PIPELINE"

FILE="$PIPELINE"/slamdunk-v0.3.4.simg

if [ -f "$FILE" ]; then
	echo "$FILE exist"
		  else 
	echo "$FILE does not exist - pulling it from docker"
	singularity pull docker://tobneu/slamdunk:v0.3.4
		
fi
							                                                                                                                                                                                                                                                                                                                                   


