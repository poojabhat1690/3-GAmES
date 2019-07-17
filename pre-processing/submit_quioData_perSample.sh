#!/usr/bin/env bash

while read index; do

	sbatch /groups/ameres/Pooja/Projects/refiningPipeline/pipeline/pre-processing/beforeMapping.new.sh -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -i /scratch/pooja/inputQuio_new/"$index" -o /scratch/pooja/SLAMannotation/dr/quio_jan17_perSample/"$index" -g /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa -t 2 -u /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr11/refSeq_dr11_GRCz11_2019Sep20/ -e /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr11/ensembl_dr11_Ensembl_Genes_93/ -m p -c "$index"


done < /groups/ameres/Pooja/Projects/otherstuff/wt_Ski7annotation_Quio/sampleInfo_EMB.txt



while read index; do

	        sbatch /groups/ameres/Pooja/Projects/refiningPipeline/pipeline/pre-processing/beforeMapping.new.sh -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -i /scratch/pooja/inputQuio_new/"$index" -o /scratch/pooja/SLAMannotation/dr//quio_jan17_perSample/"$index" -g /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa -t 2 -u /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr11/refSeq_dr11_GRCz11_2019Sep20/ -e /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr11/ensembl_dr11_Ensembl_Genes_93/ -m S -c "$index"


	done < /groups/ameres/Pooja/Projects/otherstuff/wt_Ski7annotation_Quio/sampleInfo_Ooc.txt



