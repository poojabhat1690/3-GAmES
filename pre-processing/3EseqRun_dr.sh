#!/bin/bash





timepoints=(1dpf 256cell 2cell 2dpf 4dpf bud dome oocyte sphere testis)
  
  for i in "${timepoints[@]}"
	   
	   do
		    
		     sbatch /groups/ameres/Pooja/Projects/refiningPipeline/pipeline/pre-processing/beforeMapping.new.sh -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -i /groups/ameres/Pooja/Projects/zebrafishAnnotation/dataFromAndi/for_Pooja/2016_12_QuantSeq_Annotaton/raw_data/"$i" -o /scratch/pooja/SLAMannotation/drQseq/"$i" -g /groups/ameres/bioinformatics/references/danio_rerio/dr11/danRer11.fa -t 2 -u /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr11/refSeq_dr11_GRCz11_2019Sep20/ -e /groups/ameres/Pooja/Projects/zebrafishAnnotation/dr11/ensembl_dr11_Ensembl_Genes_93/ -m p -c "$i"

		      done
