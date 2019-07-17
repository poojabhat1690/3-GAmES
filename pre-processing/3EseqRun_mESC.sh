#!/bin/bash

### running 3Eseq for the mESC data


sbatch /groups/ameres/Pooja/Projects/refiningPipeline/pipeline/pre-processing/beforeMapping.new.sh -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -i /scratch/pooja/quantSeq_mESC/SR50/ -o /scratch/pooja/SLAMannotation/mESC/SR50/ -g /groups/ameres/bioinformatics/references/mmu/mm10/chr/mmu_wholegenome.fa -t 10 -u /groups/ameres/Pooja/Projects/refiningPipeline/pipeline/mm10/refSeq_mm10_GRCm38_06-march-2017/ -e /groups/ameres/Pooja/Projects/refiningPipeline/pipeline/mm10/ensembl_mm10_Ensembl_Genes_87_06-march-2017/ -m S -c SR50_combined 


sbatch /groups/ameres/Pooja/Projects/refiningPipeline/pipeline/pre-processing/beforeMapping.new.sh -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG -i /scratch/pooja/quantSeq_mESC/SR100/ -o /scratch/pooja/SLAMannotation/mESC/SR100/ -g /groups/ameres/bioinformatics/references/mmu/mm10/chr/mmu_wholegenome.fa -t 10 -u /groups/ameres/Pooja/Projects/refiningPipeline/pipeline/mm10/refSeq_mm10_GRCm38_06-march-2017/ -e /groups/ameres/Pooja/Projects/refiningPipeline/pipeline/mm10/ensembl_mm10_Ensembl_Genes_87_06-march-2017/ -m S -c SR100_combined





