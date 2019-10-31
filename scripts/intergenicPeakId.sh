#!/bin/bash

### intergenic peak identification
awk '{if($14 != "-1") print}'  $QUANT_INTERGENIC/toExtend_longestEnsembl_refSeq_n100_sorted_distances.bed | awk -v OFS='\t' '{print $1, $2, $3, $4,$5, $6,$7,$14}' - >  $QUANT_INTERGENIC//toExtend_longestEnsembl_refSeq_n100_sorted_distances_tmp.bed
awk '{if($6=="+") print}' $QUANT_INTERGENIC/toExtend_longestEnsembl_refSeq_n100_sorted_distances_tmp.bed > $QUANT_INTERGENIC/distances_plus.bed
awk '{if($6=="-") print}' $QUANT_INTERGENIC/toExtend_longestEnsembl_refSeq_n100_sorted_distances_tmp.bed > $QUANT_INTERGENIC/distances_minus.bed

awk -v OFS='\t' '{print $1, $3-200, $3,$4,$5,$6,$7,$8}' $QUANT_INTERGENIC/distances_plus.bed >  $QUANT_INTERGENIC//distances_plus_tmp.bed
awk -v OFS='\t' '{print $1, $2, $2+200,$4,$5,$6,$7,$8}' $QUANT_INTERGENIC/distances_minus.bed > $QUANT_INTERGENIC//distances_minus_tmp.bed

cat $QUANT_INTERGENIC/distances_minus_tmp.bed $QUANT_INTERGENIC//distances_plus_tmp.bed >  $QUANT_INTERGENIC/toExtend_longestEnsembl_refSeq_n100_sorted_distances_tmp.bed

##### getting the signal in the first 200 nts
 singularity exec /groups/ameres/Pooja/dependencies.sif bedtools multicov -split -s -bams  $ivalue/rnaseq/*.bam -bed $QUANT_INTERGENIC//toExtend_longestEnsembl_refSeq_n100_sorted_distances_tmp.bed >  $QUANT_INTERGENIC//customAnnotation_longestTranscripts_100IntoUTR.bed


 #### mean of the reads in the last 200 nucleotides...
  awk -v OFS='\t' '{sum=0; for(i=9; i<=NF; i++) sum += $i; print $1,$2,$3,$4,$5,$6,$7,$8,sum/(NF-8)}' $QUANT_INTERGENIC//customAnnotation_longestTranscripts_100IntoUTR.bed >  $QUANT_INTERGENIC//customAnnotation_longestTranscripts_100IntoUTR_mean.bed

  ### i want to remove the windows that have less that 10 reads in the last 200 nts
  awk '{if($9 > 10) print}' $QUANT_INTERGENIC//customAnnotation_longestTranscripts_100IntoUTR_mean.bed > $QUANT_INTERGENIC//customAnnotation_longestTranscripts_100IntoUTR_mean_tmp.bed



  ###### now i want to create the next window..
  for ((i=1;i<2000;i+=1))
  do

	  	awk -v s=20 -v OFS='\t' '{ if($6 == "+" ) print $1,$2+s,$3+s,$4,$5,$6,$7, $8,$9 ; else print $1,$2-s,$3-s,$4,$5,$6,$7,$8,$9 }' $QUANT_INTERGENIC//customAnnotation_longestTranscripts_100IntoUTR_mean_tmp.bed > $QUANT_INTERGENIC//customAnnotation_longestTranscripts_100IntoUTR_mean_"$i".bed 
			singularity exec /groups/ameres/Pooja/dependencies.sif bedtools multicov -split -s -bams  $ivalue/rnaseq/*.bam -bed $QUANT_INTERGENIC/customAnnotation_longestTranscripts_100IntoUTR_mean_"$i".bed |  awk -v OFS='\t' '{sum=0; for(i=10; i<=NF; i++) sum += $i; print $1,$2,$3,$4,$5,$6,$7,$8,$9,sum/(NF-9), (sum/(NF-9))/$9}' |  awk -v OFS='\t' '{if($11 > 0.1) print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | awk -v OFS='\t' '{if ($i*s<$8) print $0}' - > $QUANT_INTERGENIC/customAnnotation_longestTranscripts_100IntoUTR_mean_final"$i".bed
				cp  $QUANT_INTERGENIC/customAnnotation_longestTranscripts_100IntoUTR_mean_final"$i".bed  $QUANT_INTERGENIC/customAnnotation_longestTranscripts_100IntoUTR_mean_tmp.bed 

			done


			cat $QUANT_INTERGENIC/customAnnotation_longestTranscripts_100IntoUTR_mean_final* > $QUANT_INTERGENIC/totalWindows_extensions.bed

			grep noPAS  $QUANT_MAP/nonOverlapping_total.bed | awk '{if ($10<0.24) print}' > $QUANT_MAP/nonOverlapping_total_noPAS.bed
			grep -v noPAS $QUANT_MAP/nonOverlapping_total.bed | awk '{if ($10<0.36) print}' > $QUANT_MAP/nonOverlapping_total_PAS.bed

			cat $QUANT_MAP/nonOverlapping_total_noPAS.bed $QUANT_MAP/nonOverlapping_total_PAS.bed > $QUANT_MAP/nonOverlapping_total_passing.bed


			singularity exec /groups/ameres/Pooja/dependencies.sif bedtools intersect -wa -wb -b $QUANT_INTERGENIC/totalWindows_extensions.bed -a $QUANT_MAP/nonOverlapping_total_passing.bed > $ovalue/ExtendingINtergenicRegions/overlapping_windows.bed
			sort -u -k4,4 $ovalue/ExtendingINtergenicRegions/overlapping_windows.bed  > $ovalue/ExtendingINtergenicRegions/allIntergenicPeaks_n100_new.txt



