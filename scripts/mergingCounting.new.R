### creating counting windows for counting T>C conversion and for transcriptional output
library(checkmate)
library(dplyr)
#### counting windows for half life calculations
#BOut="//clustertmp/bioinfo/pooja/SLAMannotation/dr/output/"
annotation_custom = read.table(paste0(BOut, "/final90percent/ends_greater90percent_intergenic_n100.txt"),stringsAsFactors = F,sep="\t",header = T)

assertDataFrame(annotation_custom,ncols=11)

annotation_custom = annotation_custom[,c(1:8)]

colnames(annotation_custom) = paste("V",c(1:8),sep="")



allEnsembl = read.table(paste0(ensemblDir,"/transcriptStartsAndEnds_all.txt"),sep="\t",header=F,stringsAsFactors = F)


classesToInclude = c("antisense", "bidrectional_promoter_lncRNA", "lincRNA", "macro_lncRNA", "processed_transcript", "sense_intronic", "sense_overlapping","protein_coding")
allEnsembl = allEnsembl[allEnsembl$V8 %in% classesToInclude,]

### not including these anymore but just including long non-coding RNAs from ensembl






allEnsembl <- unique(allEnsembl[,c("V1", "V2", "V3", "V4", "V6", "V7")])
colnames(allEnsembl) <- c("chr", "start", "end", "gid", "strand", "tid")

allEnsembl <- cbind(allEnsembl[,1:4], 0,  allEnsembl[,5:6])

allEnsembl$V8 = "EnsemblOriginal"
colnames(allEnsembl) <- paste0("V", 1:8)



annotation_custom_positive = annotation_custom[which(annotation_custom$V6 == "+"),]
annotation_custom_negative = annotation_custom[which(annotation_custom$V6 == "-"),]




allEnsembl_positive = allEnsembl[which(allEnsembl$V6 == "+"),]
allEnsembl_negative = allEnsembl[which(allEnsembl$V6 == "-"),]











######################## positive strand ###########################
library(GenomicRanges)
annotation_custom_positive = rbind(annotation_custom_positive,allEnsembl_positive)

annotation_custom_positive$V2 = annotation_custom_positive$V3 -250
annotation_custom_positive$V2 = annotation_custom_positive$V2 + 1

annotation_custom_positive_split = split(annotation_custom_positive,f = annotation_custom_positive$V4,drop = T )

positive_ranges = lapply(annotation_custom_positive_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4)))




allAnnotations_plus_ranges_reduced = lapply(positive_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )



reducedToDf = function(reduced){
  reduced <- data.frame(seqnames=seqnames(reduced),
                        starts=start(reduced),
                        ends=end(reduced),
                        names=c(names(reduced)),
                        scores=0,strand = strand(reduced))
  return(reduced)
}

allAnnotations_plus_ranges_reduced_df = lapply(allAnnotations_plus_ranges_reduced,function(x) reducedToDf(x))

################## minus strand ############################

annotation_custom_negative = rbind(annotation_custom_negative,allEnsembl_negative)

annotation_custom_negative$V3 = annotation_custom_negative$V2 + 250

## changing to 1 based

annotation_custom_negative$V2 = annotation_custom_negative$V2 + 1

annotation_custom_negative_split = split(annotation_custom_negative,f = annotation_custom_negative$V4,drop = T )

negative_ranges = lapply(annotation_custom_negative_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4)))


allAnnotations_minus_ranges_reduced = lapply(negative_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )




allAnnotations_minus_ranges_reduced_df = lapply(allAnnotations_minus_ranges_reduced,function(x) reducedToDf(x))


allAnnotations_plus_ranges_reduced_df = do.call(rbind,allAnnotations_plus_ranges_reduced_df)
allAnnotations_minus_ranges_reduced_df = do.call(rbind,allAnnotations_minus_ranges_reduced_df)

allAnnotations = rbind(allAnnotations_plus_ranges_reduced_df,allAnnotations_minus_ranges_reduced_df)

### converting back to 0 based annotations : 

allAnnotations$starts = allAnnotations$starts -1

#valid chromosomes

write.table(allAnnotations,paste0(BOut, "/final90percent/allAnnotations.bed"),sep="\t",quote = F,row.names = F,col.names = F)

#########################################################################
## counting windows for transcriptional output and multimapping - for this we need the ensembl 3' utr annotations, refSeq 3' utr annotations
## and intergenic peaks that pass the 90% threshold
#########################################################################

refSeq = read.table(paste0(ucscDir, "//refSeq_mrna_utrsPresent.bed"),stringsAsFactors = F)
ensembl = read.delim(paste0(ensemblDir, "/proteinCoding_annotatedUTRs.bed"),stringsAsFactors = F,header = F)

allAnnotations <- cbind(allAnnotations, "250CountWindow")
colnames(allAnnotations) <- paste0("V", 1:7)
refSeq_ensembl = rbind(refSeq,ensembl,allAnnotations)
refSeq_ensembl_positive = refSeq_ensembl %>% filter(V6=="+")
refSeq_ensembl_negative = refSeq_ensembl %>% filter(V6=="-")









total_positive = refSeq_ensembl_positive

total_positive$V2 = total_positive$V2 +1
total_negative = refSeq_ensembl_negative

total_negative$V2 = total_negative$V2 + 1




### positive strand 

total_positive_split = split(total_positive,f = total_positive$V4,drop = T )

total_positive_split_ranges = lapply(total_positive_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4))
)


total_positive_reduced = lapply(total_positive_split_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )




total_positive_reduced_df = lapply(total_positive_reduced,function(x) reducedToDf(x))

total_positive_reduced_df = do.call(rbind,total_positive_reduced_df)


####### negative strand 



total_negative_split = split(total_negative,f = total_negative$V4,drop = T )

total_negative_split_ranges = lapply(total_negative_split,function(x) with(x,GRanges(V1, IRanges(start = V2,end = V3),strand = V6,score=0,names=V4))
)


total_negative_reduced = lapply(total_negative_split_ranges,function(x) unlist(reduce(split(x, elementMetadata(x)$names))) )




total_negative_reduced_df = lapply(total_negative_reduced,function(x) reducedToDf(x))

total_negative_reduced_df = do.call(rbind,total_negative_reduced_df)

countingWindowsTranscriptionalOutput = rbind(total_positive_reduced_df,total_negative_reduced_df)
countingWindowsTranscriptionalOutput$starts = countingWindowsTranscriptionalOutput$starts -1


write.table(countingWindowsTranscriptionalOutput,paste0(BOut, "/final90percent/countingWindows_transcriptionalOutput.bed"),sep="\t",quote = F,row.names = F,col.names = F)



