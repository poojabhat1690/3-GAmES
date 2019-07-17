if(require("checkmate")){
  print("checkmate is loaded correctly")
} else {
  print("trying to install checkmate")
  install.packages("checkmate")
  if(require(checkmate)){
    print("checkmate installed and loaded")
  } else {
    stop("could not install checkmate")
  }
}


if(require("Biostrings")){
  print("Biostrings is loaded correctly")
} else {
  print("trying to install Biostrings")
  BiocManager::install("Biostrings",ask=FALSE )
  if(require(Biostrings)){
    print("Biostrings installed and loaded")
  } else {
    stop("could not install Biostrings")
  }
}


if(require("ggplot2")){
  print("ggplot2 is loaded correctly")
} else {
  print("trying to install ggplot2")
  install.packages("ggplot2")
  if(require(ggplot2)){
    print("ggplot2 installed and loaded")
  } else {
    stop("could not install ggplot2")
  }
}

if(require("reshape")){
  print("reshape is loaded correctly")
} else {
  print("trying to install reshape")
  install.packages("reshape")
  if(require(reshape)){
    print("reshape installed and loaded")
  } else {
    stop("could not install reshape")
  }
}


if(require("biomaRt")){
  print("biomaRt is loaded correctly")
} else {
  print("trying to install biomaRt")
  BiocManager::install("biomaRt",ask = F)
  if(require(biomaRt)){
    print("biomaRt installed and loaded")
  } else {
    stop("could not install biomaRt")
  }
}


## this is used for making nucleotide profiles to differnetiate between false and true positive sites.

library(checkmate)
library(Biostrings)
library(ggplot2)
library(reshape)
#library(GenomicFeatures)
library(biomaRt)



#reading in the fasta and the bed file for all ends

fasta_polyApeaks_120bps = read.table(paste0(InPath, "/peaks_",threshold,"_120bps.fa"),stringsAsFactors = F)
peaks_total_modified = read.table(paste0(InPath, "/peaks_",threshold,"_120bps.bed"))


assertDataFrame(fasta_polyApeaks_120bps)
nRow_fasta = nrow(fasta_polyApeaks_120bps)/2
assertDataFrame(peaks_total_modified,nrows = nRow_fasta)

## now adding the fasta sequences to the bed file 

sequences_polyApeaks_120bps = fasta_polyApeaks_120bps[seq(2,nrow(fasta_polyApeaks_120bps),2),]
## rearranging the bed file to put the exact coordinates of the priming sites (not +/- 60nts)

peaks_total_modified$V2 = peaks_total_modified$V7
peaks_total_modified$V3 = peaks_total_modified$V8
peaks_total_modified = peaks_total_modified[,c(1:6)]

peaks_bedFile = cbind(peaks_total_modified,sequences_polyApeaks_120bps)
peaks_bedFile$peakName = paste(peaks_bedFile$V1,peaks_bedFile$V2,sep="_")
peaks_bedFile$downstreamSeq = substr(x = peaks_bedFile$sequences_polyApeaks_120bps,start = 61,stop = 80)
content_nucleotides = mapply(function(x) alphabetFrequency(DNAString(x)),peaks_bedFile$downstreamSeq)
peaks_bedFile$totalAs = content_nucleotides["A",]/20
write.table(peaks_bedFile,paste0(InPath, "/sequences_120nts.bed"),sep="\t",quote = F,row.names = F,col.names = F)





