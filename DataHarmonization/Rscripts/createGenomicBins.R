setwd("/Users/alexandra/PhD/PyCharmProjects/ALFAssay/DataHarmonization/")
library(Biobase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(QDNAseq)
library(future)
library(GenomicRanges)
library(data.table)
library(stringr)
library(QDNAseq.hg38) 
chr_arm <- read.table("data/chr_arm.csv", sep="\t", header=TRUE)
excl = read.table("data/ENCFF356LFX.bed", sep="\t")
excl = GRanges(excl[,1], IRanges(excl[,2], excl[,3]))
# bins <- createBins(bsgenome=BSgenome.Hsapiens.UCSC.hg38, binSize=5000)

# ### 5Mbp from DELFI
# mappability<- readRDS("data/mappabilityBins-5000k.rds")
# 
# nrow(mappability)
# excl = GRanges(chr_arm)
# mappability= GRanges(mappability)
# mcols(mappability)$arm <- NA
# overlaps <- findOverlaps(query = excl, subject = bins)
# mcols(mappability)$arm[subjectHits(overlaps)] <- mcols(excl)$arm[queryHits(overlaps)]
# mappability <- keepSeqlevels(mappability,c(1:22, "X"),
#               pruning.mode="coarse")
# 
# bindt<- as.data.table(mappability)
# 
# write.table(bindt,"data/5mb_bin_hg38.csv", quote=FALSE, row.names=FALSE,sep=",")

##1M###

# devtools::install_github(asntech/QDNAseq.hg38@main)

bins <- getBinAnnotations(binSize=1000, genome="hg38")@data
    bins = bins[bins$use & bins$blacklist<20 & bins$bases>80 & bins$mappability>70,];
    z = GRanges(bins$chromosome, IRanges(bins$start, bins$end))
    seqlevels(z) = paste0("chr", seqlevels(z));
    mcols(z) = bins[,c("bases", "gc", "mappability", "blacklist", "residual")]
    bins = z;
bins = bins[!overlapsAny(bins, excl)]
mcols(bins)$arm <- NA

arms = GRanges(chr_arm)
seqlevels(arms) = paste0("chr", seqlevels(arms));

overlaps <- findOverlaps(query = arms, subject = bins)
mcols(bins)$arm[subjectHits(overlaps)] <- mcols(arms)$arm[queryHits(overlaps)]
bins <- keepSeqlevels(bins,paste0("chr",1:22),
                             pruning.mode="coarse")

bindt<- as.data.table(bins)
# bindt$seqnames <- str_replace(bindt$seqnames, "chr", "")
bindt <- bindt[!is.na(bindt$arm)]
library(tidyverse)
bindt <- bindt %>% rename(Chromosome = seqnames, Start = start, End = end)

write.table(bindt,"data/1mb_bin_hg38.csv", quote=FALSE, row.names=FALSE,sep="\t")


### 5mb ###

bins <- readRDS("data/mappabilityBins-5000k.rds")
bins = bins[bins$bases>80 & bins$mappability>70,];
z = GRanges(bins$chromosome, IRanges(bins$start, bins$end))
seqlevels(z) = paste0("chr", seqlevels(z));
mcols(z) = bins[,c("bases", "gc", "mappability")]
bins = z;
bins = bins[!overlapsAny(bins, excl)]
mcols(bins)$arm <- NA

arms = GRanges(chr_arm)
seqlevels(arms) = paste0("chr", seqlevels(arms));

overlaps <- findOverlaps(query = arms, subject = bins)
mcols(bins)$arm[subjectHits(overlaps)] <- mcols(arms)$arm[queryHits(overlaps)]
bins <- keepSeqlevels(bins,paste0("chr",1:22),
                      pruning.mode="coarse")

bindt<- as.data.table(bins)
# bindt$seqnames <- str_replace(bindt$seqnames, "chr", "")
bindt <- bindt[!is.na(bindt$arm)]
library(tidyverse)
bindt <- bindt %>% rename(Chromosome = seqnames, Start = start, End = end)

write.table(bindt,"data/5mb_bin_hg38.csv", quote=FALSE, row.names=FALSE,sep="\t")


# str_length("CACAAGCTTCTCCAGCACAGCAACTGTGTCTTATTTCTCCTTGTACTCCCA")


## generate mappability ####

# library(QDNAseq)
# library(Biobase)
# 
# if(!require(BSgenome.Hsapiens.UCSC.hg38)) {
#   source("https://bioconductor.org/biocLite.R")
#   biocLite("BSgenome.Hsapiens.UCSC.hg38")
#   library(BSgenome.Hsapiens.UCSC.hg38)
# }
# 
# bigWig <- snakemake@input[[1]]
# binSize <- as.integer(snakemake@params[["binSize"]])
# 
# print(bigWig)
# print(binSize)
# 
# # Create bin annotations
# bins <- createBins(BSgenome.Hsapiens.UCSC.hg38, binSize)
# 
# # Calculate mappability bin average
# bins$mappability <- calculateMappability(bins, bigWigFile=bigWig, bigWigAverageOverBed='bigWigAverageOverBed')
# 
# saveRDS(bins, paste("rds/mappabilityBins-", binSize, "k", ".rds", sep=""))



