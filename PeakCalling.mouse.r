################################################################################
                                        # Peak Calling
################################################################################
library(CLIPanalyze)
library(GenomicRanges)
library(SummarizedExperiment)
library(rtracklayer)

#bandwidthParam <- as.integer(args[2])  ###50 75 100
#thresholdParam <- as.integer(args[3])  ###100 150 200
bandwidthParam <- 50
thresholdParam <- 10

## BAM file name
dir_path <- file.path("/","data","mayrc","zheng","Mathieu","Bam")
files <- c(file.path(dir_path,'C57.filtered.bam'),file.path(dir_path,'C57_Cre.filtered.bam'),file.path(dir_path,'Mettle3_Cre.filtered.bam'))

allChrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY","chrM")

### Peak calling
peak.data <- findPeaks(bamfiles = files, sample.names = "Mathieu_all", chrs = allChrs,
                       stranded = TRUE, paired.end = FALSE, diffanalysis = FALSE,
                       annot.order=c("utr3", "exon", "utr5", "intron", "utr5*", "utr3*"),
                       utr5.extend = 5000, utr3.extend = 5000,
                       genome = "mm10", counting = FALSE, count.exons.only = FALSE,
                       bandwidth = bandwidthParam, count.threshold = thresholdParam,
                       nthreads = 6)

### save results
peakFile <- sprintf("/data/mayrc/zheng/Mathieu/PeakCalling/RDS/All.rds")
saveRDS(object = peak.data, file = peakFile)
export.bed(peak.data$peaks, sprintf("./All.3UTRPeaks.new.bed"))

