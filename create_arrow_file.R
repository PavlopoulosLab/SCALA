library("optparse")
'%!in%' <- function(x,y)!('%in%'(x,y))

# Example: Rscript create_arrow_file.R -s mm10 -a example -f fragments.ucsc.tsv.gz -t 4 -m 1000 -p 1 
# Will output example.arrow
 
option_list = list(
  make_option(c("-s", "--species"), type="character", default=NULL, 
              help="Species genome build. Examples: hg19, mm10", metavar="character"),
  make_option(c("-a", "--arrow_prefix"), type="character", default="myProject", 
              help="Arrow file prefix. Default name will be set as myProject", metavar="character"),
  make_option(c("-f", "--fragmentfile"), type="character", default="fragments.tsv.gz", 
              help="Fragments file path. Default file path will be set fragments.tsv.gz", metavar="character"),
  make_option(c("-t", "--minTSS"), type="integer", default=4, 
              help="The minimum numeric transcription start site (TSS) enrichment
          score required for a cell to pass filtering for use in
          downstream analyses. Default value will be set to 4.", metavar="integer"),
  make_option(c("-m", "--minFrags"), type="integer", default=1000, 
              help="The minimum number of mapped ATAC-seq fragments required per
          cell to pass filtering for use in downstream analyses. Default value will be set to 1000.", metavar="integer"),
  make_option(c("-p", "--threads"), type="integer", default=8, 
              help="The number of threads to be used for parallel computing. Default value will be set to 8.", metavar="integer")
)
 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
if (is.null(opt$species) || opt$species %!in% c("mm9","mm10","hg19","hg38")){
  stop("Please provide a valid UCSC genome build file. Valid options include only UCSC human and mouse genome builds hg19, hg38, and mm9, mm10 respectively.")
}
if (!file.exists(opt$fragmentfile)){
  stop(sprintf("Fragment file %s does not exist, please provide a valid file path.",opt$fragmentfile))
}
if(.Platform$OS.type == "windows"){
  opt$threads <- 1
}

library(ArchR)

addArchRGenome(opt$species)
addArchRThreads(threads = opt$threads) 
inputFiles<-opt$fragmentfile
names(inputFiles)<-opt$arrow_prefix

print("Creating arrow file, please wait ...")
ArrowFiles <- createArrowFiles(
	inputFiles = inputFiles,
	sampleNames = names(inputFiles),
	minTSS = opt$minTSS,
	minFrags = opt$minFrags, 
	threads = opt$threads, 
	force = TRUE
)
print("Creation of arrow file completed!!")
