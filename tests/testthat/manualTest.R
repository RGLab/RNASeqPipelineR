#' Testing for RNASeqPipeline.R
#'
#' Manual testing due to time considerations.
#'
library(RNASeqPipelineR)
library(data.table)
library(Biobase)

##source("/home/cmurie/cmurie_working/git/RNASeqPipelineR/R/RNASeqPipelineR.R")

## local instance of reference genome
utils_dir <- "/shared/silo_researcher/Gottardo_R/10_ref_files/Reference_Genome/Homo_sapiens/UCSC/hg38/"

## get temp directory for project creation
tmp <- tempdir()

createProject(project_name = "test", path=tmp, load_from_immport = FALSE)
loadProject(project_dir=tmp, name="test")

## load test fastq files stored in extdata
## contains test1_R1.fastq, test1_R2.fastq, test2_R1.fastq, test2_R2.fastq,
fpath <- system.file("extdata", "testData.zip", package="RNASeqPipelineR")
unzip(zipfile=fpath, exdir=paste0(tmp, "/test/FASTQ/"))
fpath2 <-  system.file("extdata", "testPheno.csv", package="RNASeqPipelineR")
file.copy(from=fpath2, to=paste0(tmp, "/test/RAW_ANNOTATIONS/anno.csv"))

## reference should already be built. test this function separately
buildReference(path=utils_dir, gtf_file="UCSCDec2016.gtf", fasta_file="hg38.fa", isoformsFile="UCSCKnownIsoformsDec2016.txt", doSTAR=TRUE, threads=6, name="hg38")

AlignmentSTAR(parallel_threads=2, star_threads=2, paired=TRUE, paired_pattern=c("_R1.fastq", "_R2.fastq"))
RSEMCalculateExpression(parallel_threads=1,bowtie_threads=1,paired=TRUE, nchunks=1,slurm=FALSE, fromBAM=TRUE, fromSTAR=TRUE)
RSEMAssembleExpressionMatrix(force=TRUE)
##BioCAnnotate("TxDb.Hsapiens.UCSC.hg38.knownGene",force = TRUE,lib = 'org.Hs.eg.db', na.rm=TRUE)
annotateUCSC(genome="hg38", force=TRUE)  ## replaces BioCAnnotate
runFastQC(ncores=2)
qc_matrix <- QualityControl(paired=TRUE)

mData <- mergeData()

## multiple cluster ids can map to same gene symbol so sum them and return ExpressionSet
count_eset <- sumDuplicates(mData$counts, mData$featureData, mData$annoData)
tpm_eset <- sumDuplicates(mData$tpms, mData$featureData, mData$annoData)

save(count_eset, tpm_eset, file=paste0(tmp, "/test/OUTPUT/test.RData"))


## remove temp directory and all contents (just to be safe)
unlink(tmp, recursive=TRUE, force=TRUE)

testSummation(count_eset, tmp)


## manually check that the read values match between the generated read
## count table and the original FASTQ file.
testSummation <- function(counts=count_eset, prefix=tmp) {
  
  PREFIX <- prefix 
  
  dat <- exprs(count_eset)
  
  rsem <- read.table(paste0(PREFIX, "/test/RSEM/test1Aligned.toTranscriptome.out.genes.results"), sep="\t", header=TRUE, row.names=NULL)
  counts <- rsem[,c("gene_id", "expected_count")]
  
  annoFile <- RNASeqPipelineR:::ucschg38Table
  syms <- as.character(unique(annoFile$hg38.kgXref.geneSymbol))
  
  tol = .Machine$double.eps ^ 0.5
  
  for(g in syms) {
    
    ids <- unique(annoFile[annoFile$hg38.kgXref.geneSymbol==g,]$hg38.knownIsoforms.clusterId)
    read <- dat[rownames(dat)==g,1]
    
    rDat <- 0
    for(id in ids) {
       rDat <- rDat + counts[counts$gene_id==id,]$expected_count
    }
    
    if(rDat - read > tol) {
      cat(g, ": ",rDat, " ", read, "\n")
    }
    
    rDat - read  > tol
    
  } ## end for g
  
} ## end testSummation

