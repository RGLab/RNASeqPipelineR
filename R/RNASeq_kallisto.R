#' kallisto methods
#' 
#' Implement kallisto alignment and quantification 
#' 
#' Methods to apply kallisto alignment and quantification to RNASeq data.
#' 


#' Build transcript reference index for kallisto
#' 
#' Build kallisto transcript reference if it is not already built. To force the build then
#' set force to TRUE.
#'
#' @param path path to directory containing transcript fasta file
#' @param fasta_file name of fasta file to be used
#' @param name name of index to be built
#' @param force if TRUE build index if it already exists else don't.
#' @export
kallistoBuildTranscriptIndex <- function(path=NULL,
                                         fasta_file="Homo_sapiens.GRCh38.cdna.all.fa",
                                         name="GRCH38_R91.idx",
                                         force=FALSE) {
  
  ## error checking
  if (is.null(fasta_file)){
    stop("'fasta_file' is NULL. You must provide a fasta_file")
  }
  
  if (is.null(name)){
    stop("'name' is NULL. You must provide a name for index")
  }
  
  if(!file.exists(path)) {
    stop(paste("'path'", path, "does not exist"))
  }
  
  if(!file.exists(file.path(path, fasta_file))) {
    stop(paste("'fasta_file'", fasta_file, "does not exist"))
  }
  
  if (is.null(path) & !file.exists(getDir("Utils"))) {
    stop("Please specify where to build the reference genome")
  }
  
  ## get and assign path to kallisto index directory
  if(!is.null(path)) {
    subdirs <- getConfig()[["subdirs"]]
    subdirs[["Utils"]] <- path
    assignConfig("subdirs", subdirs)
  }
  path <- getDir("Utils")
  
  fasta_file_path <- file.path(path, fasta_file)
  index_file_path <- file.path(path, name)
  
  ## exit function if index is already built and not forced or build index
  if(!force & file.exists(index_file_path)) {
    message("index is already built. Set 'force' to TRUE to rebuild index")
  } else {
    command = paste("kallisto index -i", index_file_path, fasta_file_path)
    system(command)
  } ## end if !force
  
  assignConfig("reference_genome_name", name)
  
} ## end kallistoBuildTranscriptIndex



#' Use the Kallisto tool to align and quantify RNASeq FASTQ files.
#'
#' Uses Kallisto to pseudo-align reads in FASTQ files against the reference genome 
#' transcripts. Optionally you can specify paired end reads. The code assumes
#' paired reads have different suffixes defined by 'paired_pattern' (i.e. 
#' sampleA_read_R1.fastq, sampleA_read_R2.fastq). Read 1 is assumed to be 
#' upstream and read 2 is assumed to be downstream. FASTQ files must have the 
#' suffix 'fastq' or fq'.
#' If running locally the number of available nodes should be at least the number of
#' kallisto cores plus one times the number of parallel threads.
#' 
#' @param kallisto_threads specify how many threads kallisto 
#' should use. This will be input as the '-t' parameter in the kallisto command.
#' Default is 3 as the majority of gizmo nodes have 4 processors.
#' @param paired specify whether you have paired reads or not.
#' @param minutes_requested number of minutes requested for the job 
#' (when submitting a slurm job). Generally kallisto takes less than 5 minutes to 
#' complete. Scicomp recommended 10m. Ignored if slurm = FALSE
#' @param slurm if \code{TRUE} job is submitted as a slurm batch job, otherwise it's run 
#' on the local machine. Slurm jobs will honour the nchunks and minutes_requested arguments.
#' @param parallel_threads Number of threads to run concurrent kallisto jobs if not using 
#' SLURM. Total number of threads required are parallel_threads*(kallisto_threads+1).
#' @param slurm_partition the slurm partition to submit to. Defaults
#' to campus.Ignored if slurm=FALSE
#' @param force If TRUE then align all FASTQ files else only align those
#'  files that have not been aligned. If force is FALSE then only fastq files without
#'  a directory in the Kallisto directory will be aligned.
#' @param paired_pattern specify the suffix of the paired-end fastq 
#' file names. Each suffix must be of the same length.
#' @param fastqPath specify path to FASTQ files if different than default (FASTQ directory)
#' @param mail email address to send failure message to, if desired.
#' @export
kallistoAlign <- function(kallisto_threads=3, paired=TRUE,
                          minutes_requested=10, slurm=TRUE,
                          parallel_threads=1,
                          slurm_partition=NULL, #slurm_partition="gottardo_r",
                          force=FALSE,
                          paired_pattern=c("_1.fastq", "_2.fastq"),
                          fastqPath="", mail=NULL){

  ## make sure we have enough threads if running locally
  if((kallisto_threads+1)*parallel_threads > parallel::detectCores() &! slurm){
    stop("The number of kallisto threads is more than the number of cores detected by detectCores() on the local machine for non-slurm jobs")
  }

  ## check that paired suffixes are of same length
  if(nchar(paired_pattern[1]) != nchar(paired_pattern[2])) {
    stop("The paired suffixes are of different lengths")
  }

  ## read fastq files from user inputted directory or use default
  fastq_dir <- ifelse(fastqPath=="",  getDir("FASTQ"),  fastqPath)
  if(!file.exists(fastq_dir)) {
    stop(paste(fastq_dir, "doesn't exist\n"))
  }

  ## create necessary file paths
  kallisto_dir <-  getDir("KALLISTO")
  ref_full_path = file.path(getDir("Utils"), getConfig()[["reference_genome_name"]])

  ## get names of fastq files
  orig <- list.files(fastq_dir, pattern="fastq|fq")
  sufLength <- nchar(paired_pattern[1])
  sampleNames <- unique(substr(orig, start=1, stop=nchar(orig) - (sufLength)))

  ## parse Kallisto output to get names of fastq files which have already been aligned.
  done <- list.dirs(path=kallisto_dir, full.names=FALSE, recursive=FALSE)
  remove <- c("error", "log")
  done <- done[!done %in% remove]  ## remove logging directories

  ## obtain names of fastq files that have not been aligned.
  keep <- sampleNames
  if(!force) {
    keep <- setdiff(sampleNames, done)
    if(length(keep)==0){
      message("All FASTQ files already aligned. Use 'force=TRUE' to force realignment")
      return()
    }
  } ## end if force

  ## Write slurm preamble to a shell script and execute for all fastq files
  if(slurm){
    
      ## set up email section if email address is given
      mail_string <- mail_address <- ""
      if(!is.null(mail)){
          mail_string <- "#SBATCH --mail-type=FAIL"
          mail_address <- paste0("#SBATCH --mail-user=", mail)
      }
    
      con <- file(file.path(getDir("FASTQ"), "batchSubmitJob.sh"), open = "w")
      writeLines(c("#!/bin/bash",
                   "sbatch <<EOT",
                   "#!/bin/bash",
                   "#SBATCH --job-name kallisto_$1",
                   paste0("#SBATCH -t 0-00:", minutes_requested,
                          " # Runtime in D-HH:MM"),
                   ifelse(is.null(slurm_partition), "#SBATCH -p campus",
                          paste0("#SBATCH -p ",slurm_partition," # Partition to submit to")),
                   "#SBATCH --nodes=1",
                   "#SBATCH --ntasks=1",
                   "#SBATCH --tasks-per-node=1",
                   paste0("#SBATCH --cpus-per-task=", kallisto_threads+1),
                   paste0("#SBATCH -o ", file.path(kallisto_dir, "log/kallisto_$1.out")," # File to which STDOUT will be written"),
                   paste0("#SBATCH -e ", file.path(kallisto_dir, "error/kallisto_$1.err"), " # File to which STDERR will be written"),
                   mail_string,
                   mail_address),
                 con=con)

      if(!paired){
        stop("single end read alignment not implemented yet.")
      }else{
        run_command <-  paste0('parallel -j ', ncores,
                               ' fastqc {} -o "',getConfig()[['subdirs']][['FASTQC']],
                               '" -q ::: "', paste(notrun$file, collapse='" "'), '"')
        command<- paste0("cd ", fastq_dir, " && kallisto quant -i ", ref_full_path,
                         " -t ", kallisto_threads, " -o ", paste0(kallisto_dir, "/$1"),
                         " $1", paired_pattern[1],
                         " $1", paired_pattern[2])
        writeLines(command, con=con)
      } ## !paired
     
      writeLines("EOT", con=con)
      close(con)

    ## run kallisto on a single set of fastq files
    runKallistoSLURM <- function(fastqName) {
      kalliCommand <- paste0("bash ", file.path(fastq_dir,"batchSubmitJob.sh "),
                             fastqName)
      system(kalliCommand)
    } ## end runKallisto

    ret <- apply(matrix(keep, ncol=1), 1, runKallistoSLURM)

  } else { 
      ## run kallisto for all fastq files on workstation
      runKallisto <- function(fastqName, f_dir=fastq_dir, r_path=ref_full_path,
                              k_threads=kallisto_threads, k_dir=kallisto_dir,
                              p_pattern=paired_pattern) {
  
          kalliCommand <- paste0("cd ", f_dir, " && kallisto quant -i ", r_path,
                                 " -t ", k_threads, " -o ",
                                 paste0(k_dir, "/", fastqName, " "),
                                 fastqName, p_pattern[1], " ",
                                 fastqName, p_pattern[2], " 2>&1 | tee ", k_dir,
                                 "/error/kallisto_", fastqName, ".err")
          system(kalliCommand)
      } ## end runKallisto
  
      ## run kallisto on all fastq files
      nothing <- mclapply(keep, runKallisto, mc.cores=parallel_threads)
  
  } ## end if else

} ## end kallistoAlign


#' Assemble kallisto quality control information
#' 
#' Read kallisto error files and extract and write read count, aligned read count, and 
#' estimated average length of reads. Write 'kallisto_QC.csv' to OUTPUT directory.
#' 
#' @param paired If true then paired data else single read data. Default is
#' TRUE.
#' @export
kallistoReadQC <- function(paired=TRUE) {
  
  ## get list of kallisto output directories
  kallisto_dir <-  getDir("KALLISTO")
  kalRes <- list.dirs(path=kallisto_dir, full.names=FALSE, recursive=FALSE)
  remove <- c("error", "log")
  kalRes <- kalRes[!kalRes %in% remove]  ## remove logging directories
  errDir <- paste0(kallisto_dir, "/error")
  
  ## read one file in error directory and extract #reads, #aligned reads, and
  ## estimated read length
  readErrDir <- function(name, dir=kallisto_dir, isPaired=paired) {
    errFile <- paste(dir, name, "run_info.json", sep="/")
    filey <- readLines(errFile)
    if(isPaired) {
      lines <- str_split(filey, " ")
      outie <- c(as.numeric(str_match(lines[[4]][2], "[0-9]+")),  ## num reads
                 as.numeric(str_match(lines[[5]][2], "[0-9]+")),  ## num aligned
                 as.numeric(str_match(lines[[7]][2], "[0-9]+")))  ## percent aligned
    } else {
      stop("Reading single end RNASeq data not implemented")
    }
    return(outie)
  } ## end readErrDir
  
  
  ## read one file in error directory and extract #reads, #aligned reads, and
  ## estimated read length
  readDir <- function(name, dir=errDir, isPaired=paired) {
    errFile <- paste0(dir, "/kallisto_", name, ".err")
    filey <- readLines(errFile)
    if(isPaired) {
      lines <- str_split(filey, " ")
      outie <- c(lines[[11]][c(3,5)])                 ## num reads and aligned reads
      outie <- c(outie, lines[[12]][6])               ## avg estimated read length
      outie <- as.numeric(gsub(",", "", outie))       ## remove commas
      outie <- c(outie, signif(outie[2]/outie[1], 3)) ## percent aligned reads
    } else {
      stop("Reading single end RNASeq data not implemented")
    }
    return(outie)
  } ## end readDir
  
  
  qcMat <- t(apply(matrix(kalRes, ncol=1), 1, readDir))
  qcMat <- cbind(kalRes, qcMat)
  dimnames(qcMat) <- list(kalRes, c("Sample", "numReads", "numAlignedReads", 
                                    "estAvgReadLength", "PercentAligned"))
  write.csv(qcMat, paste0(kallisto_dir, "/kallisto_QC.csv"), quote=FALSE)
  
} ## end kallistoReadQC


#' Assemble kallisto output into a single table.
#' 
#' Combine all count and tpm results from kallisto output into a single table. Also 
#' read error files and extract read count, aligned read count, and estimated average
#' length of reads. Write 'kallisto_counts.csv' and 'kallisto_tpm.csv' to OUTPUT directory.
#' 
#' @param paired If true then paired data else single read data. Default is
#' TRUE.
#' @export
kallistoAssembleOutput <- function(paired=TRUE) {
  
  ## get list of kallisto output directories
  kallisto_dir <-  getDir("KALLISTO")
  kalRes <- list.dirs(path=kallisto_dir, full.names=FALSE, recursive=FALSE)
  remove <- c("error", "log")
  kalRes <- kalRes[!kalRes %in% remove]  ## remove logging directories
  
  ## get and write kallisto QC information
  kallistoReadQC(paired=paired)
  
  ## get ensembl gene ids
  geneIDs <- fread(file.path(kallisto_dir, kalRes[1], "abundance.tsv"),
                   sep="\t", header=TRUE, select="target_id")
  
  ## remove version suffix from geneIDs
  geneIDs[, target_id:=tstrsplit(target_id, ".", fixed=TRUE, keep=1)]
  
  ## read in one kallisto output file
  readFiles <- function(file, path=kallisto_dir) {
    fullPath <- file.path(path, file, "abundance.tsv")
    tmp <- fread(fullPath, header=TRUE, sep="\t", select=c("est_counts", "tpm"))
  }
  
  ## read in all files and extract counts and tpm
  ret <- lapply(kalRes, readFiles)
  counts <- data.table(sapply(ret, "[[", 1))
  tpms <- data.table(sapply(ret, "[[", 2))
  colnames(counts) <- colnames(tpms) <- kalRes
  
  counts[, EnsemblIDs:=geneIDs$target_id]
  tpms[,EnsemblIDs:=geneIDs$target_id]
  
  ## map gene symbols to ensembl IDs
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  genes <- data.table(getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id','hgnc_symbol'),
                            mart = ensembl))
  
  countTranscripts <- merge(counts, genes, by.x="EnsemblIDs", by.y="ensembl_transcript_id", all.x=TRUE)
  tpmTranscripts<-  merge(counts, genes, by.x="EnsemblIDs", by.y="ensembl_transcript_id", all.x=FALSE)
  
  ## get list of ensembl_transcript_ids that do not exist in database.
  missingIDs <- countTranscripts[is.na(ensembl_gene_id),]
  countTranscripts <- countTranscripts[!is.na(ensembl_gene_id),]
  
  ## aggregate on gene symbol
  geneCounts <- countTranscripts[,lapply(.SD, sum()), by=hgnc_symbol, 
                                 .SDcols=kalRes]
  geneTPMs <- tpmTranscripts[,lapply(.SD, sum()), by=hgnc_symbol, 
                             .SDcols=kalRes]
  
  ## write gene counts and tpms
  outPath <- getDir("OUTPUT")
  write.csv(geneCounts, paste0(outPath, "/kallisto_counts.csv"), row.names=FALSE)
  write.csv(geneTPMs, paste0(outPath, "/kallisto_tpms.csv"), row.names=FALSE)
  
} ## end kallistoAssembleOutput


#' Collect kallisto and fastqc information.
#' 
#' Assemble kallisto and fastqc QC information and potentially combine with 
#' annotation data. Write 'kallisto_Pheno.csv' to OUTPUT directory.
#' 
#' @param paired if TRUE then paired end data else single read data
#' @param doAnnotation if TRUE then add annotation data to qc info. The 
#' annotation file must be the sole .csv file in the RAW_ANNOTATIONS directory.
#' @param phenoMatchLabel column name that contains the fastq file names without
#' the paired suffix label in the annotation data. Used to match information between 
#' the phenotype annotation file, kallisto QC, and FASTQC information. Defaults to 
#' 'Sample'. Not used if doAnnotation is FALSE.
#' @export
kallistoAssembleQC <- function(paired=TRUE, doAnnotation=TRUE, 
                               phenoMatchLabel="Sample") {
  
  ## get list of kallisto output directories
  kallisto_dir <-  getDir("KALLISTO")
  kalRes <- list.dirs(path=kallisto_dir, full.names=FALSE, recursive=FALSE)
  remove <- c("error", "log")
  kalRes <- kalRes[!kalRes %in% remove]  ## remove logging directories
  
  ## get kallisto QC
  kallQC <- read.csv(paste0(kallisto_dir, "/kallisto_QC.csv"), row.names=1,
                     stringsAsFactors=FALSE, header=TRUE)
  
  ## parse fastqc files for one set of paired reads
  parseFASTQC <- function(name, files=fastqc) {
    fileys <- files[grep(name, files)]
    len <- length(fileys)
    if(paired & len != 2) {
      stop(paste(name, "is paired data but does not have two fastqc files"))
    }
    if(!paired & len != 1) {
      stop(paste(name, "is single read data but does not have one fastqc file"))
    }
    
    ## extract sequence quality metric
    outie <- NULL
    for(file in fileys) {
      outie <- c(outie, read.table(file, sep="\t", stringsAsFactors=FALSE)[3,1])
    }
    return(outie)
  } ## end parseFASTQC
  
  ## get FASTQC QC
  fastqc <- list.files(getDir("FASTQC"), "summary.txt", full.names=TRUE, 
                       recursive=TRUE)
  fastqcQC <- t(apply(matrix(kalRes, ncol=1), 1, parseFASTQC))
  fastqcQC <- cbind(fastqcQC, kalRes)
  dimnames(fastqcQC) <- list(kalRes, c("AvgSeqQuality_R1", "AvgSeqQuality_R2", 
                                       "Sample"))
  
  ## combine kallistoQC and fastqcQC
  allQC <- merge(kallQC, fastqcQC, by="Sample")
  
  ## merge phenotype data with QC if desired
  if(doAnnotation) {
    
      phenoFile <- list.files(getDir("RAWANNOTATIONS"), pattern=".csv", 
                              full.names=TRUE)
      if(length(phenoFile) != 1) {
        stop("There must be a single .csv phenotype file in RAW_ANNOTATIONS directory. ")
      }
      pheno <- read.csv(phenoFile, stringsAsFactors=FALSE)
      
      if(!phenoMatchLabel %in% colnames(pheno)) {
        stop(paste(phenoMatchLabel, " column not found in ", basename(phenoFile)))
      }
      
      allQC <- merge(pheno, allQC, by.x=phenoMatchLabel, by.y="Sample")
      
  } ## end doAnnotation
  
  write.csv(allQC, paste0(getDir("OUTPUT"), "/kallisto_Pheno.csv"), row.names=FALSE)
} ## end kallistoAssembleQC

