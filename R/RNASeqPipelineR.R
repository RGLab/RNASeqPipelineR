#' Simplify pipelining RNASeq data
#' 
#' RNASeqPipeline allows you to standardize the processing of your RNASeq data.
#' 
#' RNASeqPipeline simplifies the process of pipelining the processing of RNASeq
#' data by helping you organize the raw, intermediate and processed data in a consistent
#' manner. 
#' 
#' @docType package
#' @name RNASeqPipeline
#' @author Raphael Gottardo and Greg Finak
#' @import data.table
#' @import RSQLite
#' @import GEOquery
#' @import SRAdb
NULL


dir.exists<-function (x)
{
  res <- file.exists(x) & file.info(x)$isdir
  setNames(res, x)
}

#' Create a new RNASeqPipeline project
#' 
#' Create the skeleton for a new RNASeqPipeline project.
#' 
#' createProject will create a new RNASeqPipeline 
#' project under directory 'project_name' in the path specified by 'path'.
#' The function creates the directory structure libraryd by RNASeqPipeline
#' within the new project directory, including locations for fastQ files, 
#' fastQC output, RSEM quantification, and optionally GEO and SRA files if the data
#' must be downloaded or linked to a GEO accession. Standard locations will also be provided 
#' to annotate the data utilizing the Immport schema. If the project exists, it will load the configuration 
#' if it can find it, otherwise it will proceed to re-configure the project.
#' 
#' @param project_name \code{character} name (ie, subdirectory) for project
#' @param path \code{character} Root directory for project
#' @param verbose \code{logical} Should verbose output be given?
#' @param load_from_immport \code{logical} Creates a 'Tab' directory for Immport tables if the data are loaded from Immport.
#' @param name \code{character} The name of the project directory
#' @return NULL
#' @export
#' @examples
#' # construct a projects skeleton in a new folder titled "myproject".
#' createProject("myproject",path=".")
#' createProject("myproject",path=".",verbose=TRUE)
#' \dontrun{
#' createProject("myproject",path=".",load_from_immport=TRUE) #won't work unless immport is set up
#' }
createProject <- function(project_name,path=".",verbose=FALSE, load_from_immport=FALSE){
  project_dir <- file.path(path,project_name)
  success<-FALSE
  oldwarn<-options("warn")$warn
  options("warn"=-1)
  normpath<-try(normalizePath(project_dir),silent=TRUE)
  options("warn"=oldwarn)
  if(dir.exists(normpath)){
    #load the configuration and return
    message("Project exists, loading configuration.")
    success<-readConfig(normalizePath(project_dir))      
  }
  if(success)
    invisible()
  cmnd_prefix <- "mkdir -p "
  dirs <- c(SRA="SRA",FASTQ="FASTQ",RSEM="RSEM",FASTQC="FASTQC",GEO="GEO",CONFIG="CONFIG",OUTPUT="OUTPUT",RAWANNOTATIONS="RAW_ANNOTATIONS")
  if(load_from_immport)
    dirs<-c(dirs,Tab="Tab")
  subdirs <- file.path(project_dir,dirs)
  names(subdirs)<-names(dirs)
  build_skeleton_commands<-paste0(cmnd_prefix,c(project_dir,subdirs))
  results<-sapply(build_skeleton_commands,system)
  if(verbose){
    for(i in seq_along(results)){
      cat(names(results)[i], ": ",ifelse(results[i]==0,"success!","failed!"),"\n")
    }
  }else{
    if(any(results>0)){
      sapply(names(results[results>0]),function(x){
        cat("Failed at: ",x,"\n")
      })
    }else{
      cat("Project successfully created.\n")
    }
  }
  configure_project(subdirs)
  saveConfig()
} 

#' Load an RNASeqPipelineR project
#' 
#' Load and RNASeqPipelineR project
#' 
#' Load an RNASeqPipelineR project from disk.
#' @param project_dir \code{character} the path to the project
#' @param name \code{character} the name of the project
#' @export
loadProject <- function(project_dir=NULL,name=NULL){
  if(dir.exists(normalizePath(file.path(project_dir,name)))){
    #load the configuration and return
    message("Loading configuration for ",name, " from ", project_dir);
    success<-readConfig(normalizePath(file.path(project_dir,name)))      
  }
  if(success){
    invisible()
  }else{
    stop("Can't load project.")
  }
}

#global config in package namespace
rnaseqpipeliner_configuration <- list()

#' Configure the project and store info in namespace
#' 
#' Configures the project based on the directories created in the project path
#' 
#' Stores the project configuration in the rnaseqpipeliner_configuration 
#' variable in the packge namespace.
#' Generally, this function shouldn't be called by the user.
#' 
#' @param subdirectories a \code{character} vector of subdirectories passed from \code{createProject}, or a path to an existing project.
#' @param fromDisk \code{logical} specifying whether to load the configuration of an existing proejct, in which case the path to the proejct is given by \code{subdirectories} argument
configure_project <- function(subdirectories=NULL,fromDisk=FALSE){
  ns <- getNamespace("RNASeqPipelineR")
  obj<-getConfig()
  if(!fromDisk){
    obj$subdirs <- normalizePath(subdirectories)
    names(obj$subdirs)<-names(subdirectories)
    obj$immport_files<- c(file_info="file_info.txt", expsample="expsample.txt", biosample="biosample.txt",bio_to_exp="biosample_2_expsample.txt",arm_to_subject="arm_2_subject.txt", arm_to_cohort="arm_or_cohort.txt",subject="subject.txt")
  }else{
    message("Reading existing configuration from disk not yet implemented")
  }
  unlockBinding(sym = "rnaseqpipeliner_configuration",ns)
  assign("rnaseqpipeliner_configuration", obj, envir = ns)
  #detect aspera 
  detectAspera()
  invisible()
}

#' Get the configuration for the project from the namespace
#' 
#' Retrieve the configuration for the project from the namespace
#' @export
getConfig<-function(){
  ns <- getNamespace("RNASeqPipelineR")
  get("rnaseqpipeliner_configuration", ns)
}

#' Assign a new element to the RNASeqPipelineR configuration
#' 
#' Assign value to name in `rnaseqpipeliner_configuration` object
#' 
#' Convenience function to assign a new or reassign an element named \code{name} to the \code{rnaseqpipeliner_configuration} object with the value \code{value}
#' @param name \code{character} assign element with this name
#' @param value \code{value} an R object that assigned to the 'rnaseqpipeliner_configuration' \code{list} with name 'name'
#' @export
assignConfig<-function(name,value){
  ns<-getNamespace("RNASeqPipelineR")
  unlockBinding(sym = "rnaseqpipeliner_configuration",ns)
  obj<-getConfig()  
  obj[[name]]<-value
  assign("rnaseqpipeliner_configuration",obj,ns)
  invisible(NULL)
}

#' Load data from Immport tables
#' 
#' Loads the GEO accession numbers from Immport tables
#' 
#' Load data from Immport tables in the \code{Tab} directory.
#' Stores the tables in the configuration in the namespace.
#' 
#' @param \code{integer} specifies whether to print warnings. -1 suppresses warnings (default)
#' @param verbose \code{integer} how verbose should we be about what's happening. 0 (default) is quiet.
#' @export
loadImmportTables <- function(warn=-1,verbose=0){
  savewarn<-options("warn")[[1]]
  options("warn"=warn)
  immport_files <- getConfig()$immport_files
  subdirs <- getConfig()$subdirs
  #test whether this is an immport project by looking for a Tab directory.
  if(!any(grepl("^Tab$",basename(subdirs)))){
    stop("Can't load Immport data as this project is not configured for Immport")
  }
  tabdir<-subdirs[grep("^Tab$",basename(subdirs))]
  #check if files exist
  anyerror<-FALSE
  for(i in names(immport_files)){
    itexists<-file.exists(file.path(tabdir,immport_files[i]))
    if(verbose>0)
      message("File ",immport_files[i]," in ./Tab directory ",c("doesn't exist","exists")[itexists+1])
    if(!itexists)
      anyerror<-TRUE
  }
  
  if(anyerror){
    stop("Populate the Immport Tab directory before continuing.")
  }
  
  #No errors, go ahead and load the data
  #Populate a named list
  immport_tables<-vector('list')
  for(i in names(immport_files)){
    immport_tables[[i]] <- fread(file.path(tabdir,immport_files[i]),verbose=FALSE)
  }
  
  #names are file_info, expsample, biosample, bio_to_exp, arm_to_subject, arm_to_cohort, subject
  #the following are created
  immport_tables$arm <- with(immport_tables,merge(arm_to_subject, arm_to_cohort, by="ARM_ACCESSION"))
  immport_tables$GSM_table <- with(immport_tables,expsample[DESCRIPTION%like%"GSM",GSM:=gsub(".*:", "", DESCRIPTION)][DESCRIPTION%like%"GSM"])
  immport_tables$pData <- with(immport_tables,merge(GSM_table[,c("EXPSAMPLE_ACCESSION","GSM","EXPERIMENT_ACCESSION"),with=FALSE], bio_to_exp[,c("BIOSAMPLE_ACCESSION","EXPSAMPLE_ACCESSION"), with=FALSE], by="EXPSAMPLE_ACCESSION"))
  immport_tables$pData <- with(immport_tables,merge(pData, biosample[,c("BIOSAMPLE_ACCESSION","STUDY_TIME_COLLECTED","STUDY_TIME_COLLECTED_UNIT","SUBJECT_ACCESSION"), with=FALSE], by="BIOSAMPLE_ACCESSION"))
  immport_tables$pData <- with(immport_tables,merge(pData, arm[,c("SUBJECT_ACCESSION","ARM_ACCESSION","NAME"),with=FALSE], by="SUBJECT_ACCESSION"))
  immport_tables$pData <- with(immport_tables,merge(pData, subject[,c("SUBJECT_ACCESSION", "AGE_REPORTED", "GENDER", "ETHNICITY"), with=FALSE], by="SUBJECT_ACCESSION"))
  immport_tables$pData<-with(immport_tables,setnames(pData, names(pData), tolower(names(pData))))
  
  #assign immport_tables to the configuration object
  assignConfig("immport_tables",immport_tables)
  options("warn"=savewarn)
}

#' Checks if the SRAdb sqlite database is available
#' 
#' If the database files is not available, it is downloaded.
#' 
#' Downloads and stores the SRAdb sqlite database if necessary and stores it in 
#' the Utils folder.
#' @export
getAndConnectSRAdb <- function(path=NULL){
  if(is.null(path)){
    stop("You must specify a path to the SRAmetadb.sqlite file, or a path where the file will be placed.")
  }
  
  #set the path to the Utils directory, doesn't actually have to be named Utils.
  path<-normalizePath(path)
  obj<-getConfig()[["subdirs"]]
  obj[["Utils"]]<-path
  
  #store the path
  assignConfig("subdirs",obj)
  
  #create if necessary
  system(paste0("mkdir -p ",file.path(path,"SRAdb")))  
  
  #download SRA db if necessary
  if(!file.exists(file.path(getConfig()[["subdirs"]]["Utils"],"SRAdb",'SRAmetadb.sqlite'))) {
    sqlfile <- getSRAdbFile(destdir = file.path(getConfig()[["subdirs"]]["Utils"],"SRAdb"))
  }else{
    message("SRAmetadb.sqlite found")
  }
  #connect and grab the data
  sra_con <- dbConnect(SQLite(),file.path(getConfig()[["subdirs"]]["Utils"],"SRAdb",'SRAmetadb.sqlite'))
  sra_tables <- dbListTables(sra_con)
  assignConfig("sra_tables",sra_tables)
  assignConfig("sra_con",sra_con)
}

#' Detect and configure the aspera client location
#' 
#' Detect the installation location of aspera client based on expected defaults
#' 
#' Sets the location of aspera client based on expected defaults. Or pass the installation path
#' @note Needs to be tested on linux and modified to locate the identity keys for aspera on linux
#' @param path \code{character} path to the aspera client. If NULL, will try to detect the location based on the OS
#' @export
detectAspera<-function(path=NULL){
  #get OS
  sysname=Sys.info()['sysname']
  #detect installation
  if(sysname=="Darwin"){
    homedir<-Sys.getenv("HOME")
    if(file.exists(file.path(homedir,"Applications/Aspera Connect.app"))){
      path <- file.path(homedir,"Applications/Aspera Connect.app","Contents","Resources")
    }
    else if(file.exists(file.path("/Applications/Aspera Connect.app/Contents/Resources"))){
      path <- file.path("/Applications/Aspera Connect.app/Contents/Resources")
    }else{
      stop("Can't detect aspera installation for Mac OS. Pass a path variable to configure.")
    }
    message("aspera detected at ", path)
    aspera_path=path
    assignConfig("aspera_path",aspera_path)
  }else{
    if(is.null(path)){
      message("Run detectAspera with the path argument to configure the aspera client for your OS.")
    }
    stopifnot(file.exists(file.path(path,"ascp")))
    aspera_path=path
    assignConfig("aspera_path",aspera_path)
  }
}

#' Download FastQ files based on Immport data and SRA db, unless they already exist.
#'
#' Based on the SRA db, download the FastQ files, unless they already exist.
#' 
#' Conditionally download the FastQ files in the SRA db based on Immport tables.
#' If the files are already present, they will not be downloaded.
#' 
#' @export
downloadSRA <- function(GSE_accession=NULL){
  #Two branches.. if immport_tables exists and if they don't
  if(names(getConfig())%in%"immport_tables"){
    n = nrow(getConfig()[["immport_tables"]][["GSM_table"]])
  }else{
    sra_con<-getConfig()[["sra_con"]]
    if(is.null(GSE_accession)){
      stop("Must supply a GEO accession to download SRA files.")
    }
    stop("Support for SRA from non-immport data is pending.")
  }
  cond_eval <- length(list.files(getConfig()[["subdirs"]][["SRA"]], pattern=".sra")) < n
  if(cond_eval){
    if(names(getConfig())%in%"immport_tables"){
      GSM_table<-getConfig()[["immport_tables"]][["GSM_table"]]
    }else{
      sra_con<-getConfig()[["sra_con"]]
      stop("Support for SRA from non-immport data is pending.")
    }
    successes<-0
    for(file in GSM_table[,GSM]) {
      gd <- getGEO(file, destdir=getConfig()[["subdirs"]][["GEO"]])
      SRX_number <- gsub(".*=SRX", "SRX", gd@header$relation[1])
      
      # Convert to aspera address
      sra_con<-getConfig()[["sra_con"]]
      run_accession <- listSRAfile(SRX_number, sra_con, fileType = "sra" )$run
      aspera_url <- paste0("anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra", "/", substr(run_accession,1,3), "/", substr(run_accession,1,6), "/", run_accession, "/", run_accession, ".sra")
      
      out<-system(paste0('ascp -i ',gsub(" ","\\\\ ",getConfig()[["aspera_path"]]),'/asperaweb_id_dsa.openssh -k 1 -T -l200m ', aspera_url, " ",getConfig()[["subdirs"]][["SRA"]]))
      if(out==0)
        successes<-successes+1
    }
    message("Downloaded ", successes, " files")
  }else{
    message("SRA files already downloaded")
  }
}

#' Utility function to truncate the test data to `n` SRA files for testing
#' 
#' This will truncate the GSM_table to `n` files for testing purposes.
#' 
#' Since there is over 100Gb of data, we truncate to n (default 2) files for testing.
#' @param n \code{integer} the number of files to limit
#' @export
devel_truncateData <- function(n=2){
  obj <- getConfig()
  #Truncate to download only two files
  immport_tables <- obj[["immport_tables"]]
  immport_tables[["GSM_table"]] <- immport_tables[["GSM_table"]][1:n,]
  assignConfig("immport_tables",immport_tables)
  message("Truncating data to first ", n, " SRA files.")
  invisible(NULL)
}

#' Annotate the SRA file pData
#' 
#' Set the pData to reflect the downloaded SRA files
#' 
#' @export
annotateSRA <- function(){
  pData<-getConfig()[["immport_tables"]][["pData"]]
  GSM_table<-getConfig()[["immport_tables"]][["GSM_table"]]
  pData<-pData[gsm%in%GSM_table[,GSM]]
  SRX_number<-pData[gsm%in%GSM_table[,GSM]][,gsm]
  #SRX_number <- pData[,gsm]
  for(file in pData[,gsm]) {
    gd <- getGEO(file, destdir=getConfig()[["subdirs"]][["GEO"]])
    number <- gsub(".*=SRX", "SRX", gd@header$relation[1]) 
    run_accession <- listSRAfile(number, getConfig()[["sra_con"]], fileType = "sra" )$run
    SRX_number[SRX_number==file] <- run_accession
  }
  pData$srr <- SRX_number
  
  pData[,c("expsample_accession","experiment_accession","arm_accession"):=NULL]
  setnames(pData, c("name", "ethnicity"), c("arm_name", "race"))
  setcolorder(pData, neworder = c("subject_accession", "age_reported", "gender", "race", "study_time_collected", "study_time_collected_unit", "arm_name", "gsm", "srr","biosample_accession"))
  write.csv(pData, file=file.path(getConfig()[["subdirs"]][["RSEM"]],"rsem_pdata.csv"), row.names=FALSE) 
  
  #save the changes
  immport_tables<-getConfig()[["immport_tables"]]
  immport_tables[["pData"]] <- pData
  immport_tables[["GSM_table"]] <- GSM_table
  assignConfig("immport_tables",immport_tables)
  message("Wrote pData to RSEM directory")
}

#' Overview of downloaded files
#' 
#' Prints an overview of the SRA files
#' 
#' Generates a table of the downloaded files for an overview to be included
#' in reports.
#' @export
#' @importFrom knitr kable
SRAOverview <- function(){
  pData <- getConfig()[["immport_tables"]][["pData"]]
  kable(head(pData), format = "markdown")
}

#' Convert SRA files to FastQ
#' 
#' Uses CLI SRA Toolkit to convert SRA files to FastQ if they don't already exist.
#' 
#' Converts SRA files to FastQ using the SRA Toolkit. The function checks whether the FastQ files exist, and if not, performs the conversion.
#' @export
#' @param ncores \code{integer} how many threads to use
convertSRAtoFastQ <- function(ncores=8){
  convertFastQ <- length(list.files(getConfig()[["subdirs"]][["FASTQ"]]))<length(list.files(getConfig()[["subdirs"]][["SRA"]]))
  convert_command <-  paste0("parallel -j ",ncores," fastq-dump {} -O ",getConfig()[["subdirs"]][["FASTQ"]]," ::: ", file.path(getConfig()[["subdirs"]][["SRA"]],"*.sra"))
  if(convertFastQ){
    out<-system(convert_command)
    if(out==0)
      message("Finished converting ",length(list.files(getConfig()[["subdirs"]][["SRA"]])), " files.")
    else
      warning("Something went wrong when converting FastQ files.")
  }else{
    message("FASTQ Files already present.")
  }
  
}

#' Concatenate FastQ files
#'
#' Combine several fastq.gz files into one fastq.gz file for each library
#' Decompress the fastq.gz 
#' Copy the fastQ files into FASTQ folder
#' @export
#' @param infile \code{character} specify the path to the individual fastq.gz directory.
#' @param outfile \code{character} specify the path to the FASTQ directory
#' @param pattern \code{character} the string want to remove in the fastq.gz files name
concatenateFastq = function(infile, outfile, pattern)
{
  samples <- system(paste0("ls ", infile), intern = T)
  sample.name <- unique(sub(pattern, "", samples))
  fastq.name <- paste0(sample.name, ".fastq.gz")
  system(paste0("cat ", infile,  "/", sample.name, "_*.fastq.gz > ", infile, "/", fastq.name))
  system(paste0("gunzip ", infile,  "/", fastq.name))
  system(paste0("mv ", infile, "/*.fastq ", outfile))
}

#' Perform FastQC quality control
#' 
#' Runs fastqc quality control on the FASTQ files
#' 
#' Produces FASTQC reports for each FASTQ file using fastqc
#' @export
#' @param ncores \code{integer} how many threads to use
runFastQC <- function(ncores=8){
  fastQCL <- length(list.files(getConfig()[["subdirs"]][["FASTQC"]]))<length(list.files(getConfig()[["subdirs"]][["FASTQ"]]))
  run_command <-  paste0('parallel -j ', ncores,' fastqc {} -o "',getConfig()[['subdirs']][['FASTQC']],'" -q ::: "',file.path(getConfig()[['subdirs']][['FASTQ']],'*.fastq"'))                             
  if(fastQCL|length(list.files(getConfig()[["subdirs"]][["FASTQC"]]))==0){
    out<-system(run_command)
    if(out==0){
      message("Finished fastqc process for ",length(list.files(getConfig()[["subdirs"]][["FASTQ"]])), " files.")
      message("Expanding archives")
      f<-list.files(pattern="zip$",path=getConfig()[["subdirs"]][["FASTQC"]])
      sapply(f,function(x){
        out<-system(paste0(paste0("cd ",getConfig()[["subdirs"]][["FASTQC"]]),"&& unzip -o ",x),intern=FALSE,ignore.stdout=TRUE,ignore.stderr=FALSE,wait=TRUE)
      })
      if(any(out!=0)){
        warning("encountered an error unzipping fastqc reports")
      }
    }
    else
      warning("Something went wrong when running fastqc")
  }else{
    message("FASTQC already run")  
  }
}

.collectFastqcFiles <- function(pattern){
  # Put the FASTQC results together
  # List sudirs
  fastqc_subdirs <- grep("_fastqc", list.dirs(getConfig()[["subdirs"]][["FASTQC"]], recursive = FALSE), value=TRUE)
  # Find reports
  sapply(fastqc_subdirs, list.files, pattern=pattern, full.names = TRUE)

}

#' Summarize the fastQC results
#' 
#' Summarize the fastqc results
#' 
#' Function to be run after `runFastQC` that will summarize the fastQC results
#' and generate a plot.
#'
#' @return prints a plot and invisibly returns the data.table used to generate the plot
#' @export
#' @import ggplot2
summarizeFastQC <- function(){
 
  fastqc_summary_files <- .collectFastqcFiles("summary.txt")
  # Read all files and create a list of data.tables
  fastqc_list <- lapply(fastqc_summary_files, function(x, ...){
    y <- fread(x, header=FALSE);
    setnames(y, colnames(y), c("result", "test", "file"))
    return(y);})
  
  fastqc_data <- rbindlist(fastqc_list)
  fastqc_data[,qresult:=ifelse(result=="FAIL",0, ifelse(result=="PASS", 1, 1/2))]
  fastqc_data[,sum_qresult:=sum(qresult), by="file"]
  fastqc_data <- fastqc_data[order(-sum_qresult)]
  
  # Let's reorder by overall fastqc score
  fastqc_data$file <- reorder(fastqc_data$file, -fastqc_data$sum_qresult)
  print(ggplot(fastqc_data, aes(x=test, y=file, fill=result))+geom_raster()+theme(axis.text.x=element_text(angle=90,hjust=1)))
  invisible(fastqc_data)
  #tidy up
  #sapply(fastqc_subdirs,function(x)system(paste0("rm -rf ",x)))
}

##' Summarize duplication statistics from fastqc
##'
##' Prints plots of the duplication distribution and % duplication across libraries
##' @return data.table of source data, invisbly
##' @export
summarizeDuplication <- function(){
    fastqc_complete <- .collectFastqcFiles("fastqc_data.txt")
    fastqc_list <- lapply(fastqc_complete, function(x, ...){
        header <- read.table(x, skip=3, nrows=7, sep='\t', colClasses='character')
        setDT(header)
        setnames(header, c('record', 'value'))
        setkey(header, record)
        totaldup <- fread(x, skip='#Total Deduplicated Percentage', nrows=1, sep='\t', header=F, colClasses=c('character', 'numeric'))
        y <- fread(x, skip='#Duplication Level', nrows=15, sep='\t', autostart=4, header=TRUE, colClasses=c('character', 'numeric', 'numeric'))        
        setnames(y, colnames(y), c("duplication level", "pct dedup", "pct total"))
        y[,file:=header['Filename', value]]
        y[,totaldup:=100-totaldup[1,2,with=FALSE]]
        return(y)})
    fastqc_data <- rbindlist(fastqc_list, fill=TRUE)
    fastqc_data[,`duplication level`:=factor(`duplication level`, levels=c(1:9, c('>10', '>50', '>100', '>500', '>1k', '>5k')))]
    M <- data.table:::melt.data.table(fastqc_data, measure.vars=c('pct dedup', 'pct total'))
    print(ggplot(M, aes(x=`duplication level`, y=value))+geom_boxplot() + facet_wrap(~variable))
    U <- unique(fastqc_data[,list(file, totaldup)])
    print(ggplot(U, aes(x=totaldup)) + geom_density() + geom_text(aes(x=totaldup, y=0, label=file), size=2, angle=90, hjust=0))
    invisible(fastqc_data)
}

#' Build a reference genome
#' 
#' Builds a reference genome at `Utils/Reference_Genome`
#' 
#' You must specify the Utils path if it is not already defined, and have your genome in a folder titled
#' `Reference_Genome`. This function will construct the reference genome using RSEM tools.
#' The command line is the default shown in the documentation.
#' `rsem-prepare-reference --gtf gtf_file --transcript-to-gene-map knownIsoforms.txt --bowtie2 fasta_file name`
#' If the gtf_file is not give, then the transcript-to-gene-map option is not used either. A fasta_file and a name must be provided.
#' @param path \code{character} specifying an \emph{absolute path} path to the Utils directory.
#' @param gtf_file \code{character} the name of the gtf file. Empty by default. If specified the function will look for a file named `knownIsoforms.txt`
#' @param fasta_file \code{character} the name of the fasta file, must be specified
#' @param name \code{character} the name of the genome output.
#' @export
buildReference <- function(path=NULL,gtf_file="",fasta_file=NULL,name=NULL){
  if(is.null(fasta_file)|is.null(name))
    stop("You must provide a fasta_file  and a genome name")
  if(is.null(path)&inherits(try(getConfig()[["subdirs"]][["Utils"]],silent=TRUE),"try-error")){
    stop("Please specify where to build the reference genome")
  }else if(is.null(path)){
    path = file.path(getConfig()[["subdirs"]][["Utils"]],"Reference_Genome")
  } else{
    subdirs<-getConfig()[["subdirs"]]
    subdirs[["Utils"]]<-path
    assignConfig("subdirs",subdirs) #save the user-provided path
    path = file.path(path,"Reference_Genome")
  }
  if(substr(path, 1, 1) != '/') stop("'path' must be an absolute path")
  if(gtf_file==""){
    gtfopt<-""
    isoforms<-""
    isoformsOpt<-""
  }
  else{
    gtfopt <- "--gtf"
    isoforms<-"knownIsoforms.txt"
    isoformsOpt<-"--transcript-to-gene-map"
  }
  #check if the reference genome has already been built
  if(length(list.files(pattern=paste0(name,".chrlist"),path=file.path(getConfig()[["subdirs"]]["Utils"],"Reference_Genome")))==0){
    command = paste0("cd ",path," && rsem-prepare-reference ",gtfopt," ",gtf_file," ",isoformsOpt, " ",isoforms," --bowtie2 ",fasta_file," ",name," ")
    system(command)
  }else{
    message("Reference Genome Found")
  }
  #set the reference genome name
  assignConfig("reference_genome_name",name)
}

#' Use the RSEM tool to align the reads
#' 
#' Use the RSEM tool to align reads
#' 
#' Uses RSEM to align reads in FASTQ files against the reference genome. Optionally you can specify paired end reads. The code assumes 
#' paired reads have fastq files that differ by one character (i.e. sampleA_read1.fastq, sampleA_read2.fastq) and will perform
#' matching of paired fastq files based on that assumption using string edit distance. Read 1 is assumed to be upstream
#' and read 2 is assumed to be downstream. 
#' The number of parallel_threads*bowtie_threads should not be more than the number of cores available on your system.
#' @param parallel_threads \code{integer} specify how many parallel processes to spawn
#' @param paired \code{logical} specify whether you have paried reads or not.
#' @param bowtie_threads \code{integer} specify how many threads bowtie should use.
#' @export
RSEMCalculateExpression <- function(parallel_threads=2,bowtie_threads=4,paired=FALSE){
  suppressPackageStartupMessages(library(parallel))
  if(parallel_threads*bowtie_threads>detectCores()){
    stop("The number of parallel_threads*bowite_threads is more than the number of cores detected by detectCores()")
  }
  rsem_dir <- getConfig()[["subdirs"]][["RSEM"]]
  fastq_dir <- getConfig()[["subdirs"]][["FASTQ"]]
  lr <- length(list.files(path=rsem_dir,pattern=".genes.results"))
  lq <- length(list.files(path=fastq_dir,pattern=".fastq"))
  reference_genome_path <- file.path(getConfig()[["subdirs"]][["Utils"]],"Reference_Genome")
  reference_genome_name <- file.path(getConfig()[["reference_genome_name"]])
  if(lr!=lq){
    #do a set difference
    keep<-setdiff(gsub("\\.fastq","",list.files(path=fastq_dir,pattern=".fastq")),gsub("\\.genes\\.results","",list.files(path=rsem_dir,pattern=".genes.results")))    
    if(!paired){
      keep<-paste0(keep,".fastq")
      myfiles<-file.path(fastq_dir,keep)
      command <- paste0("cd ",rsem_dir," && parallel -j ",parallel_threads," rsem-calculate-expression --bowtie2 -p ",bowtie_threads," {} ",file.path(reference_genome_path,reference_genome_name)," {/.} ::: ", paste(myfiles, collapse=' '))
    }else{
      keep<-paste0(keep,".fastq")
      fastq_files<-file.path(fastq_dir,keep)
      #fastq_files<-list.files(fastq_dir,"*.fastq",full=TRUE)
      library(clue)
      D<-adist(fastq_files)
      diag(D)<-1e20
      pairs<-solve_LSAP(D)
      pairs<-matrix(fastq_files[pairs],ncol=2,byrow=TRUE)
      pairs <- t(apply(pairs,1,sort)) #Sort rows lexicographically, assumption is that they differ by numeric index 1, 2, for paired reads
      pairs<-cbind(pairs,basename(pairs[,1]))
      writeLines(t(pairs),con=file(file.path(getConfig()[["subdirs"]][["FASTQ"]],"arguments.txt")))
      command<-paste0("cd ",rsem_dir," && parallel -n 3 -j ",parallel_threads," rsem-calculate-expression --bowtie2 -p ", bowtie_threads," --paired-end {1} {2} ",file.path(reference_genome_path,reference_genome_name)," {3.} :::: ",file.path(getConfig()[["subdirs"]][["FASTQ"]],"arguments.txt"))
    }
    cat(command)
    system(command)
  }else{
    message("Expression already calculated")
  }
}

#' Assemble an expression matrix of all results
#' 
#' Put all the counts from the individual libraries into a single matrix result
#' 
#' Assemble an expression matrix from the individual libraries.
#'@export 
RSEMAssembleExpressionMatrix <- function(){
  cond_eval <- length(list.files(getConfig()[["subdirs"]][["RSEM"]], pattern="rsem_"))<4
  if(cond_eval){
    message("Assembling counts matrix")
    rsem_files <- list.files(getConfig()[["subdirs"]][["RSEM"]], pattern="genes.results", full.names = TRUE)
    # Read all files and create a list of data.tables
    rsem_list <- lapply(rsem_files, function(x, ...){
      y<-fread(x, drop=c("length", "effective_length", "FPKM")); 
      y$sample_name <- gsub(".genes.results", "", basename(x));
      return(y);})
    
    # Bind all the data.tables to create a long data.table
    rsem_data_long <- rbindlist(rsem_list)
    
    # Rename a column, so that it's a bit R friendly
    setnames(rsem_data_long, "transcript_id(s)", "transcript_ids")
    
    # Reshape the long data.table to create a matrix
    # Assume that missing entries would have a TPM value of 0
    rsem_tpm_matrix <- dcast.data.table(rsem_data_long, gene_id+transcript_ids~sample_name, value.var = "TPM", fill=0)
    rsem_count_matrix <- dcast.data.table(rsem_data_long, gene_id+transcript_ids~sample_name, value.var = "expected_count", fill=0)
    
    # Keep track of the gene_id - transcript mapping 
    rsem_txs_table <- rsem_tpm_matrix[, c("gene_id","transcript_ids"), with=FALSE]
    
    # Remove the transcript column (we create an annotation table in the next chunk)
    rsem_tpm_matrix <- rsem_tpm_matrix[, transcript_ids:=NULL][order(gene_id)]
    rsem_count_matrix <- rsem_count_matrix[, transcript_ids:=NULL][order(gene_id)]
    
    # Write to disk
    write.csv(rsem_tpm_matrix, file=file.path(getConfig()[["subdirs"]][["RSEM"]],"rsem_tpm_matrix.csv"), row.names=FALSE)
    write.csv(rsem_count_matrix, file=file.path(getConfig()[["subdirs"]][["RSEM"]],"rsem_count_matrix.csv"), row.names=FALSE)
    write.csv(rsem_txs_table,file=file.path(getConfig()[["subdirs"]][["RSEM"]],"rsem_txs_table.csv"))
  }
}

#' Annotate the features using BioConductor annotation packages
#' 
#' Annotate the transcripts using bioConductor annotation packages. 
#' 
#' Annotates the transcripts using the bioconductor annotation package specified in the
#' annotation_library argument
#' @param annotation_library \code{character} specifying the annotation package to use. "TxDb.Hsapiens.UCSC.hg38.knownGene" by default.
#' @param force \code{logical} force the annotation step to re-run
#' @export
BioCAnnotate<-function(annotation_library="TxDb.Hsapiens.UCSC.hg38.knownGene",force=FALSE){
  featuredata_outfile<-"rsem_fdata.csv"
  if(!force&file.exists(file.path(getConfig()[["subdirs"]][["RSEM"]],featuredata_outfile))){
    message("Annotation already done, skipping. Use force=TRUE to rerun.")
    invisible(return(0))
  }
  do.call(library,(list(eval(annotation_library))))
  #save the annotation library
  assignConfig("annotation_library",annotation_library)
  library(annotate)
  library(org.Hs.eg.db)
  message("Annotating transcripts.")
  txdb <- get(annotation_library)
  rsem_txs_table <- fread(file.path(getConfig()[["subdirs"]][["RSEM"]],"rsem_txs_table.csv"))
  
  # create a table with one transcript per line
  tx_to_gid <- rsem_txs_table[,list(transcript_id=strsplit(as.character(transcript_ids),",")[[1]]),by="gene_id"]
  
  # map the transcripts to entrez gene ids
  tx_to_eid <- data.table(select(txdb, keys = tx_to_gid[,transcript_id], columns="GENEID", keytype="TXNAME"))
  setnames(tx_to_eid, c("TXNAME", "GENEID"), c("transcript_id", "entrez_id"))
  
  # Add gene symbol
  tx_to_eid[!is.na(entrez_id),gene_symbol:=getSYMBOL(tx_to_eid[!is.na(entrez_id),entrez_id], data='org.Hs.eg')]
  
  # merge all to map information to RSEM gene_ids
  tx_table <- merge(tx_to_gid, tx_to_eid, by="transcript_id")
  tx_table_reduced <- tx_table[,list(entrez_id=paste(unique(entrez_id),collapse=","), transcript_id=paste(unique(transcript_id),collapse=","), gene_symbol=paste(unique(gene_symbol),collapse=",")),by="gene_id"]
  
  # Write to disc and order by gene_id
  message("Writing feature data to rsem_fdata.csv")
  write.csv(tx_table_reduced[order(gene_id)], file=file.path(getConfig()[["subdirs"]][["RSEM"]],featuredata_outfile), row.names=FALSE)
}

#' Produce a report listing the tools and packages used by the pipeline.
#' 
#' Produce a report listing the tools and packages used by the pipeline.
#' 
#' Helps with reproducibility by producing a report listing the tools and packages used by the pipeline.
#' This needs to be run on the system that produced the results. Output to project directory OUTPUT/pipeline_report.md.
#' @export
pipelineReport<-function(){
  library(pander)
  library(plyr)
  #save the configuration, since we're probably done.
  saveConfig()
  command <- c("ascp --version",
               "fastq-dump --version",
               "bowtie2 --version",
               "fastqc --version",
               "rsem-calculate-expression --version")
  versions<-sapply(command,function(x)system(x,intern=TRUE)) 
  session<-capture.output(sessionInfo())
  sapply(versions,function(x)cat(paste0(x,"\n"),"\n"))
  cat(paste0(session,"\n"))
  #output to OUTPUT/pipeline_report.md
  f<-file(open = "a",description = (file.path(getConfig()[["subdirs"]][["OUTPUT"]],"pipeline_report.md")))
  sapply(versions,function(x)cat(paste0(x,"\n"),"\n",file = f))
  cat(paste0(session,"\n"),file = f)
  close(f)
}

#' Save the configuration for a project
#' 
#' Save the configuration for a project
#' 
#' Save the configuration from a project, stored in the CONFIG directory
#' 
#' @export
saveConfig <- function(){
  #TODO blindly stores configuration info. 
  #Wanto to use YAML eventually.
  config<-getConfig(); 
  saveRDS(config,file=file.path(config[["subdirs"]][["CONFIG"]],"configuration.rds"))
  invisible(TRUE)
} 

#' Read the configuration for a project
#' 
#' Read the configuration for a project
#' 
#' Read the configuration from a project, stored in the CONFIG directory
#' 
#' @param project \code{character} the path to the project folder
#' @export
readConfig <- function(project=NULL){
  #TODO this just blindly reads the config info. Some error checking needs to be done, ensuring the configuration matches what's in the project directory.
  #Ideally should call configure_project as well. 
  #Want to also store human-readable config information ultimately, perhaps using YAML.
  confdir <- file.path(project,"CONFIG")
  error<-length(list.files(confdir,pattern="configuration.rds"))!=1
  if(error){
    message("No configuration information found")
    return(invisible(FALSE))
  }
  obj<-readRDS(list.files(confdir,pattern="configuration.rds",full=TRUE))
  
  #store in the namespace
  ns <- getNamespace("RNASeqPipelineR")
  unlockBinding(sym = "rnaseqpipeliner_configuration",ns)
  assign("rnaseqpipeliner_configuration", obj, envir = ns) 
  invisible(TRUE)
}

#' Return and expression set of the counts or TPM values
#' 
#' Construct and return an expression set of counts or TPM values
#' 
#' Constructs an expression set of counts or tpm values depending on the value of which parameter.
#' @param which \code{character} \code{c("counts","tpm")} specifies what to return.
#' @export
getExpressionSet <- function(which="counts"){
  which<-match.arg(arg = which, c("counts","tpm"))
  if(which%in%"counts")
    mat <- fread(list.files(getConfig()[["subdirs"]][["RSEM"]],pattern="rsem_count_matrix.csv",full=TRUE))
  else
    mat <- fread(list.files(getConfig()[["subdirs"]][["RSEM"]],pattern="rsem_tpm_matrix.csv",full=TRUE))
  featuredata <- fread(list.files(getConfig()[["subdirs"]][["RSEM"]],pattern="rsem_fdata.csv",full=TRUE))
  pdata <- fread(list.files(getConfig()[["subdirs"]][["RSEM"]],pattern="rsem_pdata.csv",full=TRUE))
  
  #Construct an Eset and return
  pdata<-data.frame(pdata)
  rownames(pdata) <- pdata$srr  
  pdata<-AnnotatedDataFrame(data.frame(pdata))
  featuredata<-AnnotatedDataFrame(data.frame(featuredata))
  matr<-as.matrix(mat[,-1,with=FALSE])
  rownames(matr)<-mat[,gene_id]
  row.names(pdata)
  matr<-matr[,rownames(pdata),drop=FALSE]
  eset<-ExpressionSet(assayData = matr,phenoData = pdata, featureData = featuredata, annotation = getConfig()[["annotation_library"]])
  eset
}

#' Run MiTCR on each fastQ file
#' 
#' Run the MiTCR tool on each fastQ file
#' 
#' Runs MiTCR on each fastQ file. 
#' @param gene \code{character} c("TRB","TRA")
#' @param species \code{character} c("hs","mm")
#' @param ec \code{integer}, c(0,1,2)
#' @param pset \code{character} c("flex")
#' @param ncores \code{integer} number of cores for running in parallel
#' @param output_format \code{character} either "txt" or "cls". cls files can be viewed in the MiTCR viewer. Txt files can be parsed and used to annotate libraries.
#' @param paired \code{logical} specify whether data is paired (in which case the fastq files in PEAR directory are used). Defaults to FALSE.
#' @export
MiTCR <- function(gene="TRB",species=NULL,ec=2,pset="flex",ncores=1,output_format="text",paired=FALSE){
  pset<-match.arg(pset,"flex")
  output_format<-match.arg(output_format,c("txt","cls"))
  gene<-match.arg(gene,c("TRB","TRA"))
  species<-match.arg(species,c("hs","mm",NULL))
  if(!ec%in%c(0,1,2)){
    stop("Invalid value of ec")
  }
  if(length(system("which mitcr",intern=TRUE))==0){
    stop("mitcr can't be found on the path")
  }
  if(output_format=="txt")
    output_format<-".txt"
  else
    output_format<-".cls"
  command <-  paste0("mitcr -pset ",pset)
  command <- paste0(command," -gene ",gene," ")
  if(!is.null(species)){
    command <- paste0(command,"-species ",species," ")
  }
  if(!is.null(ec)){
    command <- paste0(command,"-ec ",ec," ")
  }
  #Run on each fastq file and output to TCR directory
  if(paired){
    fastqdir <- getConfig()[["subdirs"]][["PEAR"]]
    if(!dir.exists(fastqdir))
      stop("paired is TRUE, but PEAR directory not found. Did you run pear?")
  }else{
    fastqdir <- getConfig()[["subdirs"]][["FASTQ"]]
  }
  tcrdir <- file.path(dirname(fastqdir),"TCR")
  system(paste0("mkdir -p ",tcrdir))
  if(!paired){
    fastqfiles<-list.files(fastqdir,pattern="fastq$",full=TRUE)
  }else{
    fastqfiles<-list.files(fastqdir,pattern="\\.assembled\\.fastq$",full=TRUE)
  }
  if(ncores>1){
    outpath<-NULL
    for(i in fastqfiles){
      outpath<-c(outpath,file.path(tcrdir,paste0(gsub("fastq",gene,basename(i)),output_format)))
    }
    inargs<-paste(as.vector(apply(cbind(fastqfiles,outpath),1,function(x)paste(x,collapse=" "))),collapse=" ")
    system(paste0("parallel -j ",ncores," -n 2 ",command," {} ::: ",inargs))
  }else{
    for(i in fastqfiles){
      outpath<-file.path(tcrdir,paste0(gsub("fastq",gene,basename(i)),output_format))
      commandi<-paste0(command, i, " ", outpath)
      system(commandi)
    }
  } 
}

#' Run PEAR to assemble paired-end fastq files into one files.
#' 
#' Run the pear tool to assemble paired end fastq files into a single output fastq file.
#' 
#' Use this as a preliminary preprocessing step to running MiTCR if you have paired-end data.
#' @param ncores \code{integer} number of cores to use
#' @export
pear<-function(ncores=4){
  fastq_dir <- getConfig()[["subdirs"]][["FASTQ"]]
  pear_directory<-try(getConfig()[["subdirs"]][["PEAR"]],silent=TRUE)
  if(inherits(pear_directory,"try-error")){
    pear_directory<-file.path(dirname(getConfig()[["subdirs"]][["FASTQ"]]),"PEAR")
    dir.create(pear_directory)
    subdirs<-getConfig()[["subdirs"]]
    subdirs[["PEAR"]]<-pear_directory
    assignConfig("subdirs",subdirs)
    saveConfig()
  }
  pear_directory<-getConfig()[["subdirs"]][["PEAR"]]
  #We just assume that paired fastq files differ by one character, that there are an even number of them, and that the operating system
  #will return them to us in lexicographical order. 
  files<-(list.files(path=fastq_dir,pattern="*.fastq",full=FALSE))
  f<-file.path(pear_directory,"pear_arguments.txt")
  connection<-file(f,open="w")
  writeLines(files,con = connection)
  close(connection)
  command<-paste0("cd ", fastq_dir, " && parallel -j ",ncores, " -n2 pear -f {1} -r {2} -o ",file.path(pear_directory,"{1}")," :::: < ",file.path(pear_directory,"pear_arguments.txt"))
  system(command)
}
