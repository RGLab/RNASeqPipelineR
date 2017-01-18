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
#' @import stringr
NULL

## silence complaints about variables not found in calls that use non-standard evaluation
if(getRversion() >= "2.15.1") globalVariables(c(
                  'transcript_ids', # RSEMAssembleExpressionMatrix, BioCAnnotate
                  'transcript_id',  #BioCAnnotate
                  'entrez_id',
                  'gene_symbol',
                  'gene_id',
                  'record', #summarizeDuplication
                  'value',
                  'duplication level',
                  'totaldup',
                  'qresult', #summarizeFastQC
                  'result',
                  'sum_qresult',
                  'test'))
                  


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
  dirs <- c(SRA="SRA",FASTQ="FASTQ",RSEM="RSEM",FASTQC="FASTQC",GEO="GEO",CONFIG="CONFIG",OUTPUT="OUTPUT",RAWANNOTATIONS="RAW_ANNOTATIONS",RNASEQC="RNASEQC",TOPHAT="TOPHAT", BAM='BAM')
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
#' Retrieve the configuration for the project from the namespace. If
#' configuration not found exit with error.
#' 
#' @return Returns configuration structure or stops with error if not found
#' @export
getConfig<-function(){
  ns <- getNamespace("RNASeqPipelineR")
  config <- get("rnaseqpipeliner_configuration", ns)
  if(length(config) == 0) {
       stop("Configuration not found: Have you run loadProject or buildReference (for Utils directory)?")
  }
  return(config)
}


#' Retrieve full path to 'dir' in project configuration structure
#'
#' Retrieve full path to 'dir' in project configuration structure. If 'dir'
#' not found exit with error.
#'
#' @param dir \code{character} directory that path is required for
#' @return \code{character} full path to directory
getConfigDir <- function(dir) {

    path <- getConfig()[["subdirs"]]

    if(is.null(path)) stop("subdirs list not found. Configuration corrupted!")
    if( is.na(match(dir, names(path))) ) {
        stop(paste("directory", dir, "not found in project configuration"))
    }

    return(getConfig()[["subdirs"]][[dir]])
} ## end getConfigDir


#' Retrieve full path to 'dir/file' in project configuration structure
#'
#' Retrieve full path to 'dir/file' in project configuration structure. If
#' file not found exit with error.
#'
#' @param dir \code{character} directory containing file
#' @param file \code{character} file sought after
#' @return \code{character} full path to directory
getConfigFile <- function(dir, file) {

    ## construct path
    path <- getConfigDir(dir)
    filey <- paste(path, file, sep="/")

    ## error check
    if(!file.exists(filey)) {
        stop(paste(filey, "not found"))
    }
  
    ## return file path
    return(filey)
} ## end getConfigFile


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
#' @param warn \code{integer} specifies whether to print warnings. -1 suppresses warnings (default)
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
#' @param path \code{character} path to SRAmetadb.sqlite file exists or will be placed
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
#' @param GSE_accession (not implemented)
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
  kable(utils::head(pData), format = "markdown")
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

stripExtension <- function(finame, pattern='([.][A-Za-z]+$)'){
    sub(pattern, '', finame)
}

#' Perform FastQC quality control
#' 
#' Runs fastqc quality control on the FASTQ files
#' 
#' Produces FASTQC reports for each FASTQ file using fastqc
#' FASTQ input files must have extension \code{fastqc} (capitalization sensitive).
#' @export
#' @param ncores \code{integer} how many threads to use
runFastQC <- function(ncores=8){
  fastQCout <- list.files(getConfig()[["subdirs"]][["FASTQC"]])
  stripFQC <- unique(stripExtension(fastQCout, '(_fastqc.*$)')) #expanded fastqc generates at least 3 files per input
  FQfile <- list.files(getConfig()[["subdirs"]][["FASTQ"]], pattern="*.fastq", ignore.case=TRUE, full.names=TRUE)
  FQfile <- data.frame(file=FQfile, done=stripExtension(basename(FQfile)) %in% stripFQC, stringsAsFactors=FALSE)
  notrun <- FQfile[!FQfile$done,]
  run_command <-  paste0('parallel -j ', ncores,
                         ' fastqc {} -o "',getConfig()[['subdirs']][['FASTQC']],
                         '" -q ::: "', paste(notrun$file, collapse='" "'), '"')                             
  if(nrow(notrun) > 0){
      message('Running fastqc process on ', nrow(notrun), ' files.')
    out<-system(run_command)
    if(out==0){
      message("Finished fastqc process for ", nrow(notrun), " files.")
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
  M <- data.table::melt.data.table(fastqc_data, measure.vars=c('pct dedup', 'pct total'))
  print(ggplot(M, aes(x=`duplication level`, y=value))+geom_boxplot() + facet_wrap(~variable))
  U <- unique(fastqc_data[,list(file, totaldup)])
  print(ggplot(U, aes(x=totaldup)) + geom_density() + geom_text(aes(x=totaldup, y=0, label=file), size=2, angle=90, hjust=0))
  invisible(fastqc_data)
}


#' Build a reference genome 2
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
#' @param star_threads \code{integer} number of threads to use when building index
#' @param name \code{character} the name of the genome output.
#' @param gff3 \code{boolean} when gtf_file is gff3 format set it to TRUE to add "--sjdbGTFtagExonParentTranscript gene".
#' @param additional_param \code{character} default "".  Additional parameters to pass.
#' @export

buildGenomeIndexSTAR2 = function (path = NULL, gtf_file = "", fasta_file = NULL, star_threads = 1, name=NULL, gff3 = FALSE, additional_param = ""){
  if (is.null(fasta_file)){
    stop("You must provide a fasta_file")
  }

  if (is.null(path) & inherits(try(getConfig()[["subdirs"]][["Utils"]], silent = TRUE), "try-error")) {
    stop("Please specify where to build the reference genome")
  }
  else if (is.null(path)) {
    path <- file.path(getConfig()[["subdirs"]][["Utils"]], "Reference_Genome/STARIndex")
  }
  else {
    subdirs <- getConfig()[["subdirs"]]
    subdirs[["Utils"]] <- path
    assignConfig("subdirs", subdirs)
    path <- file.path(getConfig()[["subdirs"]][["Utils"]], "Reference_Genome/STARIndex")
  }
  
  if (substr(path, 1, 1) != "/") {
    stop("'path' must be an absolute path")
  }
  dir.create(path, showWarnings=FALSE)
  if (gtf_file == "") {
    gtfopt <- ""
  }
  else {
    gtfopt <- "--sjdbGTFfile"
    gtf.file <- file.path(paste0(getConfig()[["subdirs"]][["Utils"]], "/Reference_Genome/", gtf_file))
  }
  gffopt=""
  if(gff3==TRUE){
    gffopt<- "--sjdbGTFtagExonParentTranscript gene"
  }
  fasta.file <- file.path(paste0(getConfig()[["subdirs"]][["Utils"]], "/Reference_Genome/", fasta_file))
  if (length(fasta.file) > 1) {
    fasta.file <- paste(fasta.file, collapse = ",")
  }
  if (length(list.files(pattern="SAindex", path=path)) == 0) {
    command = paste0("STAR --runThreadN ", star_threads, " --runMode genomeGenerate --genomeDir ", path, " --genomeFastaFiles ", fasta.file, " ", gtfopt, " ", gtf.file, " ", gffopt , " ", additional_param  )
    system(command)
  }
  else {
    message("Reference Genome Found")
  }
  assignConfig("reference_genome_name", name)
}

#' Use the RSEM tool to annotate and quantify the reads
#'
#' Use the RSEM tool to annotate and quantify the reads
#'
#' Uses RSEM to quantify reads in FASTQ files against the reference genome. Optionally you can specify paired end reads. The code assumes
#' paired reads have fastq files that differ by one character (i.e. sampleA_read1.fastq, sampleA_read2.fastq) and will perform
#' matching of paired fastq files based on that assumption using string edit distance. Read 1 is assumed to be upstream
#' and read 2 is assumed to be downstream.
#' The number of parallel_threads*bowtie_threads should not be more than the number of cores available on your system.
#' @param parallel_threads \code{integer} specify how many parallel processes to spawn
#' @param bowtie_threads \code{integer} specify how many threads bowtie should use.
#' @param paired \code{logical} specify whether you have paried reads or not.
#' @param frag_mean \code{numeric}
#' @param frag_sd \code{numeric} For single ended reads, specifying these might make calculations of effect length more effective. Optional.
#' @param nchunks \code{integer} number of chunks to split the files for a slurm job. Ignored if slurm = FALSE
#' @param days_requested \code{integer} number of days requested for the job (when submitting a slurm job). Ignored if slurm = FALSE
#' @param slurm \code{logical} if \code{TRUE} job is submitted as a slurm batch job, otherwise it's run on the local machine. Slurm jobs will honour the nchunks and days_requested arguments.
#' @param slurm_partition \code{character} the slurm partition to submit to. Ignored if slurm=FALSE
#' @param ram_per_node \code{numeric} The number of Mb per node. Ignored if slurm=FALSE. Default of \code{parallel_threads*bowtie_threads*1000}
#' @param fromBAM \code{logical} if \code{TRUE} then RSEM will attempt to use previously aligned BAM files, in the \code{BAM} directory, rather than fastq files. The file names expected to end with \code{.transcript.bam}. See RSEM documentation for the format these files must obey.
#' @param fromSTAR \code{logical} if \code{TRUE} then RSEM will use the
#' STAR notation for the BAM files
#' @note The amount of memory requested should be set to bowtie_threads*parallel_threads*1G as this is the default requested by samtools for sorting. If insufficient memory is requested, the bam files will not be created successfully.
#' @export
RSEMCalculateExpression <- function(parallel_threads=6,bowtie_threads=1,paired=FALSE, frag_mean=NULL, frag_sd=NULL,
                                    nchunks=10,days_requested=5,slurm=FALSE, slurm_partition=NULL, #slurm_partition="gottardo_r",
                                    ram_per_node=bowtie_threads*parallel_threads*1200, fromBAM=FALSE, fromSTAR=FALSE){
  ncores<-parallel_threads*bowtie_threads
  if(ncores>parallel::detectCores()&!slurm){
    stop("The number of parallel_threads*bowite_threads is more than the number of cores detected by detectCores() on the local machine for non-slurm jobs")
  }
  #Chunking for slurm
  .chunkDataFrame<-function(df=pairs,nchunks=nchunks){
    groupsize<-nrow(df)%/%nchunks
    split(as.data.frame(df),gl(nchunks,groupsize,length=nrow(df)))
  }
  rsem_dir <- getConfig()[["subdirs"]][["RSEM"]]
  ## parse RSEM output to get names of fastq files which have already been annotated.
  done <- str_replace(list.files(path=rsem_dir,pattern="\\.genes\\.results$"), 'Aligned\\.toTranscriptome\\.out\\.genes\\.results$', '')

  #browser()
  if(!fromBAM){
    fastq_dir <- getConfig()[["subdirs"]][["FASTQ"]]
    todo <- unique(str_replace(list.files(path=fastq_dir,pattern="\\.fastq$"), '(_[12])?\\.fastq$', ''))
    extension <- if(paired) c('_1.fastq', '_2.fastq') else '.fastq'
  } else if(!fromSTAR){
    fastq_dir <- getConfig()[["subdirs"]][["BAM"]]
    orig <- list.files(path=fastq_dir,pattern="\\.bam$")
    todo <- str_replace(orig, '(\\.transcript)?\\.bam$', '')
    extension <- str_match(orig[1], '(\\.transcript)?\\.bam$')[1,1]
  }else{
    fastq_dir <- getConfig()[["subdirs"]][["BAM"]]
    orig <- list.files(path=fastq_dir,pattern="Aligned.toTranscriptome.out.bam$")
    todo <- str_replace(orig, 'Aligned.toTranscriptome.out.bam$', '')
    extension <- str_match(orig[1], 'Aligned.toTranscriptome.out.bam$')[1,1]
  }
  keep <- setdiff(todo, done)
  if(length(keep)==0){
    message("Expression already calculated")
    return()
  }
  reference_genome_path <- file.path(getConfig()[["subdirs"]][["Utils"]],"Reference_Genome")
  reference_genome_name <- file.path(getConfig()[["reference_genome_name"]])
  keep <- outer(keep, extension, paste0)
  ## Write slurm preamble to a shell script
  if(slurm){
    con<-file(file.path(getConfig()[["subdirs"]][["FASTQ"]],"batchSubmitJob.sh"),open = "w")
    ram_requested<-parallel_threads*ram_per_node
    writeLines(c("#!/bin/bash",
                 paste0("#SBATCH -n ",ncores," # Number of cores"),
                 "#SBATCH -N 1 # Ensure that all cores are on one machine",
                 paste0("#SBATCH -t ",days_requested,"-00:00 # Runtime in D-HH:MM"),
                 ifelse(is.null(slurm_partition),"", paste0("#SBATCH -p ",slurm_partition," # Partition to submit to")),
                 paste0("#SBATCH --mem=",ram_requested," # Memory pool for all cores (see also --mem-per-cpu)"),
                 paste0("#SBATCH -o ",file.path(getConfig()[["subdirs"]][["FASTQ"]],"rsem_%a.out")," # File to which STDOUT will be written"),
                 paste0("#SBATCH -e ",file.path(getConfig()[["subdirs"]][["FASTQ"]],"rsem_%a.err"), " # File to which STDERR will be written"),
                 'module add bowtie2'),con=con)
    argumentFile <- file.path(getConfig()[["subdirs"]][["FASTQ"]], "arguments_chunk_${SLURM_ARRAY_TASK_ID}.txt")
  } else{
    #If we aren't using slurm, set the number of chunks to one.
    nchunks<-1
    argumentFile <- file.path(getConfig()[["subdirs"]][["FASTQ"]], "arguments_chunk_1.txt")
  }
  if(!paired){
    #get file names, commandline for single-end
    # keep<-paste0(keep,".fastq")
    myfiles<-file.path(fastq_dir,keep)
    fragLenArg <- ''
    bamArg <- ''
    if(!is.null(frag_mean) && !is.null(frag_sd)){
      fragLenArg <- paste0(" --fragment-length-mean ", frag_mean, " --fragment-length-sd ", frag_sd)
    }
    if(fromBAM){
      command<-paste0("cd ",rsem_dir," && parallel -j ",parallel_threads*bowtie_threads," rsem-calculate-expression --no-bam-output --bam {} ",file.path(reference_genome_path,reference_genome_name)," {/.} :::: ", argumentFile)
    } else{
      command<-paste0("cd ",rsem_dir," && parallel -j ",parallel_threads," rsem-calculate-expression --bowtie2 -p ", bowtie_threads," {} ",file.path(reference_genome_path,reference_genome_name)," {/.} :::: ", argumentFile)
    }
  }else{
    #get file names, commandline for Paired end
    fastq_files<-file.path(fastq_dir,sort(keep))
    if(!is.null(frag_mean) || !is.null(frag_sd)){
      warning('`frag_mean` and `frag_sd` ignored for paired-end reads.')
    }
    if(fromBAM){
      myfiles<-file.path(fastq_dir,keep)
      command<-paste0("cd ",rsem_dir," && parallel -j ",parallel_threads*bowtie_threads," rsem-calculate-expression --no-bam-output --paired-end --bam {} ",file.path(reference_genome_path,reference_genome_name)," {/.} :::: ", argumentFile)
    } else{
      pairs<-matrix(fastq_files,ncol=2,byrow=TRUE)
      myfiles<-cbind(pairs,gsub("\\.fastq","",gsub("_[12]\\.","\\.",basename(pairs[,1]))))
      command<-paste0("cd ",rsem_dir," && parallel -n 3 -j ",parallel_threads," rsem-calculate-expression --bowtie2 -p ", bowtie_threads," --paired-end {1} {2} ",file.path(reference_genome_path,reference_genome_name)," {3.} :::: ", argumentFile)
    }
  }
  ## for both, divide split the files into chunks (maybe just 1) and write the argument chunks
  chunked<-.chunkDataFrame(data.frame(files=myfiles),nchunks)
  for(i in seq_along(chunked)){
    con2=file(file.path(getConfig()[["subdirs"]][["FASTQ"]],paste0("arguments_chunk_",i,".txt")))
    writeLines(t(chunked[[i]]),con=con2)
    close(con2)
  }
  if(slurm){
    ## add the rsem command to the shell script
    message('Using ', nchunks , 'clusters ')
    message("Sending ", command)
    writeLines(paste0(command,"\n"),con=con)
    close(con)
    slurm_command<-paste0("sbatch --array=1-",nchunks," ",file.path(getConfig()[["subdirs"]][["FASTQ"]],"batchSubmitJob.sh"))
    system(slurm_command)
  } else{
    ## or just execute it
    cat(command)
    system(command)
  }
}


#' Assemble an expression matrix of all results
#' 
#' Put all the counts from the individual libraries into a single matrix result
#' 
#' Assemble an expression matrix from the individual libraries.
#' @param force \code{logical} rerun even if output exists
#'@export 
RSEMAssembleExpressionMatrix <- function(force=FALSE){
  cond_eval <- length(list.files(getConfig()[["subdirs"]][["RSEM"]], pattern="rsem_"))<4
  if(cond_eval|force){
    message("Assembling counts matrix")
    rsem_files <- list.files(getConfig()[["subdirs"]][["RSEM"]], pattern="genes.results", full.names = TRUE)
    # Read all files and create a list of data.tables
    rsem_list <- lapply(rsem_files, function(x, ...){
      y<-fread(x, drop=c("length", "FPKM")); 
      y$sample_name <- gsub(".genes.results", "", basename(x));
      return(y);})
    
    # Bind all the data.tables to create a long data.table
    rsem_data_long <- rbindlist(rsem_list)

    ## remove 'Aligned.toTranscriptome.out' from rsem file names
    rsem_data_long[,sample_name := gsub("Aligned.toTranscriptome.out", "", sample_name)]
    
    # Rename a column, so that it's a bit R friendly
    setnames(rsem_data_long, "transcript_id(s)", "transcript_ids")
    
    # Reshape the long data.table to create a matrix
    # Assume that missing entries would have a TPM value of 0
    rsem_tpm_matrix <- dcast.data.table(rsem_data_long, gene_id+transcript_ids~sample_name, value.var = "TPM", fill=0)
    rsem_count_matrix <- dcast.data.table(rsem_data_long, gene_id+transcript_ids~sample_name, value.var = "expected_count", fill=0)
    rsem_effective_length_matrix<-dcast.data.table(rsem_data_long,gene_id+transcript_ids~sample_name,value.var="effective_length",fill=0)
    # Keep track of the gene_id - transcript mapping 
    rsem_txs_table <- rsem_tpm_matrix[, c("gene_id","transcript_ids"), with=FALSE]
    
    # Remove the transcript column (we create an annotation table in the next chunk)
    rsem_tpm_matrix <- rsem_tpm_matrix[, transcript_ids:=NULL][order(gene_id)]
    rsem_count_matrix <- rsem_count_matrix[, transcript_ids:=NULL][order(gene_id)]
    rsem_effective_length_matrix <- rsem_effective_length_matrix[, transcript_ids:=NULL][order(gene_id)]
    
    # Write to disk
    write.csv(rsem_tpm_matrix, file=file.path(getConfig()[["subdirs"]][["RSEM"]],"rsem_tpm_matrix.csv"), row.names=FALSE)
    write.csv(rsem_effective_length_matrix, file=file.path(getConfig()[["subdirs"]][["RSEM"]],"rsem_effective_length_matrix.csv"), row.names=FALSE)
    write.csv(rsem_count_matrix, file=file.path(getConfig()[["subdirs"]][["RSEM"]],"rsem_count_matrix.csv"), row.names=FALSE)
    write.csv(rsem_txs_table,file=file.path(getConfig()[["subdirs"]][["RSEM"]],"rsem_txs_table.csv"))
  }
}

##' Summarize the number of mapped, unmapped and ambiguous reads from RSEM output
##'
##' @param dir Directory containing RSEM cnt files.  If missing, use the RSEM directory from project.
##' @param log return log of reads?
##' @param plot print a ggplot summary
##' @return returns a data table of # reads per library, invisibly
RSEMSummarizeMapping <- function(dir, log=TRUE, plot=TRUE){
    if(missing(dir)) dir <- getConfig()[["subdirs"]][["RSEM"]]
    cnt_files <- list.files(dir, pattern="*.cnt", full.names = TRUE, recursive=TRUE)
    fastqc_list <- lapply(cnt_files, function(x, ...){
        y <- fread(x, skip=3)
        setnames(y, c('nmap', 'reads'))
        y[,curead:=cumsum(reads)]
        file <- str_replace(basename(x), fixed('.cnt'), '')
        total <- y[nmap==Inf,curead]
        mapped <- total-y[nmap==0,reads]
        unique <- max(y[nmap==1,reads], 0)
        z <- data.table(file, total, mapped, unique)
        return(z)})
    fql <- rbindlist(fastqc_list)
    fqlM <- data.table::melt.data.table(fql, id.vars='file')
    if(log){
        fqlM[,value:=log10(value+1)]
    }
    fqlM[,outlier:={
        bs <- boxplot.stats(value)
        value<bs$stats[1] | value > bs$stats[5]
    },keyby='variable'
     ]

    ggp <- ggplot(fqlM, aes(x=variable, y=value))+geom_boxplot() + geom_text(aes(label=ifelse(outlier, file, "")), size=2) + ylab(if(log) 'log10(reads+1)' else 'reads')
    if(plot) print(ggp)
    invisible(fqlM)
}

#' Produce a report listing the tools and packages used by the pipeline.
#' 
#' Produce a report listing the tools and packages used by the pipeline.
#' 
#' Helps with reproducibility by producing a report listing the tools and packages used by the pipeline.
#' This needs to be run on the system that produced the results. Output to project directory OUTPUT/pipeline_report.md.
#' @export
pipelineReport<-function(){
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
  obj<-readRDS(list.files(confdir,pattern="configuration.rds",full.names=TRUE))
  
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
    mat <- fread(list.files(getConfig()[["subdirs"]][["RSEM"]],pattern="rsem_count_matrix.csv",full.names=TRUE))
  else
    mat <- fread(list.files(getConfig()[["subdirs"]][["RSEM"]],pattern="rsem_tpm_matrix.csv",full.names=TRUE))
  featuredata <- fread(list.files(getConfig()[["subdirs"]][["RSEM"]],pattern="rsem_fdata.csv",full.names=TRUE))
  pdata <- fread(list.files(getConfig()[["subdirs"]][["RSEM"]],pattern="rsem_pdata.csv",full.names=TRUE))
  
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
    fastqfiles<-list.files(fastqdir,pattern="fastq$",full.names=TRUE)
  }else{
    fastqfiles<-list.files(fastqdir,pattern="\\.assembled\\.fastq$",full.names=TRUE)
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
  files<-(list.files(path=fastq_dir,pattern="*.fastq",full.names=FALSE))
  f<-file.path(pear_directory,"pear_arguments.txt")
  connection<-file(f,open="w")
  writeLines(files,con = connection)
  close(connection)
  command<-paste0("cd ", fastq_dir, " && parallel -j ",ncores, " -n2 pear -f {1} -r {2} -o ",file.path(pear_directory,"{1}")," :::: < ",file.path(pear_directory,"pear_arguments.txt"))
  system(command)
}


#' Download SRA files from SRX accessions
#' 
#' Download SRA files from SRX accessions. Downloads asynchronously. Won't message you when complete, so can't be run in 
#' batch mode at the moment.
#' @param x \code{character} a vector of SRX accession numbers
#' @export
getDataFromSRX<-function(x=NULL){
  if(is.null(x)){
    stop("Please pass a vector of SRX numbers.")
  }
  sra_con<-getConfig()[["sra_con"]]
  run_accession <- listSRAfile(x, sra_con, fileType = "sra" )$run
  aspera_url <- paste0("anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra", "/", substr(run_accession,1,3), "/", substr(run_accession,1,6), "/", run_accession, "/", run_accession, ".sra")
  out<-paste0('ascp -i ',gsub(" ","\\\\ ",getConfig()[["aspera_path"]]),'/asperaweb_id_dsa.openssh -k 1 -T -l200m ', aspera_url, " ",getConfig()[["subdirs"]][["SRA"]])
  for(i in out){
    system(i,wait=FALSE)
  }
  message("Files are downloading. Wait and check your downloads before proceeding")
}

#'Quick and dirty annotations from SRAmetadb
#'
#'Gets annotations for SRR files based on SRAmetadb contents
#'@param x \code{character} vector of SRX accessions
#'@export
annotationsFromSRX<-function(x){
  samples<-data.table(listSRAfile(x,getConfig()[["sra_con"]]))
  setnames(samples,"sample","sample_accession")
  
  src<-src_sqlite(getConfig()[["sra_con"]]@dbname, create = F)
  sample_table<-tbl(src, sql("SELECT * FROM sample"))
  
  pdata<-merge(data.table(as.data.frame(filter(sample_table,sample_accession%in%samples$sample_accession))),samples,by="sample_accession")
  pdata<-select(pdata,sample_accession,experiment,run,sample_attribute)
  
  attributes<-strsplit(pdata$sample_attribute,"\\|\\|")
  names(attributes)<-pdata$run
  attributes<-lapply(attributes,function(x)t(cbind(x)))
  annotations<-plyr::ldply(lapply(attributes,function(x){
    column_names<-gsub("(^.+?):.*","\\1",x)
    matrix_entries<-gsub(" +$","",gsub("^ +","",gsub("^.+?:(.*)","\\1",x)))
    colnames(matrix_entries)<-column_names
    matrix_entries
  }))
  setnames(annotations,".id","run")
  annotations<-merge(pdata,annotations,by="run")
  annotations[,sample_attribute:=NULL]
  setnames(annotations,"run","srr")
  write.csv(annotations, file=file.path(getConfig()[["subdirs"]][["RSEM"]],"rsem_pdata.csv"), row.names=FALSE) 
  message("Wrote pData to RSEM/rsem_pdata.csv")
  annotations
}

#' Use the Tophat tool to align reads
#' 
#' Uses Tophat to align reads in FASTQ files against the reference genome hg38 from UCSC database. Optionally you can specify paired end reads. The code assumes 
#' paired reads have fastq files that differ by one character (i.e. sampleA_read1.fastq, sampleA_read2.fastq) and will perform
#' matching of paired fastq files based on that assumption using string edit distance. Read 1 is assumed to be upstream
#' and read 2 is assumed to be downstream. 
#' 
#' The number of parallel_threads*tophat_threads should not be more than the number of cores available on your system.
#' @param path \code{character} specifying an \emph{absolute path} path to the iGenome directory.
#' @param parallel_threads \code{integer} specify how many parallel processes to spawn
#' @param tophat_threads \code{integer} specify how many threads bowtie should use.
#' @param paired \code{logical} specify whether you have paried reads or not.
#' @param nchunks \code{integer} number of chunks to split the files for a slurm job. Ignored if slurm = FALSE
#' @param days_requested \code{integer} number of days requested for the job (when submitting a slurm job). Ignored if slurm = FALSE
#' @param slurm \code{logical} if \code{TRUE} job is submitted as a slurm batch job, otherwise it's run on the local machine. Slurm jobs will honour the nchunks and days_requested arguments. 
#' @param slurm_partition \code{character} the slurm partition to submit to. Ignored if slurm=FALSE
#' @param ram_per_node \code{numeric} The number of Mb per node. Ignored if slurm=FALSE. Default of \code{parallel_threads*bowtie_threads*1000}
#' @note The amount of memory requested should be set to bowtie_threads*parallel_threads*1G as this is the default requested by samtools for sorting. If insufficient memory is requested, the bam files will not be created successfully.
#' @export
sequenceAlignmentTopHat = function(path="/shared/silo_researcher/Gottardo_R/jingyuan_working/iGenome/Mus_musculus/UCSC/mm10", parallel_threads=1,tophat_threads=6, paired=FALSE, nchunks=10,days_requested=5,
                                   slurm=FALSE,slurm_partition="gottardo_r",ram_per_node=tophat_threads*parallel_threads*1200)
{
  ncores<-parallel_threads*tophat_threads
  if(ncores>parallel::detectCores()&!slurm){
    stop("The number of parallel_threads*bowite_threads is more than the number of cores detected by detectCores() on the local machine for non-slurm jobs")
  }
  if(!slurm){
    #If we aren't using slurm, set the number of chunks to one.
    nchunks<-1
  }
  #Chunking for slurm
  .chunkDataFrame<-function(df=pairs,nchunks=nchunks){        
    groupsize<-nrow(df)%/%nchunks
    split(as.data.frame(df),gl(nchunks,groupsize,length=nrow(df)))
  }
  
  subdirs<-getConfig()[["subdirs"]]
  subdirs[["iGenome"]]<-path
  assignConfig("subdirs",subdirs)
  tophat_dir <- getConfig()[["subdirs"]][["TOPHAT"]]
  fastq_dir <- getConfig()[["subdirs"]][["FASTQ"]]
  genome.gtf <- file.path(getConfig()[["subdirs"]][["iGenome"]],"genes.gtf")
  genome.index <- file.path(getConfig()[["subdirs"]][["iGenome"]],"genome")
  fastq.files <- list.files(fastq_dir, ".fastq$", full.names = TRUE, recursive = TRUE)
  if(paired){
    samples <- unique(sub("_..fastq", "", basename(fastq.files)))
    f1 <- fastq.files[match(paste0(samples,"_1.fastq"), basename(fastq.files))]
    f2 <- fastq.files[match(paste0(samples,"_2.fastq"), basename(fastq.files))]
    pairs <- cbind(f1, f2, samples)
    chunked<-.chunkDataFrame(pairs,nchunks)
    for(i in seq_along(chunked)){
      con=file(file.path(getConfig()[["subdirs"]][["FASTQ"]],paste0("arguments_chunk_",i,".txt")))
      writeLines(t(chunked[[i]]),con=con)
      close(con)
    }
    if(slurm){
      #With chunking enabled this is meant to be submitted as a slurm batch job with the --array option set to 1-nchunks
      con<-file(file.path(getConfig()[["subdirs"]][["FASTQ"]],"batchSubmitJob.sh"),open = "w")
      command<-paste0("cd ",tophat_dir," && parallel -n 3 -j ",parallel_threads," tophat -p ", tophat_threads," -G ", genome.gtf, " --keep-fasta-order --no-novel-junc -o {3} ",
                      genome.index, " {1} {2} :::: ",file.path(getConfig()[["subdirs"]][["FASTQ"]],paste0("arguments_chunk_${SLURM_ARRAY_TASK_ID}.txt","\n")))
      #Before writing the command we need to  preamble for the slurm job
      # We'll request 1Gb of RAM per parallel thread per node
      ram_requested<-parallel_threads*ram_per_node
      writeLines(c("#!/bin/bash\n",
                   paste0("#SBATCH -n ",ncores,"                 # Number of cores"),
                   "#SBATCH -N 1                 # Ensure that all cores are on one machine",
                   paste0("#SBATCH -t ",days_requested,"-00:00           # Runtime in D-HH:MM"),
                   paste0("#SBATCH -p ",slurm_partition,"    # Partition to submit to"),
                   paste0("#SBATCH --mem=",ram_requested,"            # Memory pool for all cores (see also --mem-per-cpu)"),
                   paste0("#SBATCH -o ",file.path(getConfig()[["subdirs"]][["FASTQ"]],"tophat_%a.out"),"      # File to which STDOUT will be written"),
                   paste0("#SBATCH -e ",file.path(getConfig()[["subdirs"]][["FASTQ"]],"tophat_%a.err"),"      # File to which STDERR will be written")),con=con)
      
      writeLines(paste0(command,"\n"),con=con)      
      close(con)
      #Now the command to launch the job is sbatch batchSubmitJob.sh
      command<-paste0("sbatch --array=1-",nchunks," ",file.path(getConfig()[["subdirs"]][["FASTQ"]],"batchSubmitJob.sh"))
    }else{
      # Command for a non-slurm job
      command<-paste0("cd ",tophat_dir," && parallel -n 3 -j ",parallel_threads," tophat -p ", tophat_threads," -G ", genome.gtf, " --keep-fasta-order --no-novel-junc -o {3} ", 
                      genome.index, " {1} {2} :::: ",file.path(getConfig()[["subdirs"]][["FASTQ"]],paste0("arguments_chunk_1.txt","\n")))        
    }
  }else{
    myfiles <- cbind(fastq.files, gsub(".fastq", "", basename(fastq.files)))
    chunked<-.chunkDataFrame(myfiles,nchunks)
    for(i in seq_along(chunked)){
      con=file(file.path(getConfig()[["subdirs"]][["FASTQ"]],paste0("arguments_chunk_",i,".txt")))
      writeLines(t(chunked[[i]]),con=con)
      close(con)
    }
    if(slurm){
      #With chunking enabled this is meant to be submitted as a slurm batch job with the --array option set to 1-nchunks
      con<-file(file.path(getConfig()[["subdirs"]][["FASTQ"]],"batchSubmitJob.sh"),open = "w")
      command<-paste0("cd ",tophat_dir," && parallel -n 2 -j ",parallel_threads," tophat -p ", tophat_threads," -G ", genome.gtf, " --keep-fasta-order --no-novel-junc -o {2} ",
                      genome.index, " {1} :::: ",file.path(getConfig()[["subdirs"]][["FASTQ"]],paste0("arguments_chunk_${SLURM_ARRAY_TASK_ID}.txt","\n")))
      #Before writing the command we need to  preamble for the slurm job
      # We'll request 1Gb of RAM per parallel thread per node
      ram_requested<-parallel_threads*ram_per_node
      writeLines(c("#!/bin/bash\n",
                   paste0("#SBATCH -n ",ncores,"                 # Number of cores"),
                   "#SBATCH -N 1                 # Ensure that all cores are on one machine",
                   paste0("#SBATCH -t ",days_requested,"-00:00           # Runtime in D-HH:MM"),
                   paste0("#SBATCH -p ",slurm_partition,"    # Partition to submit to"),
                   paste0("#SBATCH --mem=",ram_requested,"            # Memory pool for all cores (see also --mem-per-cpu)"),
                   paste0("#SBATCH -o ",file.path(getConfig()[["subdirs"]][["FASTQ"]],"tophat_%a.out"),"      # File to which STDOUT will be written"),
                   paste0("#SBATCH -e ",file.path(getConfig()[["subdirs"]][["FASTQ"]],"tophat_%a.err"),"      # File to which STDERR will be written")),con=con)
      
      writeLines(paste0(command,"\n"),con=con)      
      close(con)
      #Now the command to launch the job is sbatch batchSubmitJob.sh
      command<-paste0("sbatch --array=1-",nchunks," ",file.path(getConfig()[["subdirs"]][["FASTQ"]],"batchSubmitJob.sh"))
    }else{
      # Command for a non-slurm job
      command<-paste0("cd ",tophat_dir," && parallel -n 2 -j ",parallel_threads," tophat -p ", tophat_threads," -G ", genome.gtf, " --keep-fasta-order --no-novel-junc -o {2} ", 
                      genome.index, " {1} :::: ",file.path(getConfig()[["subdirs"]][["FASTQ"]],paste0("arguments_chunk_1.txt","\n")))        
    }
  }
  cat(command)
  system(command)
}

#' Create Sample file required as the input of RNASeQC
#' 
#' The sample file is the tab-delimited description of samples and their bams. The header of the file should be: SampleID BamFile Notes
#' The BamFile should be the path to the input file
.createSampleFile = function(){
  bam.files <- list.files(getConfig()[["subdirs"]][["TOPHAT"]], "accepted_hits_added.bam$", full.names = TRUE, recursive = TRUE)
  sampleid <- sapply(strsplit(bam.files, "/"), function(x) x[length(x)-1])
  note <- sapply(strsplit(sampleid, "_"), function(x) x[1])  
  dat <- data.frame("Sample ID" = sampleid, "Bam File" = bam.files, "Notes" = note, stringsAsFactors = F)
  write.table(dat, file = paste0(getConfig()[["subdirs"]][["RAWANNOTATIONS"]],"/samplefile"), col.names=TRUE, row.names=FALSE, quote=FALSE, sep = "\t")
  dat
}

#' Create a dict (dictionary) file for the refrence genome required as the input of RNASeQC
#' 
#' This dict file should be within the same folder as the refrence genome
#' @param fasta_file \code{character} the name of the fasta file, must be specified
#' @param picard.path \code{character} specifying an \emph{absolute path} path to the PICARD software directory.
.createGenomeDict=function(picard.path="/home/jdeng/bin/picard-tools-1.110/")
{
  genome.fa <- file.path(getConfig()[["subdirs"]][["iGenome"]], "genome.fa")
  cmd=paste("java -Xmx2048m -Xms2048m -jar ",picard.path,"/CreateSequenceDictionary.jar REFERENCE=",genome.fa," OUTPUT=",gsub(".fa.*",".dict",genome.fa),sep="")
  print(cmd)
  system(cmd)
}

#' Run RNASeQC quality control
#' 
#' Pre-run Checklist
#' 1. Are the contig names consistent betwen BAM, Refrence and GTF file? 
#' 2. Is the BAM indexed? we index the bam files in the prepareBamFiles
#' 3. Is you reference Indexed? We index the refrence genome before runing the tophat
#' 4. Dose the reference genome have a dict file? we generated the dict file in the createGenomeDict
#' 
#' @param rna_seqc.path \code{character} specifying an \emph{absolute path} path to the RNASeQC directory.
#' @param picard.path \code{character} specifying an \emph{absolute path} path to the PICARD software directory.
#' @param paired \code{logical} specify whether data is paired. Defaults to FALSE.
#' @param ncores \code{integer} number of cores for running in parallel
#' @export
qualityRNASeQC = function(picard.path="/home/jdeng/bin/picard-tools-1.110/", ncores=6, rna_seqc.path="/home/jdeng/bin/RNA-SeQC_v1.1.7.jar", paired=FALSE){
  print("Step 1: Prepare bam files")
  genome.fa <- file.path(getConfig()[["subdirs"]][["iGenome"]], "genome.fa")
  bam.files <- list.files(getConfig()[["subdirs"]][["TOPHAT"]], "accepted_hits.bam$", full.names = TRUE, recursive = TRUE)
  mclapply(bam.files, function(s){
    # add read group to bam
    cmd <- paste0("nohup java -Xmx2048m -Xms2048m -jar ", picard.path, "AddOrReplaceReadGroups.jar I=", s, " O=", sub(".bam", "_added.bam", s), " SORT_ORDER=coordinate RGID=1 RGLB=1 RGPL=illumina RGPU=barcode RGSM=1")
    system(cmd)
    # index bam
    cmd <- paste0("samtools index ", sub(".bam", "_added.bam", s))
    system(cmd)
  }, mc.preschedule = F, mc.cores = ncores)
  
  print("Step2: Create Sample file\n")
  .createSampleFile()
  
  print("step3: Greate Genome dict")
  .createGenomeDict(picard.path)
  
  print("step4: Run RNASeQC")
  rnaseqc.dir <- getConfig()[["subdirs"]][["RNASEQC"]]
  genome.gtf <- file.path(getConfig()[["subdirs"]][["iGenome"]],"genes.gtf")
  samplefile <- file.path(getConfig()[["subdirs"]][["RAWANNOTATIONS"]],"samplefile")
  if(paired){
    cmd <- paste0("nohup java -Xmx2048m -Xms2048m -jar ", rna_seqc.path, " -o ", rnaseqc.dir, " -r ", genome.fa, " -s ", samplefile, " -t ", genome.gtf)                        
  }else{
    cmd <- paste0("nohup java -Xmx2048m -Xms2048m -jar ", rna_seqc.path, " -o ", rnaseqc.dir, " -r ", genome.fa, " -s ", samplefile, " -singleEnd -t ", genome.gtf)
  }
  system(cmd)
}


#' Generate QualityControl Matric.
#'
#' calculate qc statistics: nGenesOn, librarySize, alignment percentage to reference genome, exon mapping rate,
#' and per base sequence quality output by FASTQC software
#'
#' @param paired \code{logical} specify whether you have paried reads or not.
#' @export
QualityControl <- function(paired=FALSE){
  countMatrix <- fread(list.files(getConfig()[["subdirs"]][["RSEM"]],pattern="rsem_count_matrix.csv",full.names=TRUE))[ ,-1,with=FALSE]
  sample <- gsub("Aligned.toTranscriptome.out", "", colnames(countMatrix))
  result <- data.table(Sample=sample, nGeneOn = colSums(countMatrix>5))
  ### alignment rate
  bam.dir <- getConfig()[["subdirs"]][["BAM"]]
  rsem.dir <- getConfig()[["subdirs"]][["RSEM"]]
  starAlignSummary <- list.files(bam.dir, "Log.final.out$", full.names=TRUE, recursive=FALSE)
  rsemAlignSummary <- list.files(rsem.dir, ".cnt$", full.names=TRUE, recursive=TRUE)
  rsemFileNames <- gsub("Aligned.toTranscriptome.out.cnt", "", basename(rsemAlignSummary))
  starFileNames <- gsub("Log.final.out", "", basename(starAlignSummary))
  rsemAlignSummary <- rsemAlignSummary[match(starFileNames, rsemFileNames)]
  res <- t(sapply(seq(starAlignSummary), function(i) {
    print(i)
    starIn <- read.delim(starAlignSummary[i], stringsAsFactors=FALSE)
    libSize <- as.numeric(starIn[4,2])
    alignRate <- (as.numeric(sub("%","",starIn[8,2]))+as.numeric(sub("%","",starIn[23,2])))/100
    rsemIn <- scan(rsemAlignSummary[i], what = "", nlines=1)
    exonRate <- as.numeric(rsemIn[2])/libSize
    cbind(libSize, alignRate, exonRate)}))
  
  res <- data.table(Sample=starFileNames, libSize=res[,1], alignRate=res[,2], exonRate=res[,3])
  result <- merge(result, res, by="Sample")
  ###### Fastqc
  fastqc <- list.files(getConfig()[["subdirs"]][["FASTQC"]], "summary.txt", full.names=TRUE, recursive=TRUE)
  
  ## remove R[12]...fastqc suffix from file name. Initial (.*) causes right hand semantics so non-greedy matches and
  ## captures the fle name without the suffix.
  ## Extra complicated regexp just in case someone has inserted an R1 or R2 in file name before the read pair identifier
  ## added ? after R as it is possible for fastq files to not have the R. (? matches 0 or 1)
  fastqcFileNames <- sub("(.*)_R?[12].*_fastqc$", "\\1", sapply(strsplit(fastqc, "/"), function(x) x[length(x)-1]), perl=TRUE)

   if(!paired){
    res_fastqc <- sapply(fastqc, function(i) read.delim(i, header=FALSE, sep="\t", stringsAsFactors=FALSE)[2,1])
    res_fastqc[res_fastqc == "WARN"] <- "PASS"
    res_fastqc <- data.table(Sample=fastqcFileNames, fastqc=res_fastqc)
    res_fastqc$Sample <- gsub("_fastqc", "", res_fastqc$Sample)
    result <- merge(result, res_fastqc, by="Sample")
  }else{
    res_fastqc <- sapply(fastqc, function(i) read.delim(i, header=FALSE, sep="\t", stringsAsFactors=FALSE)[2,1])
    res_fastqc[res_fastqc == "WARN"] <- "PASS"
    res_fastqc <- data.table(Sample=fastqcFileNames, fastqc=res_fastqc)
    res_fastqc <- res_fastqc[,.SD[,paste(fastqc, collapse="_"), by=Sample]]
    res_fastqc$Sample <- gsub("_fastqc", "", res_fastqc$Sample)
    result <- merge(result, res_fastqc, by="Sample")
  }
  write.csv(result, file=file.path(getConfig()[["subdirs"]][["OUTPUT"]],"quality_control_matrix.csv"), row.names=FALSE)
  result
}



#' Build a genome index for STAR
#'
#' Builds a genome Index at `Utils/Reference_Genome/STARIndex`
#'
#' You must specify the Utils path if it is not already defined, and have your genome and annotation file in a folder titled
#' `Reference_Genome`. This function will construct the genome index using STAR tools.
#' The command line is the default shown in the documentation.
#' A fasta_file must be provided.
#' gtf_file is highly recommended, and if not provided, STAR will run without annotations.
#' @param path \code{character} specifying an \emph{absolute path} path to the Utils directory.
#' @param gtf_file \code{character} the name of the gtf file.
#' @param fasta_file \code{character} the name of the fasta file, must be specified
#' @param star_threads \code{integer} specify how many threads star should use.
#' @param name \code{character} the name of the genome output
#' @export

buildGenomeIndexSTAR = function (path = NULL, gtf_file = "", fasta_file = NULL, star_threads = 1, name=NULL){
  if (is.null(fasta_file)){
    stop("You must provide a fasta_file")
  }
  if (is.null(path) & inherits(try(getConfig()[["subdirs"]][["Utils"]], silent = TRUE), "try-error")) {
    stop("Please specify where to build the reference genome")
  }
  else if (is.null(path)) {
    path <- file.path(getConfig()[["subdirs"]][["Utils"]], "Reference_Genome/STARIndex")
  }
  else {
    subdirs <- getConfig()[["subdirs"]]
    subdirs[["Utils"]] <- path
    assignConfig("subdirs", subdirs)
    path <- file.path(getConfig()[["subdirs"]][["Utils"]], "Reference_Genome/STARIndex")
  }
  if (substr(path, 1, 1) != "/") {
    stop("'path' must be an absolute path")
  }
  dir.create(path, showWarnings=FALSE)
  if (gtf_file == "") {
    gtfopt <- ""
  }
  else {
    gtfopt <- "--sjdbGTFfile"
    gtf.file <- file.path(paste0(getConfig()[["subdirs"]][["Utils"]], "/Reference_Genome/", gtf_file))
  }
  fasta.file <- file.path(paste0(getConfig()[["subdirs"]][["Utils"]], "/Reference_Genome/", fasta_file))
  if (length(fasta.file) > 1) {
    fasta.file <- paste(fasta.file, collapse = ",")
  }
  if (length(list.files(pattern="SAindex", path=path)) == 0) {
    command = paste0("STAR --runThreadN ", star_threads, " --runMode genomeGenerate --genomeDir ", path, " --genomeFastaFiles ", fasta.file, " ", gtfopt, " ", gtf.file)
    system(command)
  }
  else {
    message("Reference Genome Found")
  }
  assignConfig("reference_genome_name", name)
}



#' Use the STAR tool to align the reads 
#'
#' Use the STAR tool to align reads and generate transcriptome bam files
#'
#' Uses STAR to align reads in FASTQ files against the reference genome. Optionally you can specify paired end reads. The code assumes
#' paired reads have fastq files that differ by one character (i.e. sampleA_read1.fastq, sampleA_read2.fastq) and will perform
#' matching of paired fastq files based on that assumption using string edit distance. Read 1 is assumed to be upstream
#' and read 2 is assumed to be downstream.
#' The number of parallel_threads*star_threads should not be more than the number of cores available on your system.
#' @param parallel_threads \code{integer} specify how many parallel processes to spawn
#' @param star_threads \code{integer} specify how many threads star should use.
#' @param paired \code{logical} specify whether you have paried reads or not.
#' @param paired_pattern \code{character} specify the suffix of the paried-end fastq file names.
#' @export
AlignmentSTAR <- function(parallel_threads=1, star_threads=6, paired=TRUE, paired_pattern=c("_1.fastq", "_2.fastq")){
  
  ncores<-parallel_threads*star_threads
  if(ncores>parallel::detectCores()){
    stop("The number of parallel_threads*bowite_threads is more than the number of cores detected by detectCores() on the local machine for non-slurm jobs")
  }
  .chunkDataFrame<-function(df=pairs){
    groupsize<-nrow(df)
    split(as.data.frame(df),gl(1,groupsize,length=nrow(df)))
  }
  star_dir <- getConfig()[["subdirs"]][["BAM"]]
  done <- str_replace(list.files(path = star_dir, pattern = "\\Aligned.toTranscriptome.out.bam$"), "Aligned.toTranscriptome.out.bam", "")
  fastq_dir <- getConfig()[["subdirs"]][["FASTQ"]]
  if(paired){
    todo <- unique(str_replace(list.files(path=fastq_dir,pattern="\\.fastq$"), paired_pattern, ''))
  }else{
    todo <- unique(str_replace(list.files(path=fastq_dir,pattern="\\.fastq$"), "\\.fastq$", ''))
  }
  extension <- if(paired) paired_pattern else '.fastq'
  keep <- setdiff(todo, done)
  if(length(keep)==0){
    message("Expression already calculated")
    return()
  }
  index <- file.path(getConfig()[["subdirs"]][["Utils"]], "Reference_Genome/")
  keep <- outer(keep, extension, paste0)
  argumentFile <- file.path(getConfig()[["subdirs"]][["FASTQ"]], "arguments_chunk_1.txt")
  if(!paired){
    #get file names, commandline for single-end
    myfiles<-cbind(file.path(fastq_dir,keep), gsub(".fastq", "",keep))
    command<-paste0("cd ",star_dir," && parallel -n 2 -j ",parallel_threads," STAR --runThreadN ", star_threads, " --outFileNamePrefix {2} --genomeDir ", index, " --readFilesIn {1} --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM :::: ", argumentFile)
  }else{
    #get file names, commandline for Paired end
    fastq_files<-file.path(fastq_dir,sort(keep))
    pairs<-matrix(fastq_files,ncol=2,byrow=TRUE)
    myfiles<-cbind(pairs,gsub(paired_pattern[1], "", basename(pairs[,1])))
    command<-paste0("cd ",star_dir," && parallel -n 3 -j ",parallel_threads," STAR --runThreadN ", star_threads, " --outFileNamePrefix {3} --genomeDir ", index, " --readFilesIn {1} {2} --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM :::: ", argumentFile)
  }
  ## for both, divide split the files into chunks (maybe just 1) and write the argument chunks
  chunked<-.chunkDataFrame(data.frame(files=myfiles))
  con2=file(file.path(getConfig()[["subdirs"]][["FASTQ"]], "arguments_chunk_1.txt"))
  writeLines(t(chunked[[1]]),con=con2)
  close(con2)
  cat(command)
  system(command)
}





#' Sum all read counts for multiple instances of a single gene symbol.
#' 
#' A gene symbol can map to multiple UCSC gene cluster ids so in those
#' instances all read counts mapped to a single gene symbol are summed.
#' The summed read counts are returned in an Expression Set along with
#' the feature data table and optionally an experimental design table.
#' @param counts \code{matrix} table of read counts or tpms. Rows are genes and columns are
#' samples
#' @param fDat \code{data.frame} feature data - rows map to counts rows
#' @param cDat \code{data.frame}. experimental design or phenotype data - rows map to counts columns
#' @return ExpressionSet containing the summed counts data (counts), the feature data which will map to the summed counts data (fDat), and optionally, the experimental data (cDat)
#' @export
sumDuplicates <- function(counts, fDat, cDat=NULL) {

    ## apply basic QC
    if (!is.numeric(as.matrix(counts))) 
        stop("`counts` must be numeric")
    if(!is.null(cDat)) {
        if (ncol(counts) != nrow(cDat)) 
            stop("`cDat` must contain as many rows as `counts` columns")
    }
    if (nrow(counts) != nrow(fDat)) 
        stop("`fDat` must contain as many rows as `counts` rows")
 
    ## sum read counts across identical gene symbols by sample
    sumGenes <- function(x, indy=fInd) {
         sums <- tapply(x, indy, sum, na.rm=TRUE)
    }

    ## extract index of gene symbols
    fInd <- as.character(fDat$gene_symbol)
    
    ## sum read counts of identical gene symbols for each sample. Orders
    ## by gene symbol
    newDat <- apply(counts, 2, sumGenes)

    ## remove duplicate gene symbols
    newfDat <- fDat[!duplicated(fDat$gene_symbol),]
    
    ## order feature data by gene symbol
    newfDat <- newfDat[order(newfDat$gene_symbol),]
    rownames(newfDat) <- newfDat$gene_symbol
   
    ## check that order is correct for both data objects
    if(sum(rownames(newDat) != newfDat$gene_symbol) != 0) {
        stop("sumDuplicates: feature and assay data not in correct order")
    }

    ## construct ExpressionSet for return object
    if(!is.null(cDat)) {
        newEset <- ExpressionSet(assayData=newDat,
                                 phenoData=AnnotatedDataFrame(cDat),
                                 featureData=AnnotatedDataFrame(newfDat))
    } else {
        newEset <- ExpressionSet(assayData=newDat,
                                 featureData=AnnotatedDataFrame(newfDat))
    }
   
    return(newEset)
    
} ## end sumDuplicates



#' Build a reference genome
#' 
#' Builds a reference genome at `path/Reference_Genome/`
#' 
#' You must specify the Utils path if it is not already defined, and have your genome in a folder titled
#' `Reference_Genome`. This function will construct the reference genome using RSEM tools.
#' The command line is the default shown in the documentation.
#' `rsem-prepare-reference --gtf gtf_file --transcript-to-gene-map knownIsoforms.txt --bowtie2 fasta_file name`
#' or if doSTAR=TRUE
#'  `rsem-prepare-reference --gtf gtf_file --transcript-to-gene-map knownIsoforms.txt --star --star-path starPath fasta_file name`
#' If the gtf_file is not give, then the transcript-to-gene-map option is not used either. A fasta_file and a name must be provided.
#' @param path \code{character} specifying an \emph{absolute path} path to the Reference_Genome directory.
#' @param gtf_file \code{character} the name of the gtf file. Empty by default. If specified the function will look for the named file in the 'path' directory
#' @param fasta_file \code{character} the name of the fasta file, must be specified and located in the 'path' directory
#' @param isoformsFile \code{character} the name if the known isoforms file (optional). If specified it must be located in the 'path' file.
#' @param name \code{character} the prefix of the genome output.
#' @param doSTAR \code{logical} logical showing whether to build genome with STAR index.
#' @param starPath \code{character}  full path to the 'star' executable.
#' @param threads \code{integer} number of threads to use
#' @export
buildReference <- function(path=NULL, gtf_file=NULL, fasta_file=NULL, isoformsFile=NULL, name=NULL, doSTAR=TRUE, starPath=NULL, threads=1){
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

  ## check that gtf file is provided with known isoforms file
  if(is.null(gtf_file) & !is.null(isoformsFile)) {
     stop("Must provide a .gtf file if knownIsoforms file is provided")
  }

  ## set up gtf and known isoforms arguments if they exist
  gtfopt <- ""
  isoformsOpt <- ""

  if(!is.null(gtf_file)) {
      gtfopt <- "--gtf"
  }

  if(!is.null(isoformsFile)) {
      isoformsOpt<-"--transcript-to-gene-map"
  }

  if(length(fasta_file)>1){
      fasta_file <- paste(fasta_file,collapse=",")
  }

  ## build genome with STAR
  if(doSTAR) {

      ## get path to STAR if not provided
      if(is.null(starPath)) {
          starPath <- dirname(Sys.which("STAR"))
      }
      
      ## check if the reference genome has already been built
      if(length(list.files(pattern="SAindex", path=file.path(getConfig()[["subdirs"]]["Utils"],"Reference_Genome")))==0){
          command = paste("cd ",path," && rsem-prepare-reference ", gtfopt , gtf_file, isoformsOpt, isoformsFile, "--p", threads, "--star --star-path", starPath, fasta_file, name)
          print(command)
          system(command)
      }else{
          message("Reference Genome Found")
      }

  ## bulid genome with RSEM
  } else { 
      
      if(length(list.files(pattern=paste0(name,".chrlist"),path=file.path(getConfig()[["subdirs"]]["Utils"],"Reference_Genome")))==0){
          command = paste0("cd ",path," && rsem-prepare-reference ",gtfopt," ",gtf_file," ",isoformsOpt, " ",isoformsFile, " --bowtie2 ",fasta_file," ",name," ")
          system(command)
      }else{
          message("Reference Genome Found")
      }
  } ## end if doSTAR
      
  #set the reference genome name
  assignConfig("reference_genome_name",name)

} ## end BuildReference



#' Annotate features: map UCSC gene cluster id to gene symbols
#'
#' Annotate features by mapping UCSC gene cluster ids to gene symbols.
#' 
#' @param genome /code{character} specify the genome and version of annotation
#' @param force /code{logical} if TRUE force the annoation even if the feature file already exists
#' @export
annotateUCSC <- function(genome="hg38", force=TRUE) {

    ## read in UCSC annotation file mapping cluster id to gene symbol
    if(genome == "hg38") {
        ##genFile <- paste0("ucsc", genome, "Table.txt")
        ##genFile <- sysdata.rda
        ##fpath <- system.file
        ##fpath <- system.file("extdata", genFile, package="RNASeqPipelineR")
        ##annoFile <- read.table("/shared/silo_researcher/Gottardo_R/cmurie_working/docs/ucsc38Table.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=NULL)
        ##annoFile <- read.table(fpath, sep="\t", stringsAsFactors=FALSE, header=TRUE, row.names=NULL)
        annoFile <- ucschg38Table 
    } else {
        stop(cat("Genome", genome, "is not supported."))
    }
    
   featuredata_outfile <- "rsem_fdata.csv"
    if (!force & file.exists(file.path(getConfig()[["subdirs"]][["RSEM"]], 
        featuredata_outfile))) {
        message("Annotation already done, skipping. Use force=TRUE to rerun.")
        invisible(return(0))
    }
    
    message("Annotating transcripts.")

    ## get experimental annotation file which maps cluster ids to count rows
    rsem_txs_table <- fread(file.path(getConfig()[["subdirs"]][["RSEM"]], 
        "rsem_txs_table.csv"))

    colnames(annoFile)[4] <- "gene_id"

    ## map UCSC cluster ids from rsem to gene symbols in anno file
    joiny <- plyr::join(rsem_txs_table, annoFile, by="gene_id", type="left",
                        match="first")
    
    ## remove irrelevant columns
    joiny$V1 <- joiny$hg38.knownGene.name <- joiny$hg38.knownGene.alignID <- NULL

    ## make human friendly column labels
    colnames(joiny) <- c("gene_clusterID", "transcript_ids","gene_symbol")
    
    message("Writing feature data to rsem_fdata.csv")
    write.csv(joiny, file = file.path(getConfig()[["subdirs"]][["RSEM"]], featuredata_outfile), row.names = FALSE)
    
} ## end annotateUCSC


#' Map the annotation file to the count and tpm data
#'
#' Match the count, tpm, and feature data to the annotation file. Remove
#' samples in count and tpm that are not found in the annotation file and
#' ensure that the count and tpm columns match the rows of the annoation
#' file. Also attach the results of the quality_control table to the
#' annotation file.
#'
#' @param annotationMatch \code{character} Column name of annotation file
#' that maps to the column names of the count and tpm data sets (fastq
#' file names with the suffixes removed including paired end identifiers
#' if they exist). Default column name used is 'Sample'.
#' @return \code{list} list containing the count, tpm, feature, and
#' annotation data. list names are counts, tpms, featureData, annoData.
#'
#' The annotation file, which contains the experimental design information,
#' must be constructed by the user and placed in the RAW_ANNOTATIONS
#' directory. It must be a .csv file and be the only .csv file in the
#' RAW_ANNOTATION directory. By default 'mergeData' will look for a column
#' called 'Sample' that maps each row of the annotation file to the columns
#' of the counts and tpm columns.
#' The quality control matrix is generated by running 'runFASTQ' and then
#' 'QualityControl'
#' @export
mergeData <- function(annotationMatch=NULL) {

    ## read feature, counts, and tpm data files
    fd <- fread(getConfigFile("RSEM", "rsem_fdata.csv"))
    count <-  as.data.frame(fread(getConfigFile("RSEM", "rsem_count_matrix.csv")))[fd$gene_clusterID, ]
    tpm <-  as.data.frame(fread(getConfigFile("RSEM", "rsem_tpm_matrix.csv")))[fd$gene_clusterID, ]
    
    ## double check everyting is ordered correctly
    if(sum(fd$gene_clusterID != count[,1]) != 0 | sum(fd$gene_clusterID != tpm[,1]) != 0) {
        stop("count or tpm matrices are not ordered the same as feature data matrix")
    }

    ## remove gene_id columns
    count <- count[,-1]
    tpm <- tpm[,-1]

    fd <- data.frame(fd)
    rownames(fd) <- rownames(count) <- rownames(tpm) <- fd$gene_ClusterID
 
    anno <- fread(list.files(getConfigDir("RAWANNOTATIONS"), pattern="*.csv$", full.names=TRUE))
    qc <- fread(getConfigFile("OUTPUT", "quality_control_matrix.csv"))

    ## rename column that maps to fastq file prefixes
    if(!is.null(annotationMatch)) {
         setnames(anno, annotationMatch, "Sample")
    }
 
    ## remove samples with descrepancies between fastq and anno
    inter <- intersect(anno$Sample, qc$Sample)
    anno <- anno[Sample %in% inter,]
    qc <- qc[Sample %in% inter,]
    count <- count[,colnames(count) %in% inter]
    tpm <- tpm[,colnames(tpm) %in% inter]
    samples <- colnames(count)

    ## merge
    anno <- merge(anno, qc, by="Sample")
    anno <- anno[match(samples, Sample)]
    cd <- data.frame(wellKey=samples, anno)
    rownames(cd) <- cd$wellKey
    
    ## match sample order with annotation file
    count <- count[,cd$wellKey]
    tpm <- tpm[,cd$wellKey]
    if(sum(colnames(count) == cd$wellKey) != dim(count)[[2]]) {stop("ERROR: colnames don't match")}
    if(sum(colnames(tpm) == cd$wellKey) != dim(tpm)[[2]]) {stop("ERROR: colnames don't match")}

    ## return the updated data
    return(list(counts=count, tpms=tpm, featureData=fd, annoData=cd))

} ## end mergeData
