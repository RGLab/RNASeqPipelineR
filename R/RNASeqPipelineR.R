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
#' @concept Raphael Gottardo
#' @import data.table
#' @import GEOquery
#' @import SRAdb
NULL

#' Create a new RNASeqPipeline project
#' 
#' Create the skeleton for a new RNASeqPipeline project.
#' 
#' createProject will create a new RNASeqPipeline 
#' project under directory 'name' in the path specified by 'path'.
#' The function creates the directory structure required by RNASeqPipeline
#' within the new project directory, including locations for fastQ files, 
#' fastQC output, RSEM quantification, and optionally GEO and SRA files if the data
#' must be downloaded or linked to a GEO accession. Standard locations will also be provided 
#' to annotate the data utilizing the Immport schema.
#' 
#' @param name \code{character} The name of the project directory
#' @param path \code{character} The path under which to construct the project
#' @param verbose \code{logical} Should verbose output be given?
#' @param load_from_immport \code{logical} Creates a 'Tab' directory for Immport tables if the data are loaded from Immport.
#' @return NULL
#' @export
#' @examples
#' # construct a projects skeleton in a new folder titled "myproject".
#' createProject("myproject",path=".")
#' createProject("myproject",path=".",verbose=TRUE)
#' createProject("myproject",path=".",load_from_immport=TRUE)
createProject <- function(project_name,path=".",verbose=FALSE, load_from_immport=FALSE){
  project_dir <- file.path(path,project_name)
  cmnd_prefix <- "mkdir -p "
  dirs <- c(SRA="SRA",FASTQ="FASTQ",RSEM="RSEM",FASTQC="FASTQC",GEO="GEO")
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
  immport_tables<-vector('list',length(immport_files))
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
getSRAdb <- function(path=NULL){
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
  system(paste0("mkdir -p ",path))  
  
  #download SRA db if necessary
  if(!file.exists(file.path(getConfig()[["subdirs"]]["Utils"],'SRAmetadb.sqlite'))) {
    sqlfile <- getSRAdbFile(destdir = getConfig()[["subdirs"]]["Utils"])
  }
  
  #connect and grab the data
  sra_con <- dbConnect(SQLite(),file.path(getConfig()[["subdirs"]]["Utils"],'SRAmetadb.sqlite'))
  sra_tables <- dbListTables(sra_con)
  assignConfig("sra_tables",sra_tables)
  assignConfig("sra_con",sra_con)
}

#' Detect and configure the aspera client location
#' 
#' Detect the installation location of aspera client based on expected defaults
#' 
#' Sets the location of aspera client based on expected defaults. Or pass the installation path
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
downloadFastQ <- function(){
  cond_eval <- length(list.files("./SRA/", pattern=".sra"))==0
  
  if(!cond_eval){
    GSM_table<-getConfig()[["immport_tables"]][["GSM_table"]]
    for(file in GSM_table[,GSM]) {
      gd <- getGEO(file, destdir=getConfig()[["subdirs"]][["GEO"]])
      SRX_number <- gsub(".*=SRX", "SRX", gd@header$relation[1])
      
      # Convert to aspera address
      sra_con<-getConfig()[["sra_con"]]
      run_accession <- listSRAfile(SRX_number, sra_con, fileType = "sra" )$run
      aspera_url <- paste0("anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra", "/", substr(run_accession,1,3), "/", substr(run_accession,1,6), "/", run_accession, "/", run_accession, ".sra")
      
      system(paste0('ascp -i ',gsub(" ","\\\\ ",getConfig()[["aspera_path"]]),'/asperaweb_id_dsa.openssh -k 1 -T -l200m ', aspera_url, " ",getConfig()[["subdirs"]][["SRA"]]))
    }
  }
}
