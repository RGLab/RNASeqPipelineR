% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RNASeqPipelineR.R
\name{runFastQC}
\alias{runFastQC}
\title{Perform FastQC quality control}
\usage{
runFastQC(ncores = 4)
}
\arguments{
\item{ncores}{\code{integer} how many threads to use. Must be less than or equal to number
of cores on machine.}
}
\description{
Runs fastqc quality control on the FASTQ files. Parallel processing is run using fastqc
internal threading (-t option). Each thread requires 250MB of memory.
}
\details{
Produces FASTQC reports for each FASTQ file using fastqc. FASTQ input files must have 
extension \code{fastqc} or \code{fq} (capitalization sensitive). FASTQC output will be
written to FASTQC directory.
}
