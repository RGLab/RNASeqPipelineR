RNASeqPipelineR
===============

Streamline the processing of RNASeq data

## System requirements
You need the following R packages:
data.table
GEOquery
RSQLite
SRAdb

The following command line utilities:  
SRA Toolkit (from NCBI http://www.ncbi.nlm.nih.gov/books/NBK158900/)  
ascp (Aspera scp client, distributed with Aspera Connect)  
RSEM (http://deweylab.biostat.wisc.edu/rsem/)  
bowtie2 (http://www.nature.com/nmeth/journal/v9/n4/full/nmeth.1923.html)  
BioConductor (http://www.bioconductor.org/)  
MiTCR  (http://mitcr.milaboratory.com/)
fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
PEAR (http://sco.h-its.org/exelixis/web/software/pear/) for paired end assembly
GEOQuery from BioConductor

## Development Notes

Still under development.

## Why? 

The basic idea behind the package is to held organize and manage the pipeline of processing `fastq` files, whether they come from data on GEO, or whether they are provided by collaborators.

The package creates a standard directory structure for each project and tracks where files are stored, how they are processed,and generates a report describing what tools are used in the processing. 

Output files are annotated and stored in standard locations and convenience functions are provided to construct `ExpressionSet`  or `SingleCellAssay` objects with full feature and phenotypic (pData) annotations. 

When working collaboratively it's important not to duplicate effort. Different people working on the same data set need to ensure that they are working with data that has been processed the same way, and if there are multiple analyses, having that data located and processed in a predictable way and stored in a standard locations becomes enables the writing of reusable code.

## Example

**The example needs to be pared down to a manageable size.**

Untar the Tab.tar.gz file under `tests` into `tests/Tab`. This will allow the
`test-all.R` script to run. 

When setting up a new project the package will download `SRAmetadb.sqlite` unless you point it to a current location.

SRA files are downloaded if not present. This is time consuming, obviously. The files are available on our computer server storage, and can be copied from there for development. 

Basic functionality exits, needs work to generalize to other tools.
