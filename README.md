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
bowtie2

## Development Notes

Still under development.

**The example needs to be pared down to a manageable size.**

Untar the Tab.tar.gz file under `tests` into `tests/Tab`. This will allow the
`test-all.R` script to run. 

When setting up a new project the package will download `SRAmetadb.sqlite` unless you point it to a current location.

SRA files are downloaded if not present. This is time consuming, obviously. The files are available on our computer server storage, and can be copied from there for development. 

Much of the proposed API still needs to be implemented. The project is definitely not yet ready for use.

