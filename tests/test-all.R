#create a project
require(RSQLite)
require(SRAdb)
require(data.table)
require(GEOquery)
createProject("myproject",path="./tests/",load_from_immport=TRUE)

#Should fail since we have no immport tables
system("rm -rf ./tests/myproject/Tab/*")
loadImmportTables()

#copy the immport tables
system("cp -r tests/Tab/* tests/myproject/Tab/")

#Load the immport data
loadImmportTables()

#Download the SRAdb database (if necessary)
getSRAdb(path="tests/Utils")

#Detect aspera location
detectAspera()

#truncate for testing
devel_truncateData(n=4)
getConfig()[["immport_tables"]][["GSM_table"]]

#download SRA files 
downloadSRA()

#annotate the pData with SRA provenance and write out to RSEM
annotateSRA()

#Convert SRA files to FASTQ
convertSRAtoFastQ(ncores=4)

#run fastQC
runFastQC(ncores=4)

#summarize FastQC results
summarizeFastQC()

buildReference(gtf_file="UCSC.gtf",fasta_file="hg38.fa",name="hg38")
