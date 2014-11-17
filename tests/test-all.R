#create a project
require(RSQLite)
require(SRAdb)
require(data.table)
require(GEOquery)
require(RNASeqPipelineR)
createProject("myproject",path="./tests/",load_from_immport=TRUE)

#copy the immport tables
system("cp -r tests/Tab/* tests/myproject/Tab/")

#Load the immport data
loadImmportTables()

#Download the SRAdb database (if necessary)
getSRAdb(path="tests/Utils")

#truncate for testing
devel_truncateData(n=4)

#expect 4 rows
nrow(getConfig()[["immport_tables"]][["GSM_table"]])==4

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

#build the reference genome
buildReference(gtf_file="UCSC.gtf",fasta_file="hg38.fa",name="hg38")

#Align and compute expression counts
RSEMCalculateExpression(ncores=4)

#Assemble an expression matrix of counts and tpm and save them to output files.
RSEMAssembleExpressionMatrix()

#Annotate using bioconductor
BioCAnnotate(annotation_library="TxDb.Hsapiens.UCSC.hg38.knownGene",force=FALSE)

#output version info
pipelineReport()

# For existing project load and grab data
loadProject(project_dir = "tests",name="myproject")
eset<-getExpressionSet()
