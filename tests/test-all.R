#create a project
require(SQLite)
require(SRAdb)
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

conf <- getConfig()
