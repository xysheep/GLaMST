# Load required packages
library(alakazam)
library(igraph)
library(dplyr)

# Select clone from example database
data(ExampleDb)
sub_db <- subset(ExampleDb, CLONE == 3138)

clone <- makeChangeoClone(sub_db, text_fields=c("SAMPLE", "ISOTYPE"),num_fields="DUPCOUNT")

# Show combined annotations
clone@data[, c("SAMPLE", "ISOTYPE", "DUPCOUNT")]

dnapars_exec <- "lib/dnapars.exe"
graph <- buildPhylipLineage(clone, dnapars_exec, rm_temp=TRUE)
