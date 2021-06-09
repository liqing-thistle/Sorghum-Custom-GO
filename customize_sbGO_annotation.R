# Construct custom GO annotation database for sorghum
# 04/27/2021
# This file aims to construct a custom GO annotation database for sorghum using the genome annotation information from phytozome (v3.1.1)
# In the genome annotation information file (Sbicolor_454_v3.1.1.annotation_info.txt), for each sorghum gene, the GO terms and the best hit in rice and arabidopsis were annotated if there is.
# For sorghum gene with annotated GO terms, those GO terms will be used for constructing GO annotation database. 
# For sorghum gene without annotated GO terms, the GO terms of its best hit in rice will be used as alternatives. If there is no best hit in rice or best hit in rice without GO terms, the GO terms of its best hit in arabidopsis will be used.


# Load sorghum genome annotation file for total 34129 genes 
setwd()
allGenes <- read.csv(file="primarytranscript_annotation.csv")


# Extract sorghum genes with GO terms in "GO" column and save these gene ID to a new column called sbID_withGO.
# Extract best rice and arabidopsis hits of genes without GO terms, and save related best hits to a new column respectively. 
sbGoNotNull = allGenes$GO != ""
allGenes$sbID_withGO = ""   
allGenes$riceID_sbGONull = ""
allGenes$arabiID_sbGONull = ""
allGenes$sbID_withGO[which(sbGoNotNull)] = allGenes$transcriptName[which(sbGoNotNull)]
allGenes$riceID_sbGONull[which(!sbGoNotNull)] = allGenes$Best.hit.rice.name[which(!sbGoNotNull)]
allGenes$arabiID_sbGONull[which(!sbGoNotNull)] = allGenes$Best.hit.arabi.name[which(!sbGoNotNull)]


# Retrieve GO terms for genes in "sbID_withGO" column from phytozome via biomart
library("biomaRt")                                                                                                                   
mart <- useMart(biomart = "phytozome_mart", 
                dataset = "phytozome", 
                host = "https://phytozome.jgi.doe.gov")
sbGO <- getBM(attributes = c("gene_name1","transcript_name1", "go_id"), 
              filters = "transcript_name_filter", 
              values = allGenes$sbID_withGO, 
              mart = mart)


# Retrieve GO terms for genes in "riceID_sbGONull" column from phytozome via biomart 
## Because some sorghum genes have the same best hit in rice, so duplicated rice hits should be removed at first.
UniqRiceID_forGO <- unique(allGenes$riceID_sbGONull)
riceGOraw <- getBM(attributes = c("transcript_name1", "go_id"), 
                 filters = "transcript_name_filter", 
                 values = UniqRiceID_forGO, 
                 mart = mart)
## code below can be used for checking all marts of phytozome
## host="phytozome.jgi.doe.gov"
## listMarts(host=host)

# Remove blank and NA rows in riceGOraw (because some rice genes don't have GO terms)
riceGO <- riceGOraw[!(is.na(riceGOraw$go_id) | riceGOraw$go_id==""), ]


# Get rice ID having GO terms, and save it to a new column called riceID_withGO
riceIDwithGO <- as.data.frame(unique(riceGO$transcript_name1))
allGenes$riceID_withGO = ""
allGenes$riceID_withGO[which(allGenes$riceID_sbGONull%in%riceIDwithGO[, 1])] = allGenes$riceID_sbGONull[which(allGenes$riceID_sbGONull%in%riceIDwithGO[, 1])]


# Get arabidopsis hit of sorghum genes without GO terms or best rice hit not having GO terms, and save them to a new column called arabiID_forGO.
riceGoNotNull = allGenes$riceID_withGO != ""
allGenes$arabiID_forGO = ""
allGenes$arabiID_forGO[which(!riceGoNotNull)] = allGenes$arabiID_sbGONull[which(!riceGoNotNull)]


# Retrieve GO terms for genes in "arabiID_forGO" column 
## Because some sorghum genes have the same best hit in arabidopsis, so duplicated arabidopsis hits should be removed at first.
UniqArabiID_forGO <- unique(allGenes$arabiID_forGO)
arabiGOraw <- getBM(attributes = c("transcript_name1", "go_id"), 
                 filters = "transcript_name_filter", 
                 values = UniqArabiID_forGO, 
                 mart = mart)

# Remove blank and NA rows in arabiGOraw (because some arabidopsis genes don't have GO terms)
arabiGO <- arabiGOraw[!(is.na(arabiGOraw$go_id) | arabiGOraw$go_id==""), ]


# Get arabidopsis ID having GO terms, and save it to a new column called arabiID_withGO
arabiIDwithGO <- as.data.frame(unique(arabiGO$transcript_name1))
allGenes$arabiID_withGO = ""
allGenes$arabiID_withGO[which(allGenes$arabiID_forGO%in%arabiIDwithGO[, 1])] = allGenes$arabiID_forGO[which(allGenes$arabiID_forGO%in%arabiIDwithGO[, 1])]

# Export allGenes file for future use
setwd()
write.csv(allGenes, file = "primarytranscript_with_customGOannotation.csv")


# Add related sorghum ID to riceGO 
sb_rice <- allGenes[, c(1,2,20)]
colnames(sb_rice) <- c("GeneID", "transcriptName", "transcript_name1")
library(dplyr)
sb_riceGO <-  left_join(riceGO, sb_rice, by = "transcript_name1")
## the row of sb_riceGO is 2017 while the row of riceGO is 1486. The difference is caused by the fact that some sorghum genes have the same rice hit but when we retrieve rice GO term, we used unique rice ID.

# Add related sorghum ID to arabiGO
sb_arabi <- allGenes[, c(1,2,22)]
colnames(sb_arabi) <- c("GeneID", "transcriptName", "transcript_name1")
sb_arabiGO <-  left_join(arabiGO, sb_arabi, by = "transcript_name1")

  
# Combine GO terms of sorghum and best rice hit and arabidopsis hit together to get custom GO annotation database for sorghum.
colnames(sbGO) <- c("GeneID", "transcriptName", "go_id")
sbGO$homolog <- ""
sb_riceGO <- sb_riceGO[, c(3,4,2,1)]
sb_arabiGO <- sb_arabiGO[, c(3,4,2,1)]
colnames(sb_riceGO) <- c("GeneID", "transcriptName", "go_id", "homolog")
colnames(sb_arabiGO) <- c("GeneID", "transcriptName", "go_id", "homolog")

sb_rice_arabi_allGO1 <- rbind(sbGO, sb_riceGO, sb_arabiGO)
sb_rice_arabi_allGO2 <- sb_rice_arabi_allGO1[, c(1,3)]

# export file
write.csv(sb_rice_arabi_allGO1, file="sb_rice_arabi_GOannot_full.csv")
write.csv(sb_rice_arabi_allGO2, file="sb_rice_arabi_GOannot_slim.csv")

















