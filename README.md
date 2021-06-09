# Sorghum-Custom-GO

This file aims to construct a custom gene ontology (GO) annotation database for sorghum bicolor using the genome annotation information from phytozome (Sorghum bicolor v3.1.1).  
  
In the genome annotation information file (Sbicolor_454_v3.1.1.annotation_info.txt), for each sorghum gene, the GO terms and the best hit in rice and arabidopsis were annotated if there is.  

For sorghum gene with annotated GO terms, those GO terms will be used for constructing GO annotation database. 
For sorghum gene without annotated GO terms, the GO terms of its best hit in rice will be used as alternatives. If there is no best hit in rice or best rice hit without GO terms, the GO terms of its best hit in arabidopsis will be used.
