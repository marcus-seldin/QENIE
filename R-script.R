#Note: This example performs QENIE on the liver-to-adipose circuit.  For purposes of replication and ease, other datasets are provided whereby similar analyses could easily be applied

#Import your data
Adipose <- read.delim('adipose.txt', check.names=F)
Liver <- read.delim('liver.txt', check.names=F)

#Using these two datasets we will proceed:

#Construct cross-tissue correlation and pvalue matrices
library(WGCNA)
liv.adip.cor = bicor(Liver, Adipose, use='pairwise.complete.obs')
liv.adip.p = corPvalueStudent(liv.adip.cor, nSamples = 106)

#Compute rowsum -log(pvalue) for ranking interactions and filter for factors proteins as secreted
ll = as.data.frame(liv.adip.p)
liv.adip.log = -log(ll)

liv.adip.log = as.data.frame(liv.adip.log)

#compute the rowsum
liv.adip.table = rowSums(liv.adip.log, na.rm = T)
liv.adip.table = as.data.frame(liv.adip.table)

#import secreted peptides
Secreted_proteins <- read.delim(QENIE/Secreted_proteins_Uniprot, header = T)


#retain only secreted peptides
liv.adip.table$m = match(row.names(liv.adip.table),Secreted_proteins$Gene.names...primary.., nomatch = 0)
liv.adip.table$mm = liv.adip.table$m >0
liv.adip.table = liv.adip.table[!grepl("FALSE",liv.adip.table$mm),]
liv.adip.table$m = NULL
liv.adip.table$mm = NULL

#Normalize the Ssec score by number of target tissue probes (In this case, adipose contains 12242 genes)
liv.adip.table$Ssec = liv.adip.table[,1]/length(liv.adip.log)
liv.adip.table$liv.adip.table = NULL

#order by "sig score"
final = liv.adip.table[order(liv.adip.table$Ssec, decreasing=T), , drop = FALSE]
write.table(final, file="Liver X Adipose ranked by sig score",row.names=T, col.names=T, sep='\t', quote=F)

#This produces a table to each secreted protein and its respective significance score across adipose transcripts
#Note that Notum is listed as the 5th with an Sssec of 4.148
#For our pipeline, we next check the tissue-specificty using BioGPS - note that this step is not necessary, but makes us more confident when conditioning the pathway enrichment.  This moves Notum up to the 2nd ranked protein  

#Condition correlation matrix on a by-gene basis for pathway enrichment - this example will focus on the protein, Notum
#remove gene of interest (Notum) from the correlation matrix
#the rownames of the file correspond to gene symbols for pathway enrichment, whereas the second column contains the bicor coefficent


#positively correlated pathways - "enhanced by protein"
bicor.data = as.data.frame(liv.adip.cor)

Notum = bicor.data["Notum",]
Notum  = t(Notum )
Notum  = as.data.frame(Notum)
target = Notum[order(Notum, decreasing=T), , drop = FALSE]
target = head(target, 500)
write.table(target, file="Positive Notum Liver X Adipose Pathways Enrichment File", col.names=F, sep='\t', quote=F)

#negatively correlated pathways - "suppressed by protein"
Notum = bicor.data["Notum",]
Notum  = t(Notum)
Notum  = as.data.frame(Notum)
target = Notum[order(Notum, decreasing=F), , drop = FALSE]
target = head(target, 500)
write.table(target, file="Negative Notum Liver X Adipose Pathways Enrichment File", col.names=F, sep='\t', quote=F)

#All pathways engaged by protein
Notum = bicor.data["Notum",]
Notum  = t(Notum)
Notum  = as.data.frame(abs(Notum))
target = Notum[order(Notum, decreasing=T), , drop = FALSE]
target = head(target, 500)
write.table(target, file="Absolute Notum Liver X Adipose Pathways Enrichment File", col.names=F, sep='\t', quote=F)
