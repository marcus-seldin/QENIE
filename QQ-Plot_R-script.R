#This script produces a QQ-plot estimating the relative Ssec scores for secreted vs non-secreted factors across tissues
#This script was graciously provided and devised by Simon Koplev

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

#read in Secreted Proteins
Secreted_proteins <- read.delim("Secreted_proteins_Uniprot.txt", header = T)


ranks = liv.adip.table

#create two vectors, one containing secreted factor scores and the other with non-secreted factor scores 
ranks$m = match(row.names(liv.adip.table),Secreted_proteins$Gene.names...primary.., nomatch = 0)
ranks$mm = ranks$m >0
ranks = ranks[!grepl("FALSE",ranks$mm),]
ranks$m = NULL
ranks$mm = NULL
sec = ranks

ranks = liv.adip.table
ranks$m = match(row.names(ranks),secreted.proteins[ ,1], nomatch = 0)
ranks$mm = ranks$m >0
ranks = ranks[!grepl("TRUE",ranks$mm),]
ranks$m = NULL
ranks$mm = NULL
non_sec=ranks

qqplotAnnot = function(x, y,
                       probs=c(0.001, 0.01, 0.1, 0.9, 0.99, 0.999),
                       ...)
{
  qqplot(
    x,
    y,
    pch=16, cex=0.5,
    ...
  )
  abline(0, 1, col="red")
  
  q_y = quantile(y, probs=probs, na.rm=TRUE)
  q_x = quantile(x, probs=probs, na.rm=TRUE)
  
  points(q_x, q_y, col="grey", cex=0.6)
  text(q_x, q_y, label=paste0(probs*100, "%"), pos=1, cex=0.8, col="grey")
}

qqplotAnnot(sec$liv.musc.table,non_sec$liv.musc.table,
            main="QQ-plot,Intestine x liver",
            xlab="Non-secreted factor Ssec (-ln p)",
            ylab="Secreted factor Ssec (-ln p)"
)
