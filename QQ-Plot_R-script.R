#This script produces a QQ-plot estimating the relative Ssec scores for secreted vs non-secreted factors across tissues
#This script was graciously provided and devised by Simon Koplev.  In addition, We would like to thank Arjun Sarathi and Manikandan (Mani) Narayanan at the Bioinformatics and Integrative Data Sciences (BIRDS) Lab at IIT Madras for their thorough troubleshooting and identifying a switch in the initially-described QQ-Plot Axes.

library(WGCNA)
#read in Secreted Proteins
Secreted_proteins <- read.delim("Secreted_proteins_Uniprot.txt", header = T, check.names = F)

#Import your data
Adipose <- read.delim('adipose.txt', check.names=F)
Liver <- read.delim('liver.txt', check.names=F)

#Using these two datasets we will proceed:

# Calculate cross-tissue correlation p-values
liv.adip.p = bicorAndPvalue(Liver, Adipose, use='pairwise.complete.obs')$p
scores = rowSums(-log(liv.adip.p))

# Create two vectors, one containing secreted factor scores and the other with non-secreted factor scores
sec = scores[names(scores) %in% Secreted_proteins$`Gene names  (primary )`]
non_sec = scores[!(names(scores) %in% Secreted_proteins$`Gene names  (primary )`)]

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

qqplotAnnot(non_sec, sec,
            main = "QQ-plot, Liver x adipose",
            xlab = "Non-secreted factor Ssec (-ln p)",
            ylab = "Secreted factor Ssec (-ln p)"
)
