# QENIE
Correlation-based approach for identification of endocrine interactions

#Required package
The Ssec ranking score and protein-specific pathway enrichment commands in R utilize functions from the package WGCNA


#External data
In addition to tissue-specific expression arrays, this pipeline also imports information from Uniprot to filter for secreted proteins.  This list has been uploaded for ease of pipeline execution, however is subjected to updates.  Revised lists can be retrieved from the following sources:

  UniProt <www.uniprot.org/> deposited annotations for Organism:mouse "Mus musculus (Mouse) [10090]" and subcellular localization:secreted :"Secreted [SL-0243]"
  
  Citation:
  UniProt: the universal protein knowledgebase Nucleic Acids Res. 45: D158-D169 (2017)  
  
#Tissue-specific expression
Filtering for tissue-specific expression was performed by manual inspection using BioGPS on the following mouse arrays <http://biogps.org/dataset/BDS_00009/>: GeneAtlas GNF1M, gcrma
While this can be easily automatied, we felt more confident with inspection due to the small number of samples per tissue and consequent high level of variation among expression 

  GEO code for arrays: GSE1133
  
  Citation:
  McClurg P, Janes J, Wu C, Delano DL, Walker JR, Batalov S, Takahashi JS, Shimomura K, Kohsaka A, Bass J, Wiltshire T, Su AI (2007) Genomewide association analysis in diverse inbred mice: power and population structure. Genetics 176(1):675-83.

#Pipeline execution
The pipeline begins with gene expression arrays for each tissue, where each gene is represented as an averaged value across probes used in the array. The arrays used were

#We subdivided the commands into the following steps:

1. Construct cross-tissue correlation and pvalue matrices
2. Compute rowsum -log(pvalue) for ranking interactions and filter for factors proteins as secreted
3. condition correlation matrix on a by-gene basis for pathway enrichment 

Any questions/comments, please contact mseldin3@gmail.com
