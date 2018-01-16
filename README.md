# QENIE (Quantitative Endocrine Network Interaction Estimation)
Correlation-based approach for identification of endocrine interactions

#Required package
The Ssec ranking score and protein-specific pathway enrichment commands in R utilize functions from the package WGCNA

#External data
In addition to tissue-specific expression arrays, this pipeline also imports information from Uniprot to filter for secreted proteins.  This list has been uploaded for ease of pipeline execution, however is subjected to updates.  Revised lists can be retrieved from the following sources:

  UniProt <www.uniprot.org/> deposited annotations for Organism:mouse "Mus musculus (Mouse) [10090]" and subcellular localization:secreted :"Secreted [SL-0243]"
  
  Citation:
  UniProt: the universal protein knowledgebase Nucleic Acids Res. 45: D158-D169 (2017)  
#These are provided as: Secreted_proteins_Uniprot
#Tissue-specific expression
Filtering for tissue-specific expression was performed by manual inspection using BioGPS on the following mouse arrays <http://biogps.org/dataset/BDS_00009/>: GeneAtlas GNF1M, gcrma
While this can be easily automatied, we felt more confident with inspection due to the small number of samples per tissue and consequent high level of variation among expression 

  GEO code for arrays: GSE1133
  
  Citation:
  McClurg P, Janes J, Wu C, Delano DL, Walker JR, Batalov S, Takahashi JS, Shimomura K, Kohsaka A, Bass J, Wiltshire T, Su AI (2007) Genomewide association analysis in diverse inbred mice: power and population structure. Genetics 176(1):675-83.

#Data pretreatment:  
Mouse expression arrays were performed on a Affymetrix HT_MG-430A 
GEO Accession for arrays: GSE64770
The pipeline begins with gene expression arrays for liver and adipose tissue, where each gene is represented as an averaged value across probes used in the array.  These aggregate matrices are also provided in this repository.  The arrays consisted of ~22,400 probes  which were aggregated to 12,242 genes to be used for cross-tissue analysis.  Therefore, each liver and adipose tissue expression matrix consists of 12,242 genes among 106 unique HMDP strains.   

#We subdivided the commands into the following steps (all listed in the R_script file):

1. Construct cross-tissue correlation and pvalue matrices
2. Compute rowsum -log(pvalue) for ranking interactions and filter for factors proteins as secreted
3. Condition correlation matrix on a by-gene basis for pathway enrichment 

Any questions/comments, please contact mseldin3@gmail.com
