rm(list=ls())
library(tidyverse)
library(data.table)
setwd('/Users/andtsouris/Documents/Lab/Hybrids_RNA-seq/RNA_analyses_v3/Data/GWAS')


#Getting the tpm values and crosses info
tpm <- read.csv('../calculating_phenotypes/all_tpm.csv')
crosses <- fread('../constant/crosses_done_final_forSQL.csv') %>% 
  dplyr::rename(cross=Cross, sample=sampleID) %>% select(-Reads_file)

# Making long tpm dataframe
tpm2 <- tpm %>% pivot_longer(-GeneID, names_to='sample', values_to='tpm') %>% 
  mutate(sample=sub('X','',sample)) %>% left_join(crosses)

# Finding which genes have tpm=0 in more than half of the hybrids
perc0_tpm <- tpm2 %>% filter(tpm==0) %>% group_by(GeneID) %>% tally() %>% 
  mutate(perc=n/length(unique(tpm2$sample)))
too_many_0 <- perc0_tpm %>% filter(perc>0.5) %>% select(GeneID)
#In total, 731 genes have tpm=0 in more than half of the hybrids and are excluded from GWAS.

# Generating one phenotype file for each gene, in total 6186 files
for(cur_gene in unique(tpm2$GeneID)){
  if(cur_gene %in% too_many_0$GeneID){
  } else{
  cur_pheno <- tpm2 %>% filter(GeneID==cur_gene) %>% select(cross, tpm)
  write.table(cur_pheno,
              paste0('./phenos_forGWAS/',cur_gene,'.phen'),
              sep='\t', row.names=F, quote=F, col.names=F)  
  }
}
