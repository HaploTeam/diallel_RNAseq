rm(list=ls())
library(tidyverse)
library(data.table)

setwd('/Users/andtsouris/Documents/Lab/Hybrids_RNA-seq/RNA_analyses_v3/Data/calculating_phenotypes')

# # # # # # # # # # # # # # # # # # # # # # # # # #  #
# Merging the 3 raw_counts files into one data frame # ----
#  # # # # # # # # # # # # # # # # # # # # # # # # # #
raw_counts <- fread('raw_counts/20220523_Diallel_run_1.counts') %>% select(-V1, -Length)
#raw_counts_common1 <- raw_counts %>% select(J10, AG5, AC17, AA3, N7) %>% colSums()
raw_counts <- raw_counts %>% select(-J10, -AC17, -AA3)
colnames(raw_counts)[2:ncol(raw_counts)] <- paste0('20220523_Diallel_run_*-',
                                                   colnames(raw_counts)[2:ncol(raw_counts)])


raw_counts2 <-  fread('raw_counts/20220523_Diallel_run_2.counts')%>% select(-V1, -Length)
#raw_counts_common2 <- raw_counts2 %>% select(J10, AG5, AC17, AA3, N7) %>% colSums()
raw_counts2 <- raw_counts2 %>% select(-J10, -AG5, -AC17, -AA3, -N7)
colnames(raw_counts2)[2:ncol(raw_counts2)] <- paste0('20220523_Diallel_run_*-',
                                                     colnames(raw_counts2)[2:ncol(raw_counts2)])

raw_counts3 <- fread('raw_counts/reassigned_counts.counts')%>% select(-V1, -Length)
colnames(raw_counts3)[2:ncol(raw_counts3)] <- paste0('reassigned-',
                                                     colnames(raw_counts3)[2:ncol(raw_counts3)])


raw_counts <- raw_counts %>% left_join(raw_counts2) %>% left_join(raw_counts3)
#rm(raw_counts2, raw_counts3, raw_counts_common1, raw_counts_common2)
# Gettinf the total counts per sample
total_counts <- colSums(raw_counts %>% select(-GeneID))
total_counts <- data.frame(total_counts)
total_counts$sample <- row.names(total_counts)
# # # # # # # # # # # # # # # # # # # # # # # #
# Reassigning the sample names to cross names # ----
# # # # # # # # # # # # # # # # # # # # # # # #
crosses_done <- read.csv('../constant/crosses_done_final_forSQL.csv')

a <- str_split(crosses_done$Reads_file, '/')
a2 <- as.data.frame(t(matrix(unlist(a), nrow=length(unlist(a[1]))))) %>% 
  mutate(col_name= case_when(
    V7=="20220523_Diallel_run_*"~paste0("20220523_Diallel_run_*-",sub('.fq.gz','',V8)),
    T~paste0("reassigned-",sub('.fq','',V8))
  ))
crosses_done$col_name <- a2$col_name
giving_sampleIDs <- data.frame(col_name=colnames(raw_counts %>% select(-GeneID))) %>% 
  left_join(crosses_done %>% select(col_name, sampleID)) %>% na.omit()

GeneID <- raw_counts[,1]
raw_counts2 <- as.data.frame(t(raw_counts)) %>% mutate(col_name=rownames(.)) %>% 
  left_join(giving_sampleIDs) %>% na.omit() %>% select(-col_name)
row.names(raw_counts2) <- raw_counts2$sampleID
raw_counts3 <- as.data.frame(t(raw_counts2 %>% select(-sampleID)))
raw_counts3$GeneID <- GeneID$GeneID
Y_genes <- raw_counts3 %>% filter(substr(GeneID,1,1)=='Y' & substr(GeneID,9,11)!='BAM') %>% 
  pivot_longer(-GeneID, names_to='sample',values_to='value' ) %>% 
  mutate(value=as.numeric(value))
write.csv(Y_genes %>% select(GeneID) %>% unique(), '../constant/Y_genes.csv')

#rm(a,a2, giving_sampleIDs, GeneID, raw_counts2, raw_counts)
#!!!!!!!!!!!!!!1 Still need to figure out what to do with the 3 duplicated names from the round 2 bams !!!!!
#!!!! J10, AC17, AA3


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Merging the extraORFs to the R64 genes whenever it's possible # ----
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
orf_extra <- fread('../constant/extraORFs_correspondance.csv') %>% 
  filter(`Annotation Name` %in% raw_counts3$GeneID) 
# Giving ortholog names to the Ygenes in the extraGenes set
orf_extra <- orf_extra %>% mutate(`Ortholog in SGD_2010`=case_when(
  substr(`Annotation Name`,6,6)=='Y'~substr(`Annotation Name`,6,12),
  T~`Ortholog in SGD_2010`
))
# Filtering out the extra genes without an ortholog name
wOrtho_names <- orf_extra %>% select(`Annotation Name`, `Ortholog in SGD_2010`) %>%
  filter(`Ortholog in SGD_2010`!= '')%>% dplyr::rename(GeneID=`Annotation Name`)
weird_Ortho <- wOrtho_names %>% filter(!(`Ortholog in SGD_2010` %in% unique(Y_genes$GeneID)))
weird_Ortho <- weird_Ortho %>% mutate(`Ortholog in SGD_2010` =substr(GeneID, 6,100)) %>% 
  filter(substr(GeneID,12,100) %in% c('C-A','C-B', 'W-A', 'W-B','C-C', 'W-C'))
wOrtho_names <- wOrtho_names %>% filter(!(GeneID %in% weird_Ortho$GeneID)) %>% 
  rbind(weird_Ortho)%>%
  dplyr::rename(GeneID2=`Ortholog in SGD_2010`)
wOrtho <- raw_counts3 %>% filter(GeneID %in% wOrtho_names$GeneID)
write.csv(wOrtho_names, '../constant/genes_wOrtho_names.csv', row.names=F)
#rm(weird_Ortho)
orf_correspondance <- orf_extra %>% select(`Annotation Name`,`Ortholog in SGD_2010`) %>%
  dplyr::rename(GeneID2=`Ortholog in SGD_2010`, GeneID=`Annotation Name`) %>% 
  filter(substr(GeneID2,1,1) == 'Y') 

# Renaming the extra_genes with their ortholog names, if the ortholog names start with 'Y'
wOrtho_counts0 <- raw_counts3 %>% left_join(wOrtho_names) %>% filter(substr(GeneID2,1,1)=='Y') 
wOrtho_counts <-wOrtho_counts0 %>% mutate(GeneID=GeneID2) %>% select(-GeneID2)
  
noOrtho_genes <- fread('../constant/extraORFs_correspondance.csv') %>% 
  filter(`Annotation Name` %in% raw_counts3$GeneID) %>% 
  filter(!(`Annotation Name` %in% orf_correspondance$GeneID)) %>% 
  filter(substr(`Annotation Name`,1,1)!='Y')



# Merging the genes with Ortho names to the R64 genes
wOrtho_names_inv <- wOrtho_names
colnames(wOrtho_names_inv) <- c('GeneID2', 'GeneID')
wOrtho_counts2 <- wOrtho_counts %>% pivot_longer(-GeneID, names_to='sample',values_to='value' ) %>% 
  mutate(extra_count=as.numeric(value)) %>% select(-value) %>% 
  group_by(GeneID, sample) %>% summarise(extra_count=sum(extra_count)) %>% ungroup() %>% 
  left_join(wOrtho_names_inv)

both_genes <- Y_genes %>% left_join(wOrtho_counts2)
both_genes$extra_count[is.na(both_genes$extra_count)] <- 0
both_genes$total_count <- both_genes$value + both_genes$extra_count
#rm(Y_genes, wOrtho_names_inv, wOrtho_counts2, wOrtho_counts, wOrtho_counts0)

# Adding the extra genes that don't have ortho names
noOrtho_counts <-  raw_counts3 %>% filter(GeneID %in% noOrtho_genes$`Annotation Name`) %>% 
  pivot_longer(-GeneID, names_to='sample',values_to='value' ) %>%
  mutate(value=as.numeric(value)) %>% 
  mutate(extra_count=0, GeneID2=NA, total_count=value)
both_genes_final <- both_genes %>% rbind(noOrtho_counts)
colnames(both_genes_final) <- c('GeneID', 'sample','count','extra_count','extra_ID','total_count')
both_genes_final <- both_genes_final %>% relocate(sample, GeneID, count, extra_ID, extra_count, total_count)
#both_genes_final$sample <- sub('.','',both_genes_final, fixed=T)
write.csv(both_genes_final, './merged_counts.csv', row.names=F)
#rm(noOrtho_counts, orf_extra, wOrtho, wOrtho_counts, wOrtho_names)

# Overall we have 6684 Y_genes, 417 genes with orthologous gene names and
# 233 without orthologous gene names.

# Initial number of genes
length(raw_counts$GeneID)
#Number of extra genes merged with their native counterpart
length(unique(both_genes_final$extra_ID))
#Number of Y_genes
length(unique(both_genes$GeneID))
#Number of extra genes that are not merged
length(unique(noOrtho_counts$GeneID))

#Final number of genes
length(unique(both_genes_final$GeneID))

raw_corresp <- fread('../constant/extraORFs_correspondance.csv') %>% mutate(GeneID=`Annotation Name`)
a <- filter(raw_counts, !(GeneID %in% c(both_genes_final$GeneID, both_genes_final$extra_ID))) %>% 
  select(GeneID) %>% unique() %>% left_join(raw_corresp)

# So, I'm removing 11 extraORFs that have multiple R64 counterparts, they are paradoxus introgressions 