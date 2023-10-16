rm(list=ls())
library(tidyverse)
library(data.table)
library(vcfR)
library(random)

setwd('~/Documents/Lab/Hybrids_RNA-seq/RNA_analyses_v3/Data')

SVs_vcf <- read.vcfR('./SV_data/out_jasmine_all.vcf')
SVs_vcf_tidy <- vcfR2tidy(SVs_vcf)$fix %>% dplyr::select(ID, CHROM, POS, SVLEN, SVTYPE, AVG_LEN)
SVs_vcf_split <- vcfR2tidy(SVs_vcf)$fix %>% 
  separate(SUPP_VEC, into = c("nope", "AAR","ABA","ABC","ABP","ABS",
  "ACG","ACK","ACS","AKE","ATE","AVI","BAK","BAM","BAN","BAP","BBS",
  "BFP","BHH","BKL","CCD","CGD","CLB","CMQ","CPG","SIG"), sep = "", remove = F) %>%
  left_join(SVs_vcf_tidy)
rm(SVs_vcf, SVs_vcf_tidy)
SVs_split2 <- SVs_vcf_split %>% select(AAR:SIG, SVTYPE) %>% 
  pivot_longer(-SVTYPE, names_to = 'isolate', values_to='presence') %>% 
  filter(presence==1)

# Overall perrcentage of each SVTYPE
ggplot(SVs_vcf_split, aes(x="", fill=SVTYPE))+
  geom_bar(stat='count', width=1, color='white')+
  coord_polar('y', start=0)+
  theme_void()
SVs_vcf_split %>% group_by(SVTYPE) %>% tally()

SVs_vcf_split %>% group_by(SVTYPE) %>% summarise(meanLength=mean(as.numeric(SVLEN)))
#Count over each parent
SVs_byParent <- SVs_split2 %>% group_by(isolate) %>% 
  summarise(SVcount=sum(as.numeric(presence)))
SVs_byParent <- SVs_byParent[order(SVs_byParent$SVcount),]

SVs_split2 %>% mutate(isolate=factor(isolate, levels=SVs_byParent$isolate)) %>% 
ggplot(aes(x=isolate, fill=SVTYPE))+
  geom_bar(stat='count')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
SVs_split2 %>% group_by(isolate) %>% tally(name='isolate_count') %>% 
  left_join(SVs_split2 %>% group_by(isolate, SVTYPE) %>% tally()) %>%
  group_by(isolate, SVTYPE) %>% summarise(SVTYPE_perc=n/isolate_count) %>% 
  mutate(isolate=factor(isolate, levels=SVs_byParent$isolate)) %>% 
  ggplot(aes(x=isolate, fill=SVTYPE, y=SVTYPE_perc))+
  geom_bar(stat='identity')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Frequency of each variant
ggplot(SVs_vcf_split, aes(x=50-abs(50-SUPP*3.846154)))+
  geom_histogram(bins=50, binwidth=4, color='black', fill='white')+
  xlim(0,50)+
  facet_wrap(~SVTYPE, scales='free_y')+
  theme_bw()
SVs_vcf_split2 <- SVs_vcf_split %>% mutate(MAF=50-abs(50-SUPP*3.846154))
ggplot(SVs_vcf_split, aes(x=50-abs(50-SUPP*3.846154)))+
  geom_density(adjust=2)+
  xlim(0,50)+
  facet_wrap(~SVTYPE)
low_freq <- SVs_vcf_split2 %>% filter(MAF<5) %>% select(ID) %>% unique() %>% nrow()
low_freq/nrow(SVs_vcf_split2)
#Variant length
SVs_vcf_split$SVLEN <- as.numeric(SVs_vcf_split$SVLEN)
ggplot(SVs_vcf_split, aes(x=abs(SVLEN)))+
  geom_density(adjust=1)+
  facet_wrap(~SVTYPE)+
  scale_x_log10()
ggplot(SVs_vcf_split, aes(x=abs(SVLEN)))+
  geom_density(adjust=1)+
#  facet_wrap(~SVTYPE, scales='free')+
  scale_x_log10()

# Pairwise genetic diversity ----
div <- fread('./SV_data/1039pairwiseDiff.tab') # You have to find this file from Victor or Anne.
parents0 <- c('AAR','ABA','ABC','ABP','ABS','ACG','ACK','ACS','AKE','ATE','AVI',
              'BAK','BAM','BAN','BAP','BBS','BFP','BHH','BKL','CCD','CGD','CLB',
              'CMQ','CPG','FY','XTRA_DJR')
div0 <- div %>% filter(strain1 %in% parents0, Strain2 %in% parents0)
write.csv(div0, './SV_data/parents_pairwise_diversity.csv')
div2 <- div %>% filter(Strain1 %in% 'FY' & Strain2 %in% parents0)
div2$Strain2[div2$Strain2=='XTRA_DJR'] <- 'SIG'
mean(div2$PercDiff)
SVs_byParent <- SVs_byParent %>% left_join(div2 %>% dplyr::rename(isolate=Strain2))
ggplot(SVs_byParent, aes(x=PercDiff, y=SVcount))+
  geom_point()+
  geom_smooth(method='lm')
lm(SVs_byParent$SVcount~SVs_byParent$PercDiff)
cor.test(SVs_byParent$SVcount,SVs_byParent$PercDiff)

# Checking Ty related SVs
tyness <- readRDS('./SV_data/tyness.Rds') %>% 
  mutate(ty=case_when(
    typercentage>0.5~'Ty',
    typercentage<=0.5~'not Ty'))
ggplot(tyness, aes(x="", fill=as.character(ty)))+
  geom_bar(stat='count', width=1, color='white')+
  coord_polar('y', start=0)+
  theme_void()
tyness %>% group_by(ty) %>% tally() %>% mutate(perc=n/nrow(tyness))
notTy <- SVs_vcf_split %>% left_join(tyness %>% select(ID, ty)) %>% filter(ty=='not Ty')
notTy %>% filter(SVTYPE!='TRA') %>% 
  ggplot(aes(x=abs(SVLEN)/1000))+
  #geom_density(adjust=1)+ 
  geom_histogram(bins=22, color='black')+
  #facet_wrap(~SVTYPE, scales='free_y')+
  scale_x_log10()
SVs_vcf_split %>% left_join(tyness %>% select(ID, ty)) %>% filter(ty=='Ty') %>% 
ggplot(aes(x=abs(SVLEN), fill=ty))+
  geom_density(adjust=0.51)+
  geom_vline(xintercept=6000)+
  #geom_histogram(bins=22, color='black', position='identity')+
    #facet_wrap(~SVTYPE, scales='free')+
  scale_x_log10()
nrow(tyness %>% filter(ty=='Ty' & SVLEN>2000))
# Checking the assembly quality thing #######
#Checking the mummer files
chr_sizes <- fread('./SV_data/R64_chrom_sizes.tsv') %>% 
  mutate('chr'=c('chrI','chrII','chrIII','chrIV','chrV','chrVI','chrVII','chrVIII',
                 'chrIX','chrX','chrXI','chrXII','chrXIII','chrXIV','chrXV','chrXVI')
         , thresh=0.9*V2)

all_files <- list.files('./SV_data/all_mummer_sdn/all_coords')


coordsfile <- paste0('./SV_data/all_mummer_sdn/all_coords/',all_files[1])
#Gauthier's function
num_contigs_cover_chr <- function(coordsfile, threshold = .9) {
  strain <- str_replace(tail(str_split(coordsfile, "/")[[1]], 1), "vsREF.coords", "")
  assembler <- tail(str_split(coordsfile, "/")[[1]], 4)[1]
  a <- data.table::fread(coordsfile, 
                    col.names = c("query_ID", "date", "len_query", "type", "ref", "ref_ID", 
                                  "start_query", "end_query", "start_ref", "end_ref", "identity_perc", 
                                  "simil_perc", "len_aln_query", "compat", "compat2", "compat3", "compat4", 
                                  "strand", "len_ref", "compat5", "compat6")) %>% 
    dplyr::select(query_ID, len_query, ref_ID, start_query, end_query, start_ref, end_ref, len_aln_query, len_ref) %>% 
    dplyr::filter(ref_ID != "ref|NC_001224|") %>% 
    dplyr::group_by(ref_ID, query_ID) %>% 
    dplyr::mutate(aligned = sum(len_aln_query),
                  chr_thresh = round(threshold*len_ref, digits = 0)) %>% 
    dplyr::distinct(query_ID, aligned, chr_thresh, len_ref) %>% 
    dplyr::arrange(ref_ID, desc(aligned)) %>% 
    dplyr::group_by(ref_ID) %>% 
    dplyr::slice(seq_len(which.max(cumsum(aligned) >= chr_thresh))) %>% 
    dplyr::group_by(ref_ID) %>% 
    dplyr::tally() %>% 
    dplyr::mutate(strain = strain,
                  assembler = assembler)
}
df <- do.call(rbind, lapply(list.files(path = "./SV_data/all_mummer_sn/all_coords", pattern = ".coords", recursive = T, full.names = T), num_contigs_cover_chr))
colnames(df)
ggplot(df, aes(x=substr(ref_ID,4,10), y=strain, fill=as.character(n)))+
  geom_tile()+
  scale_fill_viridis_d()+
  ylab('Strain')+
  xlab('Chromosome')


####################################
#SV characteristics of SV-eQTL
SV_info <- fread('./SV_data/SV_info.csv') %>% 
  select(-V1) %>% 
  dplyr::rename(SNP=SNP_ID)
eQTL <- fread('./GWAS/all_signif_not_linked.csv') %>% 
  filter(type=='SV') %>% left_join(SV_info)
write.csv(eQTL, './GWAS/for_supp_tables/SV_eQTL_for_supp.csv')
SV_eQTL <- tyness %>% filter(ID %in% eQTL$ID)
tyness %>% group_by(ty) %>% tally()


SV_eQTL %>% group_by(SVTYPE) %>% tally()
SV_eQTL %>% group_by(ty) %>% tally()

