rm(list=ls())
library(tidyverse)
library(data.table)
setwd('~/Desktop/Diallel_RNAseq_project/Herit_stuff/Data')

# Global functions ----
detach_packages <- function(pkg_list, character.only = FALSE)
{
  for (pkg in pkg_list){
  if(!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }}
}
#####
#================================================================#
# PART 1 : Description of the experiments and overall data set ====
#================================================================#
  ## Parents in the 1011 population ----
library(ape)
library(ggtreeExtra)
library(ggtree)
library(viridis)
bam_info <- fread('./constant/bam_files.csv', header=T)
#Still missing BY and SIGMA !!!!!!!!!!!!!!
my_strains <- c("AAR", "ABA", "ABC", "ABP", "ABS", "ACG", "ACK", "ACS", "AKE", "ATE", "AVI", "BAK", "BAM", "BAN", "BAP", "BBS", "BFP", "BHH", "BKL", "CCD", "CGD", "CLB", "CMQ", "CPG")
highlight_colour <- '#E74C3C'
  
load('./constant/species_tree/tree_assets.RData')
cladesdf2 <- cladesdf %>% mutate(inSubset='In subset')
cladesdf2$inSubset[!(cladesdf$STD_name %in% my_strains)] <- 'NA'

cols <- c(viridis(length(unique(cladesdf2$clades))), 'In subset'=highlight_colour,'NA'="transparent")
names(cols) <- c(unique(cladesdf2$clades),'In subset','NA')
tr2 <- ggtree(tree1011, 
              layout = 'ape', 
              open.angle = 10, 
              size = 0.2, 
              aes(color=clades)) %<+% 
  cladesdf2 +
  geom_tippoint(mapping = aes(color = inSubset), size = 5, show.legend = T)+
  scale_color_manual(values=cols)
tr2

rm(bam_info, tr2, tree1011, cladesdf, cladesdf2, cols, highlight_colour, my_strains)
detach_packages(c('ggtreeExtra', 'ggtree',  'ape', 'viridis'), T)

  ## Parental ecological and geographical origins ----
#my_strains <- c("BY","SIG","AAR", "ABA", "ABC", "ABP", "ABS", "ACG", "ACK", "ACS", "AKE", "ATE", "AVI", "BAK", "BAM", "BAN", "BAP", "BBS", "BFP", "BHH", "BKL", "CCD", "CGD", "CLB", "CMQ", "CPG")
my_strains <- fread('./constant/26_parents_info.csv')
my_strains %>% 
  ggplot(aes(x=Class, fill=Class))+
  geom_histogram(stat='count')+
  ylab('Number of strains')+
  xlab('Strains group')

my_strains %>% group_by(Eco_Origin) %>% tally() %>% 
ggplot(aes(x="", y=n, fill=Eco_Origin)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)
# Manually reorder the groups in Illustrator if you're going to merge this plot with the one above
rm(my_strains)


  ## Parents' geographical origins ----
my_strains <- fread('./constant/26_parents_info.csv') %>% dplyr::rename(Country=Geo_Origin) %>% 
  group_by(Country) %>% tally()
country_positions <- fread('./constant/species_tree/country_positions.csv') %>% 
  filter(Country %in% my_strains$Country) %>% left_join(my_strains)
world <- map_data("world")

ggplot() +
  geom_map(
    data = world, map = world,
    aes(long, lat, map_id = region),
    fill='blue', linewidth=0
  )+
  geom_point(data=country_positions, aes(x=Longitude, y=Latitude, size=n))




  ## Correlation with 969 data set ====
    ### Data import ----
library(data.table)
library(tidyverse)
diallel_tpm <- fread('./all_tpm.csv') %>% 
  pivot_longer(-GeneID, names_to='Strain', values_to='tpm') %>% dplyr::rename('systematic_name'=GeneID) %>% 
  mutate(dataset='diallel')

perc0_tpm <- diallel_tpm %>% filter(tpm==0) %>% group_by(systematic_name) %>% tally() %>% 
  mutate(perc=n/length(unique(diallel_tpm$Strain))) %>% filter(perc>=0.5)


expressed_genes <- diallel_tpm %>% group_by(systematic_name) %>% summarise(mean_tpm=mean(tpm)) %>% 
  filter(mean_tpm>0)
diallel_tpm <- diallel_tpm %>% filter(!(systematic_name %in% perc0_tpm$systematic_name))
pop_tpm <- fread('./constant/Sace_RNA_data/final_data_annotated_merged_04052022.tab') %>% 
  select(systematic_name, Strain, tpm) %>% unique() %>% mutate(dataset='nat_pop')
    ### Abundance correlation ----
dat <- rbind(diallel_tpm, pop_tpm) %>% mutate(tpm=log(tpm+1))
grouped <- dat %>% group_by(systematic_name, dataset) %>% 
  summarise(
    mean_expr=mean(tpm),
    var_expr=var(tpm)
  ) %>% na.omit()
temp <- grouped  %>% pivot_wider(-var_expr, values_from='mean_expr', names_from='dataset') %>% na.omit()
temp %>% ggplot(aes(x=nat_pop, y=diallel))+
  geom_abline(slope=1, intercept=0, linetype=2)+
  geom_point(alpha=0.1)+
  ggtitle('Mean expression 969 isolates vs diallel')+
  theme_classic()
cor.test(temp$nat_pop, temp$diallel)
a <- lm(temp$diallel ~ temp$nat_pop)
# Expression correlates well between the 2 datasets r2=0.8, lm_slope=0.93,
# on average expression is globally lower in the natural population by 0.1tpm 
# (probably because more ORFs were used to calculate the tpm in the 1011 than in the diallel)
length(unique((diallel_tpm %>% filter(substr(systematic_name,1,1)=='Y'))$systematic_name))
    ### Dispersion correlation ----
temp <- grouped  %>% pivot_wider(-mean_expr, values_from='var_expr', names_from='dataset') %>% na.omit() 
ggplot(temp, aes(x=nat_pop, y=diallel))+
  geom_point(alpha=0.1, aes(text=systematic_name))+
  geom_abline(slope=1, intercept=0, linetype=2)+
  ggtitle('Expression variance 969 isolates VS diallel')+
  scale_y_log10()+
  scale_x_log10()+
  theme_classic()
cor.test(temp$nat_pop, temp$diallel)
lm(temp$nat_pop~ temp$diallel)


### Parents absolute expression correlation ----
cross_info <- fread('./constant/bam_files.csv', header=T) %>% dplyr::select(sampleID, P1, P2) 
diallel_tpm2 <- diallel_tpm %>% dplyr::rename(sampleID=Strain) %>% left_join( fread('./constant/bam_files.csv', header=T)) %>% 
  mutate(sample_type=case_when(
    P1==P2~'Parent',
    T~'Hybrid'
  ))
parents <- unique(c(cross_info$P1, cross_info$P2))
pop_parents <- pop_tpm %>% filter(Strain %in% parents)
diallel_parents <- diallel_tpm2 %>%
  filter(P1==P2) %>% 
  filter(P1 %in% unique(pop_parents$Strain)) %>% select(systematic_name, P1, tpm, dataset)
colnames(diallel_parents) <- colnames(pop_parents)
dat_parents <- rbind(pop_parents, diallel_parents)%>% mutate(tpm=log(tpm+1)) %>% 
  pivot_wider(names_from=dataset, values_from=tpm, values_fn=mean) %>% as.data.frame() %>% 
  na.omit()

median((dat_parents %>% mutate(ratio=nat_pop/diallel))$ratio %>% na.omit())
cor.test(dat_parents$nat_pop, dat_parents$diallel)

ggplot(dat_parents, aes(x=as.numeric(unlist(nat_pop)), y=as.numeric(unlist(diallel))))+
  #geom_point()
  geom_bin2d()+
  geom_abline(slope=1, linetype=2)+
  theme_classic()+
  scale_fill_viridis_c()


#=====================================#
# PART 2 : Estimation of H2 and h2 ====
#=====================================#
  ## H2,h2 distribution ----
rm(list=ls())
library(fgsea)
dat <- readRDS( './diallel_analyses/all_diallel_fit_new.rds')
GCA <- dat %>% as.data.frame() %>% ungroup()
GCA <- GCA %>% group_by(GeneID, P1) %>% summarise(GCA=g1) %>% unique()

herit <- dat %>% select(GeneID, h2, H2) %>% unique()
herit %>%  pivot_longer(-GeneID, names_to='herit_type', values_to='value') %>% group_by(herit_type) %>% summarise(mean_herit=mean(value))
herit %>%  pivot_longer(-GeneID, names_to='herit_type', values_to='value') %>% group_by(herit_type) %>% summarise(median_herit=median(value))
herit %>% pivot_longer(-GeneID, names_to='herit_type', values_to='value') %>% 
  ggplot(aes(x=value, color=herit_type))+
  theme(
    panel.border=element_rect(color='grey', fill='transparent'),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(color = "lightgrey", size = 0.1)
  )+
  geom_density(size=2)+
  ggtitle('Distribution of heritability')+
  facet_wrap(~herit_type, scale='free')
  
rm(GCA)
  ## Ternary plot of additive, non-additive and remaining genetic variances ----
# Run the H2 h2 plots to load the data
library(ggtern)
library(plotly)

# reusable function for creating annotation object
axis <- function(txt) {
  list(
    title = txt, tickformat = ".0%", tickfont = list(size = 10)
  )
}

ternaryAxes <- list(
  aaxis = axis("Additive"), 
  baxis = axis("Non-additive"), 
  caxis = axis("Remaining")
)
gene.names <- fread('./constant/gene_annotation/ORF2GENE.csv')
x <- fread("./constant/GO_standards/GO_Biological_Process_2018.txt", sep='\n', header=F)
x <- x$V1
GO_BP <- {}
for ( i in x){
  cur_BP <- unlist(str_split(i, '\t'))[1]
  
  genes <- unlist(str_split(i, '\t'))[2:length(unlist(str_split(i, '\t')))]
  GO_BP <- GO_BP %>% rbind(data.frame(GOterm=rep(cur_BP, length(genes)),
                                      gene_name=genes))
}

GO_BP <- GO_BP %>% filter(gene_name !='')
herit_2 <- herit %>% dplyr::rename(systematic_name=GeneID) %>%
  left_join(gene.names) %>% 
  left_join(GO_BP) %>% 
  mutate(remaining=1-H2, non_add=H2-h2) #%>% 
  #left_join(diallel_expression) 
fig <- herit_2 %>%
  select(h2, non_add, remaining, gene_name, systematic_name) %>% 
  unique() %>% na.omit() %>% 
  plot_ly(
    a = ~h2, 
    b = ~non_add, 
    c = ~remaining,
    alpha=0.2,
    type = "scatterternary"
  ) %>% 
  layout(
    ternary = ternaryAxes
  )
fig

herit_2 %>%
  select(h2, non_add, remaining, gene_name, systematic_name) %>% 
  unique() %>% na.omit() %>% 
  ggtern(aes(
    x = non_add, 
    y = h2, 
    z = remaining))+
  geom_point(color='#3498DB', alpha=0.1, size=2)+
  theme_bw()+
  ylab('Additive')+
  xlab('Non additive')+
  zlab('Remaining')
    ### Violin plot og the genetic variances ----
herit %>% mutate(remaining=1-H2, non_add=H2-h2) %>% select(GeneID, remaining, non_add, h2) %>% 
  pivot_longer(-GeneID, names_to = 'GenVar', values_to='value') %>% 
  ggplot(aes(x=GenVar, y=value))+
  #geom_violin(aes(fill=GenVar))
  geom_boxplot(aes(color=GenVar), size=1)

  ## Ternary plots for accessory genes ----
accessory <- fread('./constant/extraORFs_correspondance.csv') %>% 
  select(`Ortholog in SGD_2010`, `Annotation Name`, `Origin (closest hit)`, `Origin assignment`) %>% 
  dplyr::rename(GeneID=`Annotation Name`, Gene_type=`Origin assignment`) %>% 
  mutate(Gene_type=sub('candidate ','',Gene_type)) %>% 
  mutate(GeneID=case_when(
    (`Ortholog in SGD_2010`!=''&substr(`Ortholog in SGD_2010`,1,1)=='Y')~`Ortholog in SGD_2010`,
    T~GeneID
  ),Gene_type=case_when(
    `Ortholog in SGD_2010`!=''~paste0(Gene_type,'_ortho'),
    T~Gene_type)) %>% filter(GeneID!='YHL008C') 

herit_acc0 <- herit %>%
  mutate(non_add=H2-h2, remaining=1-H2) %>% 
  left_join(accessory)
herit_acc <- herit_acc0
herit_acc$Gene_type[is.na(herit_acc$Gene_type)] <- '0_Sace_genes'
#herit_acc$Gene_type[herit_acc$Gene_type=='candidate HGT'] <- 'HGT'
#herit_acc0$Gene_type[herit_acc0$Gene_type=='candidate HGT'] <- 'HGT'
#herit_acc$Gene_type[herit_acc$GeneID %in% accessory$`Ortholog in SGD_2010`] <- 'Introgression_ortholog'
{
  herit_acc_short <- herit_acc0 %>% select(GeneID, h2,non_add, remaining) %>% unique()
  rownames(herit_acc_short) <- herit_acc_short$GeneID
  acc_terms <- split(herit_acc$GeneID, herit_acc$Gene_type)
  acc_h2 <- herit_acc_short$h2
  acc_nonAdd <- herit_acc_short$non_add
  acc_remain <- herit_acc_short$remaining 
  names(acc_h2) <- herit_acc_short$GeneID
  names(acc_nonAdd) <- herit_acc_short$GeneID
  names(acc_remain) <- herit_acc_short$GeneID
  
  fgsea_h2 <- fgsea(pathways = acc_terms, 
                    stats    = acc_h2,
                    minSize  = 5,
                    maxSize  = 7000,eps=0)
  fgsea_h2 <- data.frame(Gene_type=fgsea_h2$pathway, h2=fgsea_h2$NES)
  fgsea_nonAdd <- fgsea(pathways = acc_terms, 
                        stats    = acc_nonAdd,
                        minSize  = 15,
                        maxSize  = 7000,
                        eps=0)
  fgsea_nonAdd <- data.frame(Gene_type=fgsea_nonAdd$pathway, nonAdd=fgsea_nonAdd$NES)
  fgsea_remain <- fgsea(pathways = acc_terms, 
                        stats    = acc_remain,
                        minSize  = 15,
                        maxSize  = 7000)
  fgsea_remain <- data.frame(Gene_type=fgsea_remain$pathway, remain=fgsea_remain$NES)
  all_NES <- fgsea_h2 %>% left_join(fgsea_nonAdd) %>% left_join(fgsea_remain)
  rm(fgsea_h2, fgsea_nonAdd, fgsea_remain, acc_h2, acc_nonAdd, acc_remain, acc_terms)
}



for(cur_type in unique(herit_acc$Gene_type)){
  print(cur_type)
  plt <- herit_acc %>% mutate(isType=case_when(
    Gene_type%in%c(cur_type)~'HGT',
    T~'0_Other_genes'
  )) %>% 
    ggtern(aes(
      x = non_add, 
      y = h2, 
      z = remaining))+
    geom_point(color='grey', alpha=0.3)+
    geom_point(data=(. %>% filter(isType=='HGT')), color='#F39C12', size=3)+
    theme_bw()+
    ylab('Additive')+
    xlab('Non additive')+
    zlab('Remaining')+
    ggtitle(cur_type)
plot(plt)
}

herit_acc_summarised <- herit_acc %>% group_by(Gene_type) %>% 
  summarise(
    h2=mean(h2),
    #median_h2=median(h2),
    nonAdd=mean(non_add),
    #median_H2=median(H2),
    remain=mean(remaining),
    #median_remaining=median(remaining),
  )
all_NES2 <- all_NES %>% pivot_longer(-Gene_type, names_to='gen_comp_fraction', values_to='NES')
acc_to_plot <- herit_acc_summarised %>% 
  pivot_longer(-Gene_type, names_to='gen_comp_fraction', values_to='value') %>% 
  left_join(all_NES2) %>% na.omit() 
acc_to_plot %>% ggplot(aes(x=gen_comp_fraction, y=Gene_type))+
  geom_point(aes(size=value, color=NES))+
  scale_color_gradient2(midpoint=mean(acc_to_plot$NES), low="#E74C3C", mid="#F5B041",
                        high="#52BE80", space ="Lab" )


GO_to_plot <- c('GO:0009086', 'GO:0006412','amino acid biosynthetic process')
for(cur_GO in GO_to_plot){
  print(cur_GO)
  plt <- herit_2 %>%
    ggtern(aes(
      x = non_add, 
      y = h2, 
      z = remaining))+
    geom_point(color='lightgrey', alpha=0.3)+
    geom_point(data=(herit_2[grepl(cur_GO,herit_2$GOterm),]), color='#52BE80', size=3, alpha=0.5)+
    theme_bw()+
    ylab('Additive')+
    xlab('Non additive')+
    zlab('Remaining')+
    ggtitle(cur_GO)
  plot(plt)
}

rm(list=ls())
  ## Parents GCA difference and ASE ----
# Data import
dat <- readRDS( './diallel_analyses/all_diallel_fit_new.rds')
GCA <- dat %>% as.data.frame() %>% ungroup()
GCA <- GCA %>% group_by(GeneID, P1) %>% summarise(GCA=g1) %>% unique()
herit <- dat %>% select(GeneID, h2, H2) %>% unique()
all_ASE_bim <- fread('./diallel_analyses/all_ASE_binomial_test.csv') %>% 
  mutate(signif=case_when(padj<0.05~1, 
                          padj>=0.05~0))
# Data processing
ASE_byGene <- all_ASE_bim %>% group_by(sampleID, gene) %>% summarise(mean_dev=mean(abs(0.5-refFreq)),
                                                                     median_padj=median(padj),
                                                                     min_padj=min(padj),
                                                                     mean_P1=mean(P1_count/totalCount),
                                                                     mean_P2=mean(P2_count/totalCount)) %>% 
  mutate(signifASE=case_when(median_padj<0.05~1, 
                             median_padj>=0.05~0))
ASE_genes <- ASE_byGene %>% filter(signifASE==1) %>% ungroup() %>% select(gene, signifASE) %>% unique() %>% dplyr::rename(GeneID=gene)
herit2 <- herit %>% left_join(ASE_genes)
# Plotting
ggplot(herit2, aes(x=H2, color=as.factor(signifASE)))+
  geom_density()+
  ggtitle('Heritability of the genes in ASE')

herit2 %>% 
  mutate(isSignifASE=case_when(
    signifASE==1~'Significant ASE',
    T~'No ASE')) %>% 
ggplot(aes(x=isSignifASE, y=H2))+
  geom_boxplot()+
  ggtitle('Heritability of the genes in ASE')

ASE_lm <- ASE_byGene %>% 
  dplyr::rename(GeneID=gene, variable=sampleID) %>% 
  left_join(dat)

ASE_lm <- ASE_lm %>% mutate(g_diff=abs(g1-g2)/mu, parent_diff=abs(P1_2n-P2_2n)/mu)
ggplot(ASE_lm, aes(x=g_diff, color=as.factor(signifASE)))+
  geom_density()+
  scale_x_log10()+
  ggtitle('difference between the GCA of the parents in ASE cases')
rm(list=ls())


#================================================#
#PART 3 : Genome-wide heritability and GWAS  ====
#================================================#
  ## Basic GWAS results----
    ### Importing data ----
not_singletons <- fread('./GWAS_results/SNP_GWAS_not_singletons.csv')
signifSNPs <- fread('../../GWAS/Data/04 GWAS_results_corrected/all_signif_not_linked_SNPs.csv')
#signifSNPs <- fread('./GWAS/01 GWAS_results/20221123_SNP_GWAS_all_signif.tsv') 
merged_hits <- fread('../../GWAS/Data/01 GWAS_results/solving_linkage_groups/SNPs_GWAS/20221123_SNP_GWAS_signifSNPs_linkage_grouped.csv') %>% 
  dplyr::rename(SNP=SNP_A) %>% select(SNP, Pheno) %>% left_join(signifSNPs) 
#signifSNPs <- signifSNPs %>% na.omit()
# Number of phenotypes per eQTL
colnames(signifSNPs)
length(unique(signifSNPs$Pheno))
length(unique(signifSNPs$SNP))
eQTL <- signifSNPs %>% group_by(SNP) %>% tally()
### GWAS supplementary figures ----
#### Number of phenotypes per eQTL
ggplot(eQTL, aes(x=n))+
  geom_histogram(stat='count', width=2)+
  ylab('Number of eQTL')+
  xlab('Number of genes influenced')+
  theme_minimal()
length(unique(signifSNPs$Pheno))
#### Number of of eQTL per phenotypes
ggplot(signifSNPs %>% group_by(Pheno) %>% tally(), aes(x=n))+
  geom_histogram(stat='count')+
  xlab('Number of eQTL')+
  ylab('Number of phenotypes')+
  theme_minimal()
median((signifSNPs %>% group_by(Pheno) %>% tally())$n)

#### Plotting the positions of eQTL across the genome
signifSNPs <- signifSNPs %>% filter(!(Pheno %in% highly_linked_phenos$Pheno))
# Plotting local and distant eQTL
plot_eQTL_locations <- function(signifSNPs){
  coords <- fread('./constant/gene_annotation/all_gene_coords.csv') %>% select(chromosome, systematic_name, start, stop)
  eQTL2 <- signifSNPs %>% left_join(coords %>% dplyr::rename(Pheno=systematic_name)) %>% 
    select(SNP, Chr, ChrPos, PValue, SnpWeight, EffectSize, Pheno, chromosome, start, stop) %>% 
    na.omit()
  chr <- fread('./constant/R64_chrom_sizes.tsv') %>% 
    mutate(Chr=as.numeric(sub('chromosome','',V1)), chrLength=V3) %>% select(Chr, chrLength)
  eQTL2 <- eQTL2 %>% left_join(chr) %>% 
    mutate(cumul_ChrPos=chrLength+ChrPos) %>% select(-chrLength) %>% 
    left_join(chr %>% dplyr::rename(chromosome=Chr)) %>%
    mutate(cumul_start=start+chrLength) %>% 
    select(-chrLength)
  ggplot(eQTL2, aes(x=cumul_ChrPos, y=cumul_start))+
    geom_abline(slope=1)+
    geom_hline(yintercept=chr$chrLength, color='red', size=0.2)+
    geom_vline(xintercept=chr$chrLength, color='red',size=0.2)+
    geom_point(size=0.2)+
    ylab('Gene position')+
    xlab('eQTL position')
}
plot_eQTL_locations(signifSNPs)

### Checking accessory genes in GWAS ----
accessory <- fread('./constant/extraORFs_correspondance.csv') %>% 
  select(`Ortholog in SGD_2010`, `Annotation Name`, `Origin (closest hit)`, `Origin assignment`) %>% 
  dplyr::rename(GeneID=`Annotation Name`, Gene_type=`Origin assignment`) %>% 
  mutate(Gene_type=sub('candidate ','',Gene_type)) %>% 
  mutate(GeneID=case_when(
    (`Ortholog in SGD_2010`!=''&substr(`Ortholog in SGD_2010`,1,1)=='Y')~`Ortholog in SGD_2010`,
    T~GeneID
  ),Gene_type=case_when(
    `Ortholog in SGD_2010`!=''~paste0(Gene_type,'_ortho'),
    T~Gene_type)) %>% filter(GeneID!='YHL008C') %>% 
  dplyr::rename(Pheno=`Ortholog in SGD_2010`)
accessory$Pheno[accessory$Pheno==''] <- accessory$GeneID[accessory$Pheno=='']
accessory2 <- accessory %>% select(Pheno) %>% mutate(isAccessory=T) 
eQTL3 <- eQTL2 %>% left_join(accessory2)
eQTL3$isAccessory[is.na(eQTL3$isAccessory)] <- F
ggplot(eQTL3, aes(x=isAccessory, y=EffectSize))+
  geom_boxplot()

ggplot(eQTL3, aes(x=isAccessory, y=abs(SnpWeight)))+
  geom_boxplot()
wilcox.test(abs(eQTL3$SnpWeight[eQTL3$isAccessory==T]),abs(eQTL3$SnpWeight[eQTL3$isAccessory==F]))
wilcox.test(eQTL3$EffectSize[eQTL3$isAccessory==T],eQTL3$EffectSize[eQTL3$isAccessory==F])

eQTL3 %>% select(Pheno, isAccessory) %>% unique() %>% group_by(isAccessory) %>% tally()

GWAS_phenos <- data.frame(Pheno= list.files('../../GWAS/Data/02_Phenotype')) %>% 
  mutate(Pheno=sub('.phen','',Pheno)) %>% left_join(accessory2) %>% 
  filter(substr(Pheno,1,1)=='Y')
GWAS_phenos %>% unique() %>% group_by(isAccessory) %>% tally()

acc_eQTL <- 245
core_eQTL <- 1436
acc <- 708
core <- 5423

eQTL <- acc_eQTL+core_eQTL
phenos <- acc+core
(acc_eQTL/acc)/(eQTL/phenos)
(core_eQTL/core)/(eQTL/phenos)
acc/phenos
acc_eQTL/acc
core_eQTL/core
fisher.test(matrix(c(acc_eQTL,eQTL,acc,phenos), nrow=2))
fisher.test(matrix(c(core_eQTL,eQTL,core,phenos), nrow=2))
acc_eQTL_enrich <- data.frame(group=c('Accessory', 'Core'), enrich=c(1.262108-1,0.9657805-1))
ggplot(acc_eQTL_enrich, aes(x=group, y=enrich))+
  geom_bar(stat='identity')+
  coord_flip()
### Plotting eQTL hotspots ----
plot_eQTL_locations <- function(signifSNPs){
  coords <- fread('./constant/gene_annotation/all_gene_coords.csv') %>% select(chromosome, systematic_name, start, stop)
  eQTL2 <- signifSNPs %>% left_join(coords %>% dplyr::rename(Pheno=systematic_name)) %>% 
    select(SNP, Chr, ChrPos, PValue, SnpWeight, EffectSize, Pheno, chromosome, start, stop) %>% 
    na.omit()
  chr <- fread('/Users/andreas/Documents/Hybrids_RNA-seq/Data/constant/Sace_ref/R64_chrom_sizes.tsv') %>% 
    mutate(Chr=as.numeric(sub('chromosome','',V1)), chrLength=V3) %>% select(Chr, chrLength)
  eQTL2 <- signifSNPs %>% 
    group_by(Chr, ChrPos, SNP) %>% tally() %>% 
    left_join(chr) %>% 
    mutate(cumul_ChrPos=chrLength+ChrPos) %>% select(-chrLength)
#    left_join(chr %>% dplyr::rename(chromosome=Chr)) %>%
#    mutate(cumul_start=start+chrLength) %>% 
#    select(-chrLength)
  ggplot(eQTL2, aes(x=cumul_ChrPos, y=n))+
    geom_vline(xintercept=chr$chrLength, color='red',size=0.3)+
    geom_segment(aes(x=cumul_ChrPos, xend=cumul_ChrPos, y=0, yend=n), size=0.3)+
    theme_minimal()
}



coords <- fread('./constant/gene_annotation/all_gene_coords.csv') %>% select(chromosome, systematic_name, start, stop)
eQTL2 <- signifSNPs %>% left_join(coords %>% dplyr::rename(Pheno=systematic_name)) %>% 
  select(SNP, Chr, ChrPos, PValue, SnpWeight, EffectSize, Pheno, chromosome, start, stop) %>% 
  na.omit()
chr <- fread('./constant/R64_chrom_sizes.tsv') %>% 
  mutate(Chr=as.numeric(sub('chromosome','',V1)), chrLength=V3) %>% select(Chr, chrLength)
eQTL2 <- eQTL3 %>% left_join(chr) %>% 
  mutate(cumul_ChrPos=chrLength+ChrPos) %>% select(-chrLength) %>% 
  left_join(chr %>% dplyr::rename(chromosome=Chr)) %>%
  mutate(cumul_start=start+chrLength) %>% 
  select(-chrLength)
ggplot(eQTL2, aes(x=cumul_ChrPos, y=cumul_start))+
  geom_abline(slope=1)+
  geom_hline(yintercept=chr$chrLength, color='red', size=0.2)+
  geom_vline(xintercept=chr$chrLength, color='red',size=0.2)+
  geom_point(size=0.2)+
  ylab('Gene position')+
  xlab('eQTL position')+
  theme_classic()


    ### Local vs Distant eQTL ----
window <- 25000
local_eQTL <- eQTL2 %>% filter(cumul_ChrPos-(cumul_start-window)>0 & cumul_ChrPos-(cumul_start+window)<0)%>% 
  mutate(eQTL_type='local')
distant_eQTL<- eQTL2 %>% filter(!(cumul_ChrPos-(cumul_start-window)>0 & cumul_ChrPos-(cumul_start+window)<0)) %>% 
  mutate(eQTL_type='distant')
nrow(distant_eQTL)/nrow(eQTL2)*100

eQTL2 <- rbind(local_eQTL, distant_eQTL)
eQTL2_table <- eQTL2 %>% select(SNP, Pheno, eQTL_type) %>% unique()
signifSNPs_table <- signifSNPs %>% left_join(eQTL2_table)

ggplot(eQTL2, aes(x=eQTL_type, y=EffectSize))+
  geom_violin()+
  geom_boxplot(width=0.3, aes(color=eQTL_type))
ggplot(eQTL2, aes(x=isAccessory, y=EffectSize, fill=eQTL_type))+
  geom_boxplot()
eQTL2 %>% group_by(eQTL_type) %>% tally()
eQTL2 %>% group_by(eQTL_type) %>% summarise
wilcox.test(local_eQTL$EffectSize, distant_eQTL$EffectSize)

#

wilcox.test(eQTL2$EffectSize[eQTL2$isAccessory==T&eQTL2$eQTL_type=='distant'],
            eQTL2$EffectSize[eQTL2$isAccessory==F&eQTL2$eQTL_type=='distant'])
wilcox.test(eQTL2$EffectSize[eQTL2$isAccessory==T&eQTL2$eQTL_type=='local'],
            eQTL2$EffectSize[eQTL2$isAccessory==F&eQTL2$eQTL_type=='local'])
wilcox.test(eQTL2$EffectSize[eQTL2$isAccessory==T],
            eQTL2$EffectSize[eQTL2$isAccessory==F])
#
eQTL2 %>% select(SNP, Pheno, isAccessory) %>% unique() %>% 
  group_by(Pheno, isAccessory) %>% tally() %>%
  ggplot(aes(x=n, color=isAccessory))+
  geom_density(adjust=2)+
  scale_y_sqrt()


rm(list=ls())
### Enrichment of rare variants in the GWAS results ----
signifSNPs <- fread('../../GWAS/Data/04 GWAS_results_corrected/all_signif_not_linked_SNPs.csv')%>% 
  dplyr::rename(SNP_name=SNP)
not_singletons <- fread('./GWAS_results/SNP_GWAS_not_singletons.csv')
signifSNPs <- fread('../../GWAS/Data/04 GWAS_results_corrected/all_signif_not_linked_SNPs.csv')
# Getting the MAFs
mafs1011 <- fread('./constant/MAFs/maf1011_20230716.csv') %>% select(Frequency, CHROM, POS) %>% 
  mutate(CHROM =as.numeric(sub('chromosome','',CHROM))) %>% 
  dplyr::rename(SNP_chr=CHROM, SNP_pos=POS, MAF1011=Frequency)
#colnames(mafs1011) <- c('SNP_chr', 'SNP_pos','MAF1011')
mafs <- read.csv('./constant/MAFs/Hybrids_SubsetGWAS_with_coords.mafs') %>% 
  select(-X) %>% left_join(mafs1011) %>% left_join(signifSNPs)%>% filter(!(is.na(MAF1011)))%>% 
  filter(SNP_name %in% not_singletons$ID)
sum(mafs$MAF_diallel<0.05)
sum(mafs$MAF1011<0.05)
sum(mafs$MAF_diallel<0.01)
sum(mafs$MAF1011<0.01)
mean(mafs$MAF_diallel)
mean(mafs$MAF_diallel)
mean(mafs$MAF1011)
median(mafs$MAF_diallel)
median(mafs$MAF1011)

1205/26957
mafs$signif <- F
mafs$signif[mafs$SNP_name%in%signifSNPs$SNP_name] <- T
mafs$rare <- F
mafs$rare[mafs$MAF1011<0.05] <- T
#nrow(mafs %>% filter(rare==T))
ggplot(mafs, aes(x=rare, y=EffectSize))+
  geom_violin(aes(color=rare))+geom_boxplot(width=0.5)
wilcox.test(mafs$SnpWeight[mafs$rare==T],mafs$SnpWeight[mafs$rare==F], alternative='greater')
ggplot(mafs, aes(x=MAF1011, y=MAF_diallel, color=EffectSize))+
  geom_point()+
  geom_vline(xintercept=0.05)+
  geom_hline(yintercept=0.05)
wilcox.test(mafs$EffectSize[mafs$rare==T],mafs$EffectSize[mafs$rare==F])
# t.test(c(1,2,3,4,5), c(2,3,4,5,6))
t.test(mafs$EffectSize[mafs$rare==T],mafs$EffectSize[mafs$rare==F])

nrow(mafs %>% filter(signif==T & rare==T))/nrow(mafs %>% filter(signif==T))
make_enrichment_plot <- function(SNPs_var){
  library(ggpubr)
  library(patchwork)
  #SNPs_var <- SNPs_var %>% filter(MAF_diallel>0.05, signif==F) %>% select(signif, rare, MAF1011, MAF_diallel)
  plot1 <- ggplot()+
    geom_point(data=filter(SNPs_var, signif==F),
               aes(x=MAF1011, y=MAF_diallel, color=signif),
               color='#5DADE2', fill='#5DADE2', stroke =0.3, shape = 21, size=2, alpha=0.1)+
    geom_point(data=filter(SNPs_var, signif==T),
               aes(x=MAF1011, y=MAF_diallel), color='#E69F00', fill='#E69F00', size=2, alpha=1)+
    geom_vline(xintercept = 0.05, linetype='dashed')+
    geom_hline(yintercept = 0.05, linetype='dashed')+
    theme_pubr()+
    ylab('Allele frequency (diallel)')+
    xlab('Allele frequeny (1011 population)')
  dens1 <- ggplot(SNPs_var, aes(x = MAF1011, fill = signif)) + 
    geom_density(alpha = 0.7, color='transparent') + 
    theme_void() + 
    theme(legend.position = "none")+
    scale_y_reverse()+
    scale_fill_manual(values=c("#5DADE2","#E69F00"))
  dens2 <- ggplot(SNPs_var, aes(x = MAF_diallel, fill = signif)) + 
    geom_density(alpha = 0.7, color='transparent') + 
    theme_void() + 
    theme(legend.position = "none")+
    scale_y_reverse()+
    coord_flip()+
    scale_fill_manual(values=c("#5DADE2","#E69F00"))
  
  dens2+plot1+plot_spacer()+dens1+
    plot_layout(ncol=2, nrow=2, widths=c(1,4), heights=c(4,1))
  
}
mafs %>% filter(!(is.na(MAF1011))) %>% group_by(MAF1011<0.05) %>% tally() 
  #group_by(signif) %>% tally()
mafs %>% select(SNP_name, rare,signif) %>% unique() %>% group_by(rare, signif) %>% tally()
make_enrichment_plot(mafs)
mafs %>% group_by(rare, signif) %>% tally() %>% 
  dplyr::rename(signif_rare_count=n) %>% left_join(mafs %>%
    group_by(rare) %>% tally() %>% dplyr::rename(rare_count=n)
  ) %>% mutate(signif_freq=signif_rare_count/rare_count) %>% 
ggplot(aes(x=rare, fill=signif, y=signif_freq))+
  geom_bar(stat='identity')
mafs %>% group_by(signif, rare) %>% tally()
fisher.test(table(mafs$signif, mafs$rare))
rm(list=ls())

rare_eQTL <- 104
rare <- rare_eQTL+1037
common_eQTL <- 927
common <- 29745+common_eQTL
not_eQTL <- 29745+1037
eQTL <- rare_eQTL+common_eQTL
all <- eQTL+not_eQTL
###Rare enrichment
rare_enrich <- (rare_eQTL/rare)/(eQTL/all)
common_enrich <- (common_eQTL/common)/(eQTL/all)
eQTL_enrich <- data.frame(type=c('Rare','Common'), enrich=c(rare_enrich-1, common_enrich-1))
ggplot(eQTL_enrich, aes(x=type, y=enrich))+
  geom_bar(stat='identity')+
  coord_flip()
rare_eQTL/eQTL
rare/all








    ### Comparing SV vs SNPs ----
# signifSNPs <- fread('./GWAS/01 GWAS_results/20221123_SNP_GWAS_all_signif.tsv')
# SNP_matrix_size <- 31972
# SV_matrix_size <- 12517
# 
# folder <- './GWAS/01 GWAS_results/raw_results/SNP_GWAS'
# read_thresholds <- function(folder){
#   all_phenos <- list.files(folder)
#   all_thresh <- {}
#   for (i in all_phenos){
#     print(i)
#     if (file.exists(paste0(folder,'/',i,'/',i,'.threshold.txt'))){
#       cur_thresh <- fread(paste0(folder,'/',i,'/',i,'.threshold.txt')) %>% mutate(Pheno=i)
#       all_thresh <- all_thresh %>% rbind(cur_thresh)
#       if (nrow(cur_thresh)>0){
#         #plt <- plotManhattan(my_dir,i)
#         #ggsave(paste0('plots/',i,'.png'), device='png')
#       }
#     } else {print(i)}
#   }
#   all_thresh <- all_thresh %>% dplyr::rename(thresh=x)
#   return(all_thresh)
# }
# 
# # Correcting and filtering the sign_eQTL 
# SNP_thresh <- read_thresholds('./GWAS/01 GWAS_results/raw_results/SNP_GWAS') %>% 
#   mutate(new_thresh=thresh*SNP_matrix_size/(SNP_matrix_size+SV_matrix_size))
# SV_thresh <- read_thresholds('./GWAS/01 GWAS_results/raw_results/20221130_GWAS_SV_wIndels') %>% 
#   mutate(new_thresh=thresh*SV_matrix_size/(SV_matrix_size+SNP_matrix_size))
# SNP_sign_corrected <- fread('./GWAS/01 GWAS_results/20221123_SNP_GWAS_all_signif.tsv') %>% 
#   left_join(SNP_thresh) %>% filter(PValue<new_thresh)
# signifSNPs <- fread('./GWAS/01 GWAS_results/20221130_GWAS_SV_wIndels_all_signif.tsv')
# SNP_groups <- fread('GWAS/03 Genotype groups/SV/0.8_SV_groups.csv', header=T) %>% select(-V1)
# size_of_groups <- SNP_groups %>% group_by(group) %>% tally()
# signifSNPs <- signifSNPs %>% filter(!(SNP %in% SNP_groups$SNP)) %>% left_join(SNP_groups)
# unlinked_SVs <- 12517-nrow(SNP_groups)
# signifSNPs2 <- signifSNPs %>% mutate(group=case_when(is.na(group)~ SNP,T~group))
# SV_sign_corrected <- all_signif %>% left_join(SV_thresh) %>% filter(PValue<new_thresh) 
# # Merging SNPs and SVs
# SNP_merge <- SNP_sign_corrected %>% select(SNP, SnpWeight, EffectSize, Pheno, thresh, new_thresh) %>% 
#   mutate(type='SNP')
# SV_merge <- SV_sign_corrected %>% select(SNP, SnpWeight, EffectSize, Pheno, thresh, new_thresh) %>% 
#   mutate(type='SV')
# all_signif <- rbind(SV_merge, SNP_merge)
all_signif <- fread('../../GWAS/Data/04 GWAS_results_corrected/all_signif_not_linked.csv')
SV_sign_corrected <- fread('../../GWAS/Data/04 GWAS_results_corrected/all_signif_not_linked.csv') %>% 
  filter(type=='SV')
all_signif %>% unique() %>% group_by(type) %>% tally() %>% 
  ggplot(aes(x=1, y=n, fill=type))+
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0)
# SNP vs SV effect size
ggplot(all_signif, aes(x=type, y=EffectSize))+
  geom_violin()+
  geom_boxplot(width=0.3)
a <- all_signif %>% filter(type=='SNP')
b <- all_signif %>% filter(type=='SV')
wilcox.test(abs(a$SnpWeight), abs(b$SnpWeight))
## Genome-wide heritability ----
    ### Data import ----
#all_g2h0 <- fread('./diallel_analyses/genome_wide_heritability/all_g2h_nonFiltered_SNPs.csv') %>% select(-V1)
all_g2h0 <- read_rds('./diallel_analyses/genome_wide_heritability/raw_results/g2h_results_filtered_log/h2data_120.rds') #%>% select(-V1)
all_g2h <- all_g2h0 %>% pivot_wider(names_from = 'Genome', values_from = 'h2', values_fn=mean)
gene.names <- fread('./constant/gene_annotation/ORF2GENE.csv')
all_g2h$Condition <- sub('./GWAS/00 phenos_forGWAS//','',all_g2h$Condition)
all_g2h$Condition <- sub('.phen','',all_g2h$Condition)
all_g2h <- all_g2h %>% dplyr::rename(systematic_name=Condition) %>% left_join(gene.names)
GOterms <- fread('./constant/gene_Annotation/go_slim_mapping.tab', header=F)
colnames(GOterms) <- c('systematic_name','description','ID','GOterm_type','GOterm','GO_ID','entry_type')
all_g2h_2 <- all_g2h %>% left_join(GOterms)
GO_group_counts <- all_g2h_2 %>% group_by(GOterm) %>% tally() %>% 
  dplyr::rename(GOgroup_size=n)


GO_groups <- all_g2h_2 %>% group_by(GOterm) %>% 
  summarise(meanSNP=mean(SNP), 
            meanSV=mean(SV),
            meanBoth=mean(both)) %>% 
  left_join(GO_group_counts) %>% 
  mutate(missing_herit=1-meanBoth)
library(ggtern)
library(plotly)
    ### SV vs SNPs heritability ----
ggtern(data=GO_groups,
       aes(x=meanSV,y=1-meanBoth, z=meanSNP)) +
  geom_point()+
  theme_bw()+
  xlab('SV heritability')+
  ylab('Remaining heritability')+
  zlab('SNP heritability')


    ### SV vs SNPs heritability ----
all_g2h_2 %>% 
  ggtern(aes(x=SV,y=1-both, z=SNP)) +
  geom_point(alpha=0.03)+
  theme_bw()+
  xlab('SV heritability')+
  ylab('Remaining heritability')+
  zlab('SNP heritability')
ggplot(all_g2h0, aes(x=Genome, y=h2))+
  geom_boxplot()
all_g2h0 %>% group_by(Genome) %>% tally()

### Comparison with the diallel h2 ====
dat <- readRDS( './diallel_analyses/all_diallel_fit_new.rds')
GCA <- dat %>% as.data.frame() %>% ungroup()
GCA <- GCA %>% group_by(GeneID, P1) %>% summarise(GCA=g1) %>% unique()

herit <- dat %>% select(GeneID, h2, H2) %>% unique() %>% 
  dplyr::rename(systematic_name=GeneID) %>% 
  left_join(all_g2h) %>% na.omit()

ggplot(herit, aes(x=h2, y=both))+
  geom_point(alpha=0.3)+
  ylab('G2H SNP+SV h2')+
  xlab('Mixed Model H2')+
  scale_color_binned()+
  geom_abline(slope=1)
cor.test(herit$h2, herit$both)

