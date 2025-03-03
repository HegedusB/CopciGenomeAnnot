# library("movAPA")
library("ComplexHeatmap")
library("tidyr")
library("stringr")
library("readr")
library("tibble")
library('docstring')

work_description <- 'APA workflow on without movAPA'
options(width = 300)

DT <- format(Sys.time(), '%Y%m%d')

# set wd ---------------------------------------------------

# load PAC count data ----------------------------- 
pac_raw_df <- read_tsv('../pac_with_max_pas_readsupport_annot_pasreadsupport_noNA_annotEntry_FilteredRelevantPAC.tsv')

pacMeta <- pac_raw_df[1:27]
pacExp <- pac_raw_df[28:ncol(pac_raw_df)]

# annotation check / summary ----------------------------------------------------------------------------------------------
sum(is.na(pac_raw_df$transcript_id_nn))

length(unique(pac_raw_df$transcript_id_nn))
min(pac_raw_df$read_support)
pac_raw_df %>% group_by(transcript_id_nn) %>% summarise(pac_n = n()) %>% ungroup() -> pac_n_summary
summary(pac_n_summary$pac_n)
boxplot(pac_n_summary$pac_n)
sum(pac_n_summary$pac_n == 1)
barplot(table(cut(pac_n_summary$pac_n, breaks=c(1:20,100), include.lowest=TRUE, right=FALSE)), las=2)
table(cut(pac_n_summary$pac_n, breaks=c(1:20,100), include.lowest=TRUE, right=FALSE)) %>% as.data.frame() %>% knitr::kable()

pac_raw_df %>% group_by(transcript_id_nn) %>% summarise(overlap_index = sum(overlap)) -> pac_annot_3utr_overlap_index
table(pac_annot_3utr_overlap_index$overlap_index)

pac_raw_df %>% 
  group_by(transcript_id_nn) %>% 
  summarise(maxRsupportDist = distance_nn_correct[pas_max_read_support == max(pas_max_read_support)][[1]]) %>% 
  ungroup() -> maxRsupportDist_tbl

table(maxRsupportDist_tbl$maxRsupportDist == 0)
max(maxRsupportDist_tbl$maxRsupportDist)

table(cut(maxRsupportDist_tbl$maxRsupportDist, c(0,1,10,100,200,400,800,1000), include.lowest=TRUE, right=FALSE)) %>% knitr::kable()
barplot(table(cut(maxRsupportDist_tbl$maxRsupportDist, c(0,1,10,20,40,60,80,100,200,400,800,1000),  include.lowest=TRUE, right=FALSE)))


# Test dataset ------------------------------------------------------
only_count_sum <- sort(colSums(pacExp))

# boxplot
boxplot(only_count_sum, main="PAC read support / sample")
abline(h=5E5, lt=2)

summary(only_count_sum)

# ecdf
plot(ecdf(only_count_sum), main="ecdf")

# barplot
barplot(sort(only_count_sum, decreasing = TRUE),las=2)


quantile(only_count_sum,0.9)
only_count_sum[only_count_sum<=quantile(only_count_sum,0.1)] %>% knitr::kable()
only_count_sum[only_count_sum>=quantile(only_count_sum,0.9)] %>% sort(., decreasing = TRUE)%>% knitr::kable()

only_count_sum[str_detect(names(only_count_sum), pattern="CTGT")]

library(factoextra)
library(FactoMineR)
library(ggrepel)
library(gplots)

PA_PCA <- PCA(log2(pacExp+0.01), graph = F)

as.data.frame(PA_PCA$var$coord) %>% 
  rownames_to_column(var="exp") %>% 
  ggplot(aes(x=`Dim.1`, y=`Dim.2`, label=exp)) + 
  geom_point() + 
  geom_text_repel(size=3)

ggplot2::ggsave(str_c("PCA_PA_count_", DT, ".pdf"), device = "pdf", width = 40, height = 40, limitsize = FALSE)

group_info_tbl <- read_tsv("../copci_meta_all_pecs.tsv")


sample_name <- names(pacExp)

setdiff(sample_name, group_info_tbl$sample)
setdiff(group_info_tbl$sample, sample_name)

group_name <- group_info_tbl$treatment[match(sample_name, group_info_tbl$sample)]

colData <- data.frame(sample = sample_name, 
                      group = group_name)

colData <- mutate(colData, sample = str_replace_all(sample, pattern = "-", replacement="_"),
                  group = str_replace_all(group, pattern = "-", replacement="_"))


colnames(pac_raw_df)[28:ncol(pac_raw_df)] <- str_replace_all(colnames(pac_raw_df[28:ncol(pac_raw_df)]), pattern = "-", replacement = "_")


setdiff(colnames(pac_raw_df)[28:ncol(pac_raw_df)],colData$sample)
setdiff(colData$sample,colnames(pac_raw_df)[28:ncol(pac_raw_df)])

# create PAC database -----------------------------------

# Annotate PACs ----------------------------------------------------------------------------------------------------------


pac_raw_df <- mutate(pac_raw_df, pac_index = str_c('PAC_', pac_index))

pacMeta <- pac_raw_df[1:27]
pacExp <- pac_raw_df[28:ncol(pac_raw_df)]
rownames(pacExp) <- pacMeta$pac_index

# CPM normalisation
getSizeFactor<-function(dat, method='DESEQ', oLibsize=FALSE, thd=0.95) {
  
  method=toupper(method)
  if (is.data.frame(dat)) dat=as.matrix(dat)
  if (method!='EDGER' && method!='DESEQ' && method!='TPM' && method!='TPMX') {
    method='DESEQ'
  }
  stopifnot (method %in% c('EDGER','DESEQ','DESEQ2','CPM','TPM', 'TPMX'))
  if (method=='EDGER') {
    suppressPackageStartupMessages( library( "edgeR" ) )
    f <- calcNormFactors(dat)
  } else if (method %in% c('DESEQ','DESEQ2')) {
    suppressPackageStartupMessages( library( "DESeq2" ) )
    ds=DESeqDataSetFromMatrix(countData = dat,
                              colData = DataFrame(1:ncol(dat)),
                              design = ~ 1)
    f=sizeFactors(estimateSizeFactors(ds))
  } else if (method=='TPM') {
    f=colSums(dat)/1000000
  }else if (method=='TPMX') {
    rowsums=rowSums(dat);
    uprow=quantile(rowsums,probs=thd)
    rest=dat[rowsums<uprow,]
    libs=colSums(rest);
    minlib=min(libs)
    f=libs/minlib
  }
  if (oLibsize){
    f=as.integer(colSums(dat)/f)
  }
  return(f)
}

normBySizeFactor<-function(dat, sizeFactor) {
  tmp=matrix(rep(sizeFactor,nrow(dat)), nrow=nrow(dat), byrow=T)
  return(round(dat/tmp))
}

pacExp_norm <- pacExp
sf <- getSizeFactor(pacExp_norm, method = "TPM")
pacExp_norm <- normBySizeFactor(pacExp_norm, sf)

barplot(sort(sf, decreasing=TRUE), las=2)
sf[sf<=quantile(sf,0.1)] %>% sort(.) %>% knitr::kable() 
sf[sf>=quantile(sf,0.9)] %>% sort(.,decreasing = TRUE) %>% knitr::kable() 

barplot(colSums(pacExp_norm), las=2, ylim=c(0, 1E6))
abline(h=1e6)

# PCA norm count 
PAnorm_PCA <- PCA(log2(pacExp_norm+0.01), graph = F)

as.data.frame(PAnorm_PCA$var$coord) %>% 
  rownames_to_column(var="exp") %>% 
  ggplot(aes(x=`Dim.1`, y=`Dim.2`, label=exp)) + 
  geom_point() + 
  geom_text_repel(size=3)

ggplot2::ggsave(str_c("PCA_PA_norm_count_", DT, ".pdf"), device = "pdf", width = 40, height = 40, limitsize = FALSE)


# pool biological replicates ----------------------------------------------------------------------------------------------------------

colData_lt <- split(colData, colData$group)

pacExp_normPool_lt <- lapply(colData_lt, function(x) {
  result <- as_tibble(rowSums(pacExp_norm[x$sample]), rownames='pac_index')
  return(column_to_rownames(result, var = 'pac_index'))
})

for(i in names(pacExp_normPool_lt)) {
  names(pacExp_normPool_lt[[i]])[[1]] <- i
}

pacExp_normPool_tbl <- bind_cols(pacExp_normPool_lt)


barplot(colSums(pacExp_normPool_tbl), las=2, ylim=c(0, 4E6))
abline(h=1e6)

# PCA norm pool count 
PAnormPool_PCA <- PCA(log2(pacExp_normPool_tbl+0.01), graph = F)

as.data.frame(PAnormPool_PCA$var$coord) %>% 
  rownames_to_column(var="exp") %>% 
  ggplot(aes(x=`Dim.1`, y=`Dim.2`, label=exp)) + 
  geom_point() + 
  geom_text_repel(size=3)

ggplot2::ggsave(str_c("PCA_PA_normPool_count_", DT, ".pdf"), device = "pdf", width = 40, height = 40, limitsize = FALSE)


as.data.frame(PAnormPool_PCA$var$coord[,1:2]) %>% 
  filter(rownames(.)!="CTGT") %>%
  t(.) %>% 
  heatmap.2(scale = "row", trace = "none", col="bluered", cexCol =  1, cexRow= 1.5,margins = c(8, 8)) # a PCA 1-2 DIM lett klaszterezve

# filter genes with 3UTR APA ----------------------------------------------------------------------------------------------------------

pacMeta %>% 
  group_by(transcript_id_nn) %>% 
  mutate(utr3_length = ifelse(cds_strand == '+', pas_start-cds_pos,
                              cds_pos-pas_start)) %>% 
  ungroup() -> pacMeta

max(pacMeta$utr3_length)

pac_normPool_df <- left_join(pacMeta, rownames_to_column(pacExp_normPool_tbl, var='pac_index'), by='pac_index')

pac_normPool_df %>% 
  group_by(transcript_id_nn) %>% 
  mutate(APA_index = n()>1) %>% 
  ungroup() -> pac_normPoolAPAindex_df

dplyr::select(pac_normPoolAPAindex_df, pac_index, transcript_id_nn, utr3_length, APA_index, 29:ncol(pac_normPool_df)) %>% 
  dplyr::rename('PA'='pac_index', 'gene'='transcript_id_nn', 'three_UTR_length'='utr3_length') %>% 
  gather(key="exp", value="count", -1, -2, -3, -4) -> pac_normPoolShortAPAindex_df


write_tsv(pac_normPoolAPAindex_df, 'pac_with_max_pas_readsupport_annot_pasreadsupport_noNA_FilteredRelevantPAC_NormAPAindex.tsv')
write_tsv(pac_normPoolShortAPAindex_df, 'pac_with_max_pas_readsupport_annot_pasreadsupport_noNA_FilteredRelevantPAC_NormShortAPAindex.tsv')

pac_normPoolAPAindex_df %>% 
  filter(APA_index) %>% 
  dplyr::count(transcript_id_nn) %>% 
  summarise(min_n = min(n),
            max_n = max(n),
            tscript_n = n())


quantseq_pooled_cpm <- read_tsv("/node9_data/bhegedus/pjx_Copci_Quantseq_20220323/DEG_analysis_20220913/All_comparisons_20221021/CopciAB_Quantseq_cpm_pooled.tsv")
quantseq_pooled_cpm <- mutate(quantseq_pooled_cpm, gene_id = str_remove(gene_id, pattern = '\\.T0$')) 

setdiff(names(pac_normPool_df[29:ncol(pac_normPool_df)]), names(quantseq_pooled_cpm)[-1])
quantseq_pooled_cpm <- quantseq_pooled_cpm %>% dplyr::rename('OIDIA_1' = 'OIDIA_0h', 'OIDIA_2' = 'OIDIA_18h', 'LoKo_1' = 'LoKo_Tier2', 'LoKo_2' = 'LoKo_Tier3')
setdiff(names(pac_normPool_df[29:ncol(pac_normPool_df)]), names(quantseq_pooled_cpm)[-1])

pacQuantseqCor <- sapply(names(pac_normPool_df[29:ncol(pac_normPool_df)]), function(x) {
  pac_normPool_df[c('transcript_id_nn', x)] %>% 
    setNames(c('gene_id', 'count')) %>% 
    group_by(gene_id) %>% 
    summarise(countSum = sum(count)) %>% 
    ungroup() %>% 
    left_join(quantseq_pooled_cpm[c('gene_id', x)], by='gene_id') %>% 
    dplyr::select(-1) -> result
  
  return(cor(result[!is.na(result[[2]]),])[1,2])
})

as_tibble(pacQuantseqCor, rownames='sample') %>% arrange(desc(value)) -> pacQuantseqCor_tbl

knitr::kable(pacQuantseqCor_tbl)


# Diversity: Shannon index, Evenness ----------------------------------------------------------------------------------------------------------

pac_normPoolShortAPAindex_df

pac_normPoolShortAPAindex_df %>%
  filter(APA_index) %>% 
  group_by(exp, gene) %>%
  filter(sum(count)>=10) %>%
  ungroup() -> gene_utr_length_count_short_50_tbl

gene_utr_length_count_short_50_tbl %>% group_by(exp) %>% summarise(APA_n = length(unique(gene))) %>% ungroup() -> apa_n_50_tbl
gene_utr_length_count_short_50_tbl %>% group_by(exp,gene) %>% summarise(n=n()) %>% arrange(n)
gene_utr_length_count_short_50_tbl %>% group_by(exp,gene) %>% summarise(n=n()) %>% arrange(desc(n))
gene_utr_length_count_short_50_tbl %>% group_by(exp,gene) %>% summarise(sum_count=sum(count)) %>% arrange(sum_count)
min(gene_utr_length_count_short_50_tbl$count)


apa_n_50_tbl %>% 
  arrange(desc(APA_n)) %>% 
  mutate(exp = factor(exp, levels = unique(exp), ordered = T)) %>% 
  ggplot(aes(x=exp, y=APA_n)) +
  geom_bar(stat="identity", position="dodge") +
  theme(axis.text.x = element_text(angle = 90))


library(vegan)

# Diversity indexes
sum_count_shannon_index <- gene_utr_length_count_short_50_tbl %>% 
  group_by(exp,gene) %>% 
  summarise(sum_count=sum(count),
            shannone_index = diversity(count,"shannon"),
            PA_n = n(),
            evenness = shannone_index/log(PA_n)) %>%
  arrange(sum_count) %>% 
  ungroup()

# short check

etbl <-  filter(sum_count_shannon_index, evenness > 0.8)
lapply(1:100, function(i) {
  return(filter(gene_utr_length_count_short_50_tbl, gene==etbl$gene[i], exp==etbl$exp[i]))
})

sum_count_shannon_index %>% 
  filter(exp != "CTGT") %>% 
  ggplot(aes(x=evenness, color=exp)) +
  geom_density()

temp <- sum_count_shannon_index %>% 
  mutate(q = cut(sum_count, quantile(sum_count, breaks = c(0,0.25,0.5,0.75,1)), labels=c("Q_0.25", "Q_0.5", "Q_0.75", "Q_1"), include.lowest=T))
temp %>% 
  filter(exp != "CTGT") %>% 
  ggplot(aes(x=evenness, color=exp)) +
  geom_density() +
  facet_wrap(.~q)

# rank correlation (method="spearman")
cor_count_Sindex_tbl  <- sum_count_shannon_index %>% 
  group_by(exp) %>% 
  summarise(cor_log2count_Sindex = cor.test(log2(sum_count), shannone_index, method="spearman", exact=F)$estimate,
            cor_log2count_Sindex_pvalue = cor.test(log2(sum_count), shannone_index, method="spearman",exact=F)$p.value)

cor_count_Sindex_tbl %>% 
  filter(exp != "CTGT") %>% 
  arrange(cor_log2count_Sindex) %>% 
  mutate(exp = factor(exp, levels = exp, ordered = T)) %>% 
  ggplot(aes(x=exp,y=cor_log2count_Sindex, color=-log10(cor_log2count_Sindex_pvalue))) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))

# evenness vs gene count
sum_count_shannon_index %>% 
  filter(exp != "CTGT") %>% 
  ggplot(aes(x=log2(sum_count), y=evenness, color=exp)) +
  geom_smooth(se = FALSE)

## Fig a XXX. ------------------------------------

p1 <- sum_count_shannon_index %>% 
  filter(exp != "CTGT") %>% 
  ggplot(aes(x=log2(sum_count), y=evenness)) +
  geom_smooth(aes(group=exp), se = TRUE, color='black', size=0.2, alpha=0.1) +
  geom_smooth(se = FALSE, color='red', size=0.2) +
  theme_bw() +
  ylab('Evenness') +
  xlab('log2 gene count (CPM)') +
  theme(axis.text = element_text(size = 5, color = 'black'),
        axis.text.x = element_text(vjust = 0.5),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(vjust = -1),
        axis.title.x = element_text(vjust = +2),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.2),
        axis.ticks = element_line(size = 0.3))

ggsave(filename = 'evennessVsGeneCount.pdf', plot = p1, height = 45, width = 45, units = 'mm')


# evenness vs gene count BS_0h as an example
filter(sum_count_shannon_index, exp == "BS_0h") %>%
  ggplot(aes(x=log2(sum_count), y=evenness)) +
  geom_point(size=0.5, alpha=0.5)

# rank correlation (method="spearman")
cor_count_evenness_tbl  <- sum_count_shannon_index %>% 
  group_by(exp) %>% 
  summarise(cor_log2count_evenness = cor.test(log2(sum_count), evenness, method="spearman", exact=F)$estimate,
            cor_log2count_evenness_pvalue = cor.test(log2(sum_count), evenness, method="spearman", exact=F)$p.value)

cor_count_evenness_tbl %>% 
  filter(exp != "CTGT") %>%
  arrange(cor_log2count_evenness) %>% 
  mutate(exp = factor(exp, levels = exp, ordered = T)) %>% 
  ggplot(aes(x=exp,y=cor_log2count_evenness, color=-log10(cor_log2count_evenness_pvalue))) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))

## Fig b XXX. ------------------------------------

cor_count_evenness_tbl %>% 
  filter(exp != "CTGT") %>% 
  ggplot(aes(x='cor', y=cor_log2count_evenness)) + 
  geom_violin()


cor_count_evenness_tbl %>% 
  filter(exp != "CTGT") %>% 
  mutate(minlog10pvalue = -log10(cor_log2count_evenness_pvalue)) %>% 
  ggplot(aes(x='-log10pvalue', y=minlog10pvalue)) + 
  geom_violin()

p2 <- cor_count_evenness_tbl %>% 
  dplyr::rename('spearman_rank_cor' = 'cor_log2count_evenness') %>% 
  filter(exp != "CTGT") %>% 
  mutate(minlog10pvalue = -log10(cor_log2count_evenness_pvalue)) %>%
  dplyr::select(-cor_log2count_evenness_pvalue) %>% 
  gather(key = 'groups', value = 'values', -exp) %>% 
  mutate(groups = factor(groups, levels=c('spearman_rank_cor', 'minlog10pvalue'), ordered = TRUE)) %>% 
  ggplot(aes(x='x', y=values)) + 
  geom_violin(fill='gray50', lwd=0.3) +
  facet_wrap(.~groups, scales = 'free_y') +
  theme_bw() +
  theme(axis.text = element_text(size = 5, color = 'black'),
        axis.text.x = element_text(vjust = 0.5),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(vjust = -1),
        axis.title.x = element_text(vjust = +2),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.2),
        axis.ticks = element_line(size = 0.3),
        strip.text = element_blank()
        )

ggsave(filename = 'evennessVsGeneCountCorrelation.pdf', plot = p2, height = 45, width = 45, units = 'mm')



# PA n vs gene count
sum_count_shannon_index %>% 
  filter(exp != "CTGT") %>% 
ggplot(aes(x=log2(sum_count), y=PA_n, color=exp)) +
  geom_smooth(se = FALSE)

cor_count_PAn_tbl  <- sum_count_shannon_index %>% 
  group_by(exp) %>% 
  summarise(cor_log2count_PAn = cor.test(log2(sum_count), PA_n, method="spearman", exact=F)$estimate,
            cor_log2count_PAn_pvalue = cor.test(log2(sum_count), PA_n, method="spearman", exact=F)$p.value)

cor_count_PAn_tbl %>% 
  filter(exp != "CTGT") %>%
  arrange(cor_log2count_PAn) %>% 
  mutate(exp = factor(exp, levels = exp, ordered = T)) %>% 
  ggplot(aes(x=exp,y=cor_log2count_PAn, color=-log10(cor_log2count_PAn_pvalue))) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))

pac_normPoolShortAPAindex_df

length(unique(pac_normPoolShortAPAindex_df$PA))
length(unique(pac_normPoolShortAPAindex_df$gene))

pac_normPoolShortAPAindex_df %>%
  group_by(exp, gene) %>%
  mutate(pcent=count/sum(count), gene_count=sum(count)) %>%
  filter(count >= 10) %>%
  mutate(PA_rank=rank(-pcent, ties.method = "first")) %>%
  ungroup() -> APA_3UTR_countPcentRankAPAindex

table(filter(APA_3UTR_countPcentRankAPAindex, !APA_index)$PA_rank)
APA_3UTR_countPcentRankAPAindexF <- filter(APA_3UTR_countPcentRankAPAindex, APA_index)
nonAPA_3UTR_countPcentRankAPAindexF <- filter(APA_3UTR_countPcentRankAPAindex, !APA_index)

min(APA_3UTR_countPcentRankAPAindexF$count)
min(dplyr::count(APA_3UTR_countPcentRankAPAindexF, exp, gene)$n)

# when is the Rank_1 used -----------------------------------------------------------------------------------

APA_3UTR_countPcentRankAPAindexF %>%
  filter(PA_rank == 1) %>% 
  dplyr::count(gene, PA) %>%
  group_by(gene) %>% 
  mutate(p = n/sum(n)) %>%
  ungroup() -> APA_only_rank1APAindex

APA_only_rank1APAindex %>% group_by(gene) %>% summarise(n=n()) %>% .$n %>% max()

APA_only_rank1APAindex %>% 
  group_by(gene) %>% 
  arrange(desc(p)) %>% 
  mutate(p_order = seq_len(n())) %>% 
  ungroup() -> APA_only_rank1APAindex_p1

filter(APA_only_rank1APAindex_p1, gene=='CopciAB_1000046')

APA_only_rank1APAindex_p1 %>%
  filter(p_order == 1) %>% 
  mutate(n_group = cut(n, c(1,10,20,30,40,50,60,71),include.lowest = TRUE, right = TRUE)) -> APA_only_rank1APAindexFCut_p1 

APA_only_rank1APAindexFCut_p1 %>% 
  group_by(n_group) %>% 
  summarise(n_group_size = n(),
            p = median(p)) -> violinLabel

APA_only_rank1APAindexFCut_p1 %>%   
  ggplot(aes(y=p, x=n_group, fill=n_group)) +
  geom_violin() +
  stat_summary(fun=median, geom="point", size=2, color="red") +
  geom_text(data = violinLabel, aes(x=n_group, y=1.05, label=n_group_size)) +
  theme_bw()

## Fig xxx ---------------------------------------------------------------------------------------------------

p3 <- APA_only_rank1APAindexFCut_p1 %>% 
  filter(n >= 35) %>% 
  ggplot(aes(y=p, x='temp')) +
  geom_violin(fill='gray50', lwd=0.3) +
  theme_bw() +
  theme(axis.text = element_text(size = 5, color = 'black'),
        axis.text.x = element_text(vjust = 0.5),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(vjust = -1),
        axis.title.x = element_text(vjust = +2),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.2),
        axis.ticks = element_line(size = 0.3),
        strip.text = element_blank()
  )

ggsave(filename = 'PACrank1UseRatio.pdf', plot = p3, height = 45, width = 20, units = 'mm')


## cor PAC count ratio vs PAC count / rank groups ------------------------------------

APA_3UTR_countPcentRankAPAindexF %>% 
  filter(PA_rank <=4) %>% 
  filter(exp!="CTGT") %>% 
  ggplot(aes(x=log2(count), y=pcent, color=exp)) +
  geom_smooth(se=FALSE) +
  facet_wrap(.~PA_rank)

# PAC count ratio vs gene count
APA_3UTR_countPcentRankAPAindexF %>% 
  filter(PA_rank <=4) %>% 
  filter(exp!="CTGT") %>% 
  ggplot(aes(x=log2(gene_count), y=pcent, color=exp)) +
  geom_smooth(se=FALSE) +
  facet_wrap(.~PA_rank)

APA_3UTR_countPcentRankAPAindexF %>% 
  filter(PA_rank <=4) %>% 
  filter(exp!="CTGT") %>% 
  group_by(PA_rank, exp) %>% 
  summarise(max_gene_count = max(gene_count)) %>% 
  arrange(exp, PA_rank) %>% 
  View()


## Fig. xxx ------------------------------------------------------------------------------------

p4 <- APA_3UTR_countPcentRankAPAindexF %>% 
  filter(PA_rank <=4) %>% 
  filter(exp!="CTGT") %>% 
  ggplot(aes(x=log2(gene_count), y=pcent)) +
  geom_smooth(aes(group=exp), se=TRUE, color = 'black', size=0.2, alpha=0.1) +
  facet_wrap(.~PA_rank) +
  geom_smooth(se = FALSE, color='red', size=0.2) +
  theme_bw() +
  ylab('Percent of PAC usage') +
  xlab('log2 gene count (CPM)') +
  theme(axis.text = element_text(size = 5, color = 'black'),
        axis.text.x = element_text(vjust = 0.5),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(vjust = -1),
        axis.title.x = element_text(vjust = +2),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.2),
        axis.ticks = element_line(size = 0.3),
        strip.text = element_blank())

ggsave(filename = 'percentOfPACusageVSgeneCountCPM.pdf', plot = p4, height = 45, width = 45, units = 'mm')

# gene expression / PAC usage
cor_rank_countPcentAPAindexF <- APA_3UTR_countPcentRankAPAindexF %>% 
  filter(PA_rank <=4) %>% 
  group_by(exp, PA_rank) %>% 
  summarise(n=n(),
            cor_pcentCount = cor.test(gene_count, pcent, method="spearman", exact=F)$estimate,
            corPvalue_pcentCount = cor.test(gene_count, pcent, method="spearman", exact=F)$p.value)


library(ggbeeswarm)

cor_rank_countPcentAPAindexF %>% 
  group_by(PA_rank) %>% 
  summarise(n_group_size_mean = round(mean(n),digits = 2),
            cor_pcentCount_mean = mean(cor_pcentCount)) -> violinLabel

cor_rank_countPcentAPAindexF %>% 
  ggplot(aes(x=as.character(PA_rank), y=cor_pcentCount)) +
  geom_violin() +
  geom_beeswarm(aes(color=exp, fill=exp)) +
  geom_text(data = violinLabel, aes(x=PA_rank, y=0.4, label=n_group_size_mean))

# extra data for Fig. xxx
cor_rank_countPcentAPAindexF %>% 
  mutate(minlog10pvalue = -log10(corPvalue_pcentCount)) %>% 
  group_by(PA_rank) %>% 
  summarise(mean_rho = mean(cor_pcentCount),
            mean_minlog10pvalue = mean(minlog10pvalue))


cor_rank_countPcentAPAindexF %>% 
  filter(exp!="CTGT") %>% 
  ggplot(aes(x=exp, y=cor_pcentCount, group=PA_rank ,color=as.factor(PA_rank))) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))


# average UTR length/ranks ------------------------------------------------------

APA_3UTR_countPcentRankAPAindexF %>% 
  filter(PA_rank<=4) %>% 
  group_by(exp, PA_rank) %>% 
  summarise(mean_three_UTR_length = mean(three_UTR_length),
            n_gene = n()) %>% 
  ungroup() -> APA_3UTR_countPcentRankAPAindexF_lengthSummary

APA_3UTR_countPcentRankAPAindexF_lengthSummary %>% 
  ggplot(aes(x=exp, y=mean_three_UTR_length, group=PA_rank, color=as.factor(PA_rank))) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))


APA_3UTR_countPcentRankAPAindexF_lengthSummary %>% 
  group_by(PA_rank) %>% 
  summarise(n_group_mean_gene_n = round(mean(n_gene),digits = 2)) -> violinLabel

APA_3UTR_countPcentRankAPAindexF_lengthSummary %>% 
  ggplot(aes(x=as.character(PA_rank), y=mean_three_UTR_length)) +
  geom_violin() +
  geom_beeswarm(aes(color=exp, fill=exp)) +
  geom_text(data = violinLabel, aes(x=PA_rank, y=455, label=n_group_mean_gene_n))


APA_3UTR_countPcentRankAPAindexF %>% 
  filter(PA_rank<=4) %>% 
  ggplot(aes(x= as.character(PA_rank), y=three_UTR_length)) + 
  geom_boxplot()


nonAPA_3UTR_countPcentRankAPAindexF %>% 
  group_by(exp, PA_rank) %>% 
  summarise(mean_three_UTR_length = mean(three_UTR_length),
            n_gene = n()) %>% 
  ungroup() %>% 
  mutate(PA_rank=0)-> nonAPA_3UTR_countPcentRankAPAindexF_lengthSummary

arrange(nonAPA_3UTR_countPcentRankAPAindexF, three_UTR_length)
nonAPA_3UTR_countPcentRankAPAindexF %>% 
  ggplot(aes(x=log2(gene_count), y=three_UTR_length)) +
  geom_point()

arrange(nonAPA_3UTR_countPcentRankAPAindexF, three_UTR_length)

dplyr::count(nonAPA_3UTR_countPcentRankAPAindexF, gene) %>% arrange(desc(n)) -> nonAPA_3UTR_countPcentRankAPAindexF_GeneCount
hist(nonAPA_3UTR_countPcentRankAPAindexF_GeneCount$n)

nonAPA_expNumbFilteredGenes <- filter(nonAPA_3UTR_countPcentRankAPAindexF_GeneCount, n>= 35)$gene

nonAPA_3UTR_countPcentRankAPAindexF %>% 
  filter(gene %in% nonAPA_expNumbFilteredGenes) %>% 
  group_by(exp, PA_rank) %>%
  summarise(mean_three_UTR_length = mean(three_UTR_length),
            n_gene = n()) %>% 
  ungroup() %>% 
  mutate(PA_rank=0)-> nonAPA_3UTR_countPcentRankAPAindexF_expNumbFlengthSummary


bind_rows(nonAPA_3UTR_countPcentRankAPAindexF_lengthSummary,APA_3UTR_countPcentRankAPAindexF_lengthSummary) %>% 
  ggplot(aes(x=exp, y=mean_three_UTR_length, group=PA_rank, color=as.factor(PA_rank))) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))


# violin plot
bind_rows(nonAPA_3UTR_countPcentRankAPAindexF_lengthSummary,APA_3UTR_countPcentRankAPAindexF_lengthSummary) %>% 
  group_by(PA_rank) %>% 
  summarise(n_group_mean_gene_n = round(mean(n_gene),digits = 2)) %>% 
  mutate(PA_rank = as.character(PA_rank)) -> violinLabel


bind_rows(nonAPA_3UTR_countPcentRankAPAindexF_lengthSummary,APA_3UTR_countPcentRankAPAindexF_lengthSummary) %>% 
  ggplot(aes(x=as.character(PA_rank), y=mean_three_UTR_length)) +
  geom_violin() +
  geom_beeswarm(aes(color=exp, fill=exp)) +
  geom_text(data = violinLabel, aes(x=PA_rank, y=455, label=n_group_mean_gene_n))

# Fig. xxx ----------------------------------------------------------------------------

p5 <- bind_rows(nonAPA_3UTR_countPcentRankAPAindexF_lengthSummary, APA_3UTR_countPcentRankAPAindexF_lengthSummary) %>% 
  ggplot(aes(x=as.character(PA_rank), y=mean_three_UTR_length)) +
  geom_violin(fill='gray50', lwd=0.3) +
  geom_text(data = violinLabel, aes(x=PA_rank, y=455, label=n_group_mean_gene_n)) + 
  xlab('PAC rank') +
  ylab("Mean 3'UTR length") +
  theme_bw() +
  theme(axis.text = element_text(size = 5, color = 'black'),
        axis.text.x = element_text(vjust = 0.5),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(vjust = -1),
        axis.title.x = element_text(vjust = +2),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.2),
        axis.ticks = element_line(size = 0.3))

ggsave(filename = 'mean3UTRlengthVSPACrank.pdf', plot = p5, height = 45, width = 45, units = 'mm')
