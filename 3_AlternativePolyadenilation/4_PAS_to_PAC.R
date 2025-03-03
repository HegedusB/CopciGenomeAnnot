# load libraries ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(rtracklayer)
  library(dplyr)
  library(Biostrings)
  library(GenomicRanges)
  library(BSgenome)
  library(tibble)
  library(tidyr)
  library(docstring)
  library(ggplot2)
})


# set wd and screen parameters
work_description <- 'Bind PAS to PAC'
options(width = 300)

# load data -----------------------------------------------------------------------------------------

## filtered pas -------------------------------------------------------------------------------------

pas_f_df <- as.data.frame(import("All_PASs_f.sorted_filtered.bed", format = "bed"))
pas_r_df <- as.data.frame(import("All_PASs_r.sorted_filtered.bed", format = "bed"))
pas_all_df <- bind_rows(pas_r_df, pas_f_df)
pas_all_df <- dplyr::rename(pas_all_df, read_support=name)
pas_all_df <- dplyr::select(pas_all_df, -width)
pas_all_df <- mutate(pas_all_df, read_support = as.integer(read_support), seqnames = as.character(seqnames), strand = as.character(strand))
pas_all_gr <- makeGRangesFromDataFrame(pas_all_df, keep.extra.columns = TRUE)

pas_all_df %>% group_by(strand) %>% summarise(n=n()) # short n summary

## clustered pas (pac) ----------------------------------------------------------------------------------------

pac_f_df <- as.data.frame(import("All_PASs_f.sorted_filtered_merged.bed", format = "bed"))
pac_r_df <- as.data.frame(import("All_PASs_r.sorted_filtered_merged.bed", format = "bed"))

pac_all_df <- bind_rows(pac_r_df, pac_f_df)
pac_all_df <- pac_all_df %>% dplyr::select(-width) %>% dplyr::rename(read_support=name)
pac_all_df <- mutate(pac_all_df, read_support = as.integer(read_support), seqnames = as.character(seqnames), strand = as.character(strand))
pac_all_gr <- makeGRangesFromDataFrame(pac_all_df, keep.extra.columns = TRUE)

pac_all_df %>% group_by(strand) %>% summarise(n=n()) # number of PACs by strand
head(pac_all_df) # number of all PACs

# a pac top pas identify position  -----------------------------------------------------------------------------------------------
# representative PAS, alternative PAS

pas_pac_olap <- findOverlaps(pas_all_gr, pac_all_gr) 

pas_all_df[queryHits(pas_pac_olap),"pac_index"] <- subjectHits(pas_pac_olap)

pas_all_df %>% 
  group_by(pac_index) %>% 
  summarise(sum_read_support = sum(read_support)) -> pas_all_readsupp_summary
identical(pas_all_readsupp_summary$sum_read_support, pac_all_df$read_support)


# Number of PAS per PAC
pas_all_lt <- split(pas_all_df, pas_all_df$pac_index)

table(sapply(pas_all_lt, nrow)) 

PASperPAC <- table(sapply(pas_all_lt, nrow)) %>% as.data.frame() %>% mutate(Var1= as.integer(Var1))
PASperPAC$group <- cut(PASperPAC$Var1, c(1:10,1000), include.lowest=TRUE, right=FALSE)
ggplot(PASperPAC, aes(x=group, y=Freq)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  xlab('number of PAS sites in a PAC cluster')

summary(sapply(pas_all_lt, nrow))


# representative PAS
# pas with max read support
pac_max_pas_lt <- mclapply(pas_all_lt, function(x) {
  result <- x[x$read_support == max(x$read_support),][1,]
  result$cluster_n <- nrow(x)
  return(result)
}, mc.cores = 10)

pac_max_pas_df <- bind_rows(pac_max_pas_lt)
table(pac_max_pas_df$end - pac_max_pas_df$start) 

summary(pac_max_pas_df$read_support)
boxplot(pac_max_pas_df$read_support, breaks = 1000)
boxplot(pac_max_pas_df$read_support, breaks = 1000, ylim=c(0, 2000))

pac_all_df_p <- pac_all_df
pac_all_df_p[pac_max_pas_df$pac_index, "pas_max_read_support"] <- pac_max_pas_df$read_support
pac_all_df_p[pac_max_pas_df$pac_index, "cluster_n"] <- pac_max_pas_df$cluster_n
pac_all_df_p[pac_max_pas_df$pac_index, "pas_start"] <- pac_max_pas_df$start
pac_all_df_p[pac_max_pas_df$pac_index, "pas_end"] <- pac_max_pas_df$end

pac_all_df_p <- arrange(pac_all_df_p, nchar(seqnames), seqnames, start)

pac_all_df_p_gr <- makeGRangesFromDataFrame(pac_all_df_p, keep.extra.columns = TRUE)

export(pac_all_df_p_gr, "pac_with_max_pas_readsupport.gtf") 


# max read support plot
boxplot(pac_all_df_p$pas_max_read_support)
summary(pac_all_df_p$pas_max_read_support)
summary(pac_all_df_p$read_support)
summary(pac_all_df_p$cluster_n)

max_pas_df <- dplyr::select(pac_all_df_p, -start, -end) %>% dplyr::rename(start = pas_start, end = pas_end)
max_pas_df$row_index <- seq_len(nrow(max_pas_df))
max_pas_gr <- makeGRangesFromDataFrame(max_pas_df, keep.extra.columns = TRUE)

# bind PAC to 3'UTR border  -----------------------------------------------------------------------------------------------------

annot_gp <- read_tsv("/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/clean_annotation_20220427_dir/CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered.genePred", 
                     col_name=F, col_type=cols(X9 = col_character(), X10 = col_character()))

length(unique(str_remove(annot_gp$X1, pattern='\\.T[0-9]+$')))

annot_gp <- filter(annot_gp, str_detect(X1, pattern="T0$")) 

utr_3_border <- mutate(annot_gp, utr3_border = ifelse(X3 == "+", X5, X4+1))
utr_3_border <- dplyr::select(utr_3_border, 1:3, utr3_border)
colnames(utr_3_border) <- c("transcript_id", "seqnames", "strand", "start")
utr_3_border$end <- utr_3_border$start
utr_3_border <- arrange(utr_3_border, nchar(seqnames), seqnames, start)
utr_3_border <- makeGRangesFromDataFrame(utr_3_border, keep.extra.columns = T)

## distance from 3'UTR ------------------------------------------------------------------
dtnearest <- distanceToNearest(x = max_pas_gr, subject = utr_3_border)
max_pas_noerror_gr <- max_pas_gr[queryHits(dtnearest)]


max_pas_noerror_df <- as.data.frame(max_pas_noerror_gr)
max_pas_noerror_df$transcript_id_nn <- mcols(utr_3_border)[subjectHits(dtnearest), "transcript_id"]
max_pas_noerror_df$distance_nn <- mcols(dtnearest)$distance

max_pas_noerror_df$overlap <- FALSE

foverlaps <- findOverlaps(query = max_pas_noerror_gr, subject = utr_3_border)
max_pas_noerror_df$overlap[queryHits(foverlaps)] <- TRUE
max_pas_noerror_df <- mutate(max_pas_noerror_df, distance_nn_correct = ifelse(overlap, distance_nn, distance_nn + 1))

pac_all_df_p[max_pas_noerror_df$row_index, "transcript_id_nn"] <- max_pas_noerror_df$transcript_id_nn
pac_all_df_p[max_pas_noerror_df$row_index, "distance_nn"] <- max_pas_noerror_df$distance_nn
pac_all_df_p[max_pas_noerror_df$row_index, "overlap"] <- max_pas_noerror_df$overlap
pac_all_df_p[max_pas_noerror_df$row_index, "distance_nn_correct"] <- max_pas_noerror_df$distance_nn_correct

pac_all_df_p <- mutate(pac_all_df_p, name = str_c(transcript_id_nn, "_d", distance_nn_correct, "_max", pas_max_read_support))

pac_all_df_p_gr <- makeGRangesFromDataFrame(pac_all_df_p, keep.extra.columns = TRUE)

export(pac_all_df_p_gr, "pac_with_max_pas_readsupport_annot.gtf")

hist(pac_all_df_p$distance_nn_correct, breaks = 1000, main = "", xlab="pac max pas position distance from TTS")

summary(pac_all_df_p$distance_nn_correct)
summary(pac_all_df_p$pas_max_read_support)

distance_problems <- pac_all_df_p[pac_all_df_p$distance_nn >= 1000,]

dim(distance_problems)
dim(pac_all_df_p)

# check readsupport  ----------------------------------------------------

pas_readsupport_path <- list.files(".", pattern = "_Hisat2-mapping_3prime-ends.sorted.bed", full.names = FALSE)
pas_readsupport_names <- str_remove(pas_readsupport_path, pattern = "_Hisat2-mapping_3prime-ends.sorted.bed")

bed_col_names <- c("seqnames", "start", "end", "count", "utr_tag", "strand")
pas_readsupport_gr_lt <- lapply(pas_readsupport_path, function(x) {
  bed_file <- read_tsv(x, col_names = bed_col_names)
  bed_file <- bed_file %>% mutate(start = start + 1) %>% dplyr::select(-utr_tag)
  bed_file <- makeGRangesFromDataFrame(bed_file, keep.extra.columns = TRUE)
  return(bed_file)
})


pac_pas_olap_sum_lt <- mclapply(pas_readsupport_gr_lt, function(x) {
  pac_pas_olap <- findOverlaps(x, pac_all_df_p_gr)
  pas_olap <- x[queryHits(pac_pas_olap),]
  pac_olap <- pac_all_df_p_gr[subjectHits(pac_pas_olap)]

  pas_olap_lt <- split(as.data.frame(pas_olap), subjectHits(pac_pas_olap))
  pas_olap_countsum <- sapply(pas_olap_lt, function(x) {return(sum(x$count))})
  pas_olap_countsum <- as.data.frame(pas_olap_countsum) %>% rownames_to_column(var = "pac_index")
  return(pas_olap_countsum)
}, mc.cores = 10)

names(pac_pas_olap_sum_lt) <- pas_readsupport_names
pac_pas_olap_sum_tbl <- bind_rows(pac_pas_olap_sum_lt, .id = "sample")
pac_pas_olap_sum_tbl <- tidyr::spread(pac_pas_olap_sum_tbl, key=sample, value=pas_olap_countsum)
pac_pas_olap_sum_tbl <- arrange(pac_pas_olap_sum_tbl, as.integer(pac_index))
pac_all_df_p2_tbl <- as.data.frame(pac_all_df_p_gr)
pac_all_df_p2_tbl <- bind_cols(pac_all_df_p2_tbl, pac_pas_olap_sum_tbl)


count_check_tbl <- tibble(read_support_pac = pac_all_df_p2_tbl$read_support,
                          read_support_sumpersample = apply(pac_all_df_p2_tbl[-(1:17)],1, function(x) sum(x, na.rm = T)))

(sum(count_check_tbl$read_support_pac == count_check_tbl$read_support_sumpersample))/nrow(count_check_tbl)


plot(count_check_tbl$read_support_pac,count_check_tbl$read_support_sumpersample)
cor(count_check_tbl$read_support_pac,count_check_tbl$read_support_sumpersample)

plot(log2(count_check_tbl$read_support_pac), log2(count_check_tbl$read_support_sumpersample))
cor(log2(count_check_tbl$read_support_pac), log2(count_check_tbl$read_support_sumpersample))

cor(count_check_tbl$read_support_pac, count_check_tbl$read_support_sumpersample,method = 'spearman')

abselteres <- abs(count_check_tbl$read_support_pac - count_check_tbl$read_support_sumpersample)
summary(abselteres)
boxplot(abselteres, ylim=c(0,60))

write_tsv(pac_all_df_p2_tbl, "pac_with_max_pas_readsupport_annot_pasreadsupport.tsv")

pac_all_df_p3_tbl <- pac_all_df_p2_tbl
pac_all_df_p3_tbl[18:ncol(pac_all_df_p3_tbl)][is.na(pac_all_df_p3_tbl[18:ncol(pac_all_df_p3_tbl)])] <- 0
write_tsv(pac_all_df_p3_tbl, "pac_with_max_pas_readsupport_annot_pasreadsupport_noNA.tsv")


# Quick check ------------------------------------------------------------------------------

pac_all_df_p_lt <- split(pac_all_df_p, pac_all_df_p$transcript_id_nn)
max_pas_readsupport_lt <- lapply(pac_all_df_p_lt, function(x) return(x[x$read_support == max(x$read_support),][1,])) # itt nem a pas_max_read_support volt nézve !!!!!
max_pas_readsupport <- bind_rows(max_pas_readsupport_lt)
summary(max_pas_readsupport$distance_nn_correct) # viszonylag közel
summary(max_pas_readsupport$pas_max_read_support)

# connect pac with the annotation ------------------------------------------------------------------------------------------------------
utr_3_border 
max_pas_gr 

library(GenomicFeatures)

annot_path <- '/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/clean_annotation_20220427_dir/CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered.gtf'  
annot_gr <- import(annot_path)
annot_gr <- annot_gr[str_detect(annot_gr$gene_id, pattern = 'T0$')]


annot_TxDb <- makeTxDbFromGRanges(annot_gr)
annot_genes_gr <- genes(annot_TxDb)
annot_exon_gr <- exonsBy(annot_TxDb, by = "tx", use.names=TRUE)
annot_cds_gr <- cdsBy(annot_TxDb, by = "tx", use.names=TRUE) 
annot_introns_gr <- intronsByTranscript(annot_TxDb, use.names=TRUE)
annot_5UTR_gr <- fiveUTRsByTranscript(annot_TxDb, use.names=TRUE)
annot_3UTR_gr <- threeUTRsByTranscript(annot_TxDb, use.names=TRUE)

pacAnnotEntry_fct <- function(x) {
 

  pacCollector_lt  <- list()
  
  pacGeneOlap <- findOverlaps(query = max_pas_gr, subject = annot_genes_gr, ignore.strand=x)
  
  sum(duplicated(queryHits(pacGeneOlap)))
  annot_genes_gr[subjectHits(pacGeneOlap[duplicated(queryHits(pacGeneOlap))])]
  
  (length(unique(queryHits(pacGeneOlap))) / length(max_pas_gr))*100
  
  
  pacCollector_lt$genePAC <- unique(max_pas_gr$row_index[queryHits(pacGeneOlap)])
  pacCollector_lt$intergenePAC <- setdiff(max_pas_gr$row_index, pacCollector_lt$genePAC)
  
  pacExonOlap <- findOverlaps(query = max_pas_gr, subject = annot_exon_gr, ignore.strand=x)
  
  sum(duplicated(queryHits(pacExonOlap)))
  annot_exon_gr[subjectHits(pacExonOlap[duplicated(queryHits(pacExonOlap))])]
  
  (length(unique(queryHits(pacExonOlap))) / length(max_pas_gr))*100
  
  pacCollector_lt$exonPAC <- unique(max_pas_gr$row_index[queryHits(pacExonOlap)])
  
  pacCDSOlap <- findOverlaps(query = max_pas_gr, subject = annot_cds_gr, ignore.strand=x)
  
  sum(duplicated(queryHits(pacCDSOlap)))
  annot_exon_gr[subjectHits(pacCDSOlap[duplicated(queryHits(pacCDSOlap))])]
  
  (length(unique(queryHits(pacCDSOlap))) / length(max_pas_gr))*100
  
  pacCollector_lt$cdsPAC <- unique(max_pas_gr$row_index[queryHits(pacCDSOlap)])
  
  annot_introns_gr
  
  pacIntronOlap <- findOverlaps(query = max_pas_gr, subject = annot_introns_gr, ignore.strand=x)
  
  sum(duplicated(queryHits(pacIntronOlap)))
  annot_exon_gr[subjectHits(pacIntronOlap[duplicated(queryHits(pacIntronOlap))])]
  
  (length(unique(queryHits(pacIntronOlap))) / length(max_pas_gr))*100
  
  pacCollector_lt$intronPAC <- unique(max_pas_gr$row_index[queryHits(pacIntronOlap)])
  
  
  
  temp <- intersect(pacCollector_lt$intronPAC, pacCollector_lt$exonPAC)
  max_pas_gr[temp]
  
  length(unique(c(pacCollector_lt$intronPAC, pacCollector_lt$exonPAC))) / length(max_pas_gr)
  

  annot_5UTR_gr
  
  pac5utrOlap <- findOverlaps(query = max_pas_gr, subject = annot_5UTR_gr, ignore.strand=x)
  
  sum(duplicated(queryHits(pac5utrOlap)))
  annot_exon_gr[subjectHits(pac5utrOlap[duplicated(queryHits(pac5utrOlap))])]
  
  (length(unique(queryHits(pac5utrOlap))) / length(max_pas_gr))*100
  
  pacCollector_lt$utr5PAC <- unique(max_pas_gr$row_index[queryHits(pac5utrOlap)])
  
  annot_3UTR_gr
  
  pac3utrOlap <- findOverlaps(query = max_pas_gr, subject = annot_3UTR_gr, ignore.strand=x)
  
  sum(duplicated(queryHits(pac3utrOlap)))
  annot_exon_gr[subjectHits(pac3utrOlap[duplicated(queryHits(pac3utrOlap))])]
  
  (length(unique(queryHits(pac3utrOlap))) / length(max_pas_gr))*100
  
  pacCollector_lt$utr3PAC <- unique(max_pas_gr$row_index[queryHits(pac3utrOlap)])
  
  temp <- unique(c(pacCollector_lt$intronPAC, pacCollector_lt$cdsPAC, pacCollector_lt$utr3PAC, pacCollector_lt$utr5PAC))
  length(temp) / length(max_pas_gr)
  
  return(pacCollector_lt)
}

pacAnnotEntry_lt <- lapply(c(TRUE, FALSE), FUN = pacAnnotEntry_fct)

names(pacAnnotEntry_lt) <- c('sense_antisense', 'only_sense')

medianMaxPASresult_lt <- lapply(pacAnnotEntry_lt, function(x) {
  meanMaxPAS_lt <- lapply(x, function(y) {
    return(median(max_pas_gr$pas_max_read_support[max_pas_gr$row_index %in% y]))
  })
  medianMaxPAS_tbl <- as_tibble(unlist(meanMaxPAS_lt), rownames='annotEntry')
  return(medianMaxPAS_tbl)
})

medianMaxPASresult_tbl <- bind_rows(medianMaxPASresult_lt, .id='sense_group')


## summary -------------------------------------------------------------------------------------------------------------------------

pacAnnotEntrySummary_lt <- lapply(pacAnnotEntry_lt, function(x) {
  sapply(x, length) %>% 
    as_tibble(rownames = 'annotEntry') %>% 
    mutate(pcent = (value/length(max_pas_gr))*100) -> result
  return(result)
})

pacAnnotEntrySummary_tbl <- bind_rows(pacAnnotEntrySummary_lt, .id='sense_group')

pacAnnotEntrySummaryPlus_tbl <- left_join(pacAnnotEntrySummary_tbl, medianMaxPASresult_tbl,  by=c('sense_group', 'annotEntry'), suffix=c('.count', '.medianReadSupport'))

pacAnnotEntrySummaryPlus_tbl %>% 
  mutate(annotEntry = factor(annotEntry, levels=c('genePAC', 'exonPAC', 'cdsPAC', 'intronPAC', 'utr5PAC', 'utr3PAC', 'intergenePAC'), ordered = TRUE)) -> pacAnnotEntrySummaryPlus_tbl

length(unique(unlist(pacAnnotEntry_lt$sense_antisense)))
length(unique(unlist(pacAnnotEntry_lt$sense_antisense[c('intergenePAC', 'cdsPAC', 'intronPAC', 'utr5PAC', 'utr3PAC')])))
length(unique(unlist(pacAnnotEntry_lt$only_sense)))
length(unique(unlist(pacAnnotEntry_lt$only_sense[c('intergenePAC', 'cdsPAC', 'intronPAC', 'utr5PAC', 'utr3PAC')])))

PACgeneAnnotLocExtend_lt <- list(cdsPAC_sense = pacAnnotEntry_lt$only_sense$cdsPAC,
                                 cdsPAC_antisense = setdiff(pacAnnotEntry_lt$sense_antisense$cdsPAC, pacAnnotEntry_lt$only_sense$cdsPAC),
                                 intronPAC_sense = pacAnnotEntry_lt$only_sense$intronPAC,
                                 intronPAC_antisense = setdiff(pacAnnotEntry_lt$sense_antisense$intronPAC, pacAnnotEntry_lt$only_sense$intronPAC),
                                 utr5PAC_sense = pacAnnotEntry_lt$only_sense$utr5PAC,
                                 utr5PAC_antisense = setdiff(pacAnnotEntry_lt$sense_antisense$utr5PAC, pacAnnotEntry_lt$only_sense$utr5PAC),
                                 utr3PAC_sense = pacAnnotEntry_lt$only_sense$utr3PAC,
                                 utr3PAC_antisense = setdiff(pacAnnotEntry_lt$sense_antisense$utr3PAC, pacAnnotEntry_lt$only_sense$utr3PAC),
                                 intergenePAC = pacAnnotEntry_lt$sense_antisense$intergenePAC)

PACgeneAnnotLocExtend_lt <- lapply(PACgeneAnnotLocExtend_lt, function(x) return(as_tibble(x)))
PACgeneAnnotLocExtend_tbl <- bind_rows(PACgeneAnnotLocExtend_lt, .id='annotEntry') %>% 
  dplyr::rename('pac_index' = 'value') %>% 
  mutate(pac_index = as.character(pac_index))

PACgeneAnnotLocExtend_tbl <- left_join(PACgeneAnnotLocExtend_tbl, pac_all_df_p3_tbl[c('pac_index', 'pas_max_read_support')], by='pac_index')

PACgeneAnnotLocExtend_tbl %>% 
  group_by(pac_index) %>% 
  summarise(annotEntry_sum = str_c(annotEntry, collapse = '-')) %>% 
  mutate(pac_index = as.character(pac_index)) -> PACgeneAnnotLocExtendCollapsed_tbl

setdiff(pac_all_df_p3_tbl$pac_index, PACgeneAnnotLocExtendCollapsed_tbl$pac_index)
setdiff(PACgeneAnnotLocExtendCollapsed_tbl$pac_index, pac_all_df_p3_tbl$pac_index)

pac_all_df_p4_tbl <- left_join(pac_all_df_p3_tbl, PACgeneAnnotLocExtendCollapsed_tbl, by='pac_index') %>% 
  mutate(name = str_c(name, annotEntry_sum, sep = ':')) %>% 
  dplyr::select(-annotEntry_sum)


write_tsv(pac_all_df_p4_tbl, "pac_with_max_pas_readsupport_annot_pasreadsupport_noNA_annotEntry.tsv")

pac_all_df_p4_gr <- makeGRangesFromDataFrame(pac_all_df_p4_tbl, keep.extra.columns = TRUE)
export(pac_all_df_p4_gr, 'pac_with_max_pas_readsupport_annot_pasreadsupport_noNA_annotEntry.gtf')

## Table xxx: All PAC metadata in gtf-----------------------------------------------------------


PACgeneAnnotLocExtend_tbl %>% 
  group_by(annotEntry) %>% 
  summarise(annotEntrySize = n(),
            median_pas_max_read_support = median(pas_max_read_support)) %>% 
  ungroup() %>% 
  mutate(annotEntryPercent = round((annotEntrySize / sum(annotEntrySize))*100, digits = 2)) %>% 
  mutate(annotEntry = factor(annotEntry, levels = c('cdsPAC_sense', 'cdsPAC_antisense',
                                                    'intronPAC_sense','intronPAC_antisense', 
                                                    'utr5PAC_sense', 'utr5PAC_antisense',
                                                    'utr3PAC_sense', 'utr3PAC_antisense',
                                                    'intergenePAC'), ordered = TRUE)) %>% 
  arrange(annotEntry) -> PACgeneAnnotLocExtendSummary_tbl

## Table xxx: All PAC localization summary -----------------------------------------------------------------------

write_tsv(PACgeneAnnotLocExtendSummary_tbl, 'PAC_annot_entry.tsv')


## plot ------------------------------------------------------------------------------------------------------

# plot all

PACgeneAnnotLocExtendSummary_tbl %>% 
  ggplot(aes(x=annotEntry, y=annotEntryPercent, fill=median_pas_max_read_support)) +
  geom_bar(position = 'dodge', stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90))


PACgeneAnnotLocExtendSummary_tbl %>% knitr::kable()

# filter best pac ------------------------------------------------------------------------------------------------------

##  3'utr PAC + Intergen PAC + downstream antisense pac  -------------------------

# 3'UTR: 
utr3PAC_tbl <- tibble(row_index = pacAnnotEntry_lt$only_sense$utr3PAC, 
                      annotEntry = 'utr3PAC')

# intergene:
setdiff(pacAnnotEntry_lt$only_sense$intergenePAC, pacAnnotEntry_lt$sense_antisense$intergenePAC)
setdiff(pacAnnotEntry_lt$sense_antisense$intergenePAC, pacAnnotEntry_lt$only_sense$intergenePAC)

intergenePAC_tbl <- tibble(row_index = unique(c(pacAnnotEntry_lt$sense_antisense$intergenePAC, 
                                                pacAnnotEntry_lt$only_sense$intergenePAC)), 
                           annotEntry = 'intergenePAC')

genePAC_tbl <- bind_rows(utr3PAC_tbl, intergenePAC_tbl) %>% mutate(row_index = as.character(row_index))
pac_all_df_p3Short_tbl <- pac_all_df_p4_tbl[,1:17]

pac_all_archive_tbl <- pac_all_df_p3Short_tbl

utr_3_border_short <- as.data.frame(utr_3_border)[c('start', 'transcript_id')] %>% dplyr::rename('utr3border'='start')
pac_all_df_p3Short_tbl <- left_join(pac_all_df_p3Short_tbl, utr_3_border_short, by=c('transcript_id_nn'='transcript_id'))

setdiff(genePAC_tbl$row_index, pac_all_df_p3Short_tbl$pac_index)

genePACplus_tbl <- inner_join(pac_all_df_p3Short_tbl, genePAC_tbl, by=c('pac_index'='row_index'))
filter(genePACplus_tbl, str_detect(name, pattern='-'))


annot_cdsStrand <- as.data.frame(strand(annot_cds_gr)) %>% group_by(group_name) %>% summarise(strand = unique(value))

cdsRange_tbl <- tibble(transcript_id = names(annot_cds_gr),
                       cds_start = min(start(annot_cds_gr)),
                       cds_end = max(end(annot_cds_gr)))

cdsRange_tbl <- left_join(cdsRange_tbl, annot_cdsStrand, by=c('transcript_id'='group_name'))

cdsRange_tbl$cds_pos <- apply(cdsRange_tbl,1,function(x){
  if(x[['strand']]=="+") {
    return(as.integer(x[['cds_end']]))
  } else {
    return(as.integer(x[['cds_start']]))
  }
})


cdsRange_tbl <- dplyr::select(cdsRange_tbl, -cds_start, -cds_end) %>% rename('cds_strand' = 'strand')

genePACplus_tbl <- left_join(genePACplus_tbl, cdsRange_tbl, by=c('transcript_id_nn'='transcript_id'))
sum(is.na(genePACplus_tbl$transcript_id_nn))
genePACplus_tbl <- genePACplus_tbl[!is.na(genePACplus_tbl$transcript_id_nn),]

genePACplus_lt <- split(genePACplus_tbl,genePACplus_tbl$transcript_id_nn)

overTheGenePAC_lt <- lapply(genePACplus_lt, function(x) {
  if(x[['cds_strand']][[1]] == '+') {
    return(filter(x, pas_start <= x[['cds_pos']]))
  } else {
    return(filter(x, pas_start >= x[['cds_pos']]))
  }
})

table(sapply(overTheGenePAC_lt, nrow)) 
overTheGenePAC_tbl <- bind_rows(overTheGenePAC_lt[sapply(overTheGenePAC_lt, nrow)>=1])

nrow(overTheGenePAC_tbl) 
genePACplus_tbl <- filter(genePACplus_tbl, !pac_index %in% overTheGenePAC_tbl$pac_index)
dim(genePACplus_tbl)

filter(pac_all_df_p3Short_tbl, transcript_id_nn == 'CopciAB_465586.T0')
filter(genePACplus_tbl, pac_index == '855') 

genePACplus_tbl %>% 
  filter(annotEntry == 'intergenePAC') %>% 
  .$name %>% 
  str_split(pattern = ':') %>% 
  sapply(., '[[', 2) %>% 
  str_split(pattern = '_') %>% 
  sapply(., function(x) x == 'sense') %>% 
  sapply(., function(x) any(x)) %>% 
  table()


## 1000 bp filtering --------------------------------------------------------------

genePACplus1000_tbl <- filter(genePACplus_tbl, distance_nn_correct<=1000)
dim(genePACplus1000_tbl) 
length(unique(genePACplus1000_tbl$transcript_id_nn)) 

summary(genePACplus1000_tbl$distance_nn_correct)
boxplot(genePACplus1000_tbl$distance_nn_correct)
summary(genePACplus1000_tbl$pas_max_read_support)
boxplot(genePACplus1000_tbl$pas_max_read_support)

PACnPerGene <- genePACplus1000_tbl %>% group_by(transcript_id_nn) %>% summarise(n = n())
summary(PACnPerGene$n)
boxplot(PACnPerGene$n)

count(PACnPerGene, n) %>% knitr::kable()
count(PACnPerGene, n) %>% 
  ggplot(aes(x=n, y=nn)) +
  geom_bar(stat = 'identity', position = 'dodge')


sum(count(PACnPerGene, n)$nn*count(PACnPerGene, n)$n)

PACnPerGene_nnp <- count(PACnPerGene, n) %>% mutate(nn_p = round((nn/sum(nn))*100, digits = 2))

# Fig. xxx --------------------------------------------------------------------------------------------
PACnPerGene_nnp %>% 
  mutate(n_index = ifelse(n<=10, n, 11)) %>% 
  group_by(n_index) %>% 
  summarise(nn_p = sum(nn_p)) %>% 
  mutate(n_index = as.character(n_index)) %>% 
  mutate(n_index = ifelse(n_index!='11', n_index, '10+')) %>% 
  mutate(n_index = factor(n_index, levels=n_index, ordered = TRUE)) %>% 
  ggplot(aes(x=n_index, y=nn_p)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw() +
  xlab('Number of PACs') +
  ylab('Percent of genes') +
  theme(axis.text = element_text(size = 5, color = 'black'),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(vjust = -1),
        axis.title.x = element_text(vjust = +2),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.2),
        axis.ticks = element_line(size = 0.3))

ggsave('percentOfGenesVsPACnumber_setFontSize.pdf', device = 'pdf', 
       height = 45, width = 45, units = 'mm')

PACnPerGene_nnp %>% knitr::kable()

sum(PACnPerGene_nnp[2:5,'nn_p',drop=TRUE])
sum(PACnPerGene_nnp[-1,'nn_p',drop=TRUE])

## az 1000 bp max filtering --------------------------------------------
genePACplus1000MaxPac_tbl <- group_by(genePACplus1000_tbl, transcript_id_nn) %>% 
  filter(pas_max_read_support == max(pas_max_read_support)) %>% 
  filter(seq_len(n()) == 1) %>% 
  ungroup()

summary(genePACplus1000MaxPac_tbl$distance_nn_correct)
boxplot(genePACplus1000MaxPac_tbl$distance_nn_correct)
summary(genePACplus1000MaxPac_tbl$pas_max_read_support)
boxplot(genePACplus1000MaxPac_tbl$pas_max_read_support)

genePACplus1000MaxPac_tbl$distance_nn_correct %>% 
  cut(c(0:100,1000), include.lowest=TRUE, right=FALSE) %>% 
  table() %>%
  barplot(las=2)

genePACplus1000MaxPac_tbl$distance_nn_correct %>% 
  cut(c(0,1,5,10,20, 50, 100,1000), include.lowest=TRUE, right=FALSE) %>% 
  table() %>% 
  as.data.frame() %>% 
  setNames(c('sizeGroup','n')) %>% 
  mutate(p = round((n/sum(n))*100, digits = 2)) -> genePACplus1000MaxPacDistanceDist_tbl


genePACplus1000MaxPacDistanceDist_tbl %>% 
  ggplot(aes(x=sizeGroup, y=p)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw()


# longread data  ------------------------------------------------------------------------

CopciAB_id_dict <- read_tsv('CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered_T0isoforme_source.tsv')

genePACplus1000MaxPac_tbl %>% 
  mutate(transcript_id_nn = str_remove(transcript_id_nn, pattern = '\\.T0$')) %>% 
  left_join(CopciAB_id_dict, tID_T0, by=c('transcript_id_nn'='tID_T0')) -> genePACplus1000MaxPacSourceInfo_tbl

dim(genePACplus1000MaxPacSourceInfo_tbl)

genePACplus1000MaxPacSourceInfo_tbl %>% 
  mutate(cutIndex = cut(distance_nn_correct, c(0:100,1000), include.lowest=TRUE, right=FALSE)) -> genePACplus1000MaxPacSourceInfo_tbl

split(CopciAB_id_dict$tID_T0 %in% genePACplus1000MaxPacSourceInfo_tbl$transcript_id_nn, CopciAB_id_dict$source_index_short) %>% 
  sapply(., table)


genePACplus1000MaxPacSourceInfo_tbl %>% 
  count(source_index_short, cutIndex) %>% 
  ggplot(aes(x=cutIndex, y=n)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_grid(source_index_short~.) +
  theme(axis.text.x = element_text(angle = 90))

genePACplus1000MaxPacSourceInfo_tbl %>% 
  mutate(cutIndex = cut(distance_nn_correct, c(0,1,5,10,20, 50, 100,1000), include.lowest=TRUE, right=FALSE)) %>% 
  count(source_index_short, cutIndex) %>% 
  mutate(p = round((n/sum(n))*100, digits = 2)) -> genePACplus1000MaxPacSourceInfoDistanceDist_tbl

# Fig. xxx -----------------------------------------------------------------------------------------------------------
color_table <- tibble(source_index_short = c('alternative', 'longread'),
                      color = c('gray50', 'black'))

cols <- c('alternative' = 'black', 'longread' = 'gray50')

genePACplus1000MaxPacSourceInfoDistanceDist_tbl %>% 
  ggplot(aes(x=cutIndex, y=p, fill=source_index_short)) +
  geom_bar(stat = 'identity') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(labels=c('overlap', 
                            'from 1 to 4', 
                            'from 5 to 9', 
                            'from 10 to 19', 
                            'from 20 to 49', 
                            'from 50 to 99', 
                            'from 100 to 1000')) +
  ylab('Percent of genes') +
  xlab('Deviation from the TTS (nt)') +
  theme(axis.text = element_text(size = 5, color = 'black'),
        axis.text.x = element_text(angle = 50, vjust = 0.5),
        axis.title = element_text(size = 6),
        axis.title.y = element_text(vjust = -1),
        axis.title.x = element_text(vjust = +2),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size=0.2),
        axis.ticks = element_line(size = 0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.key.size = unit(2, 'mm'),
        legend.position=c(0.8,1),
        legend.box = "horizontal",
        legend.box.margin = margin(0, 0, 0, 12))

ggsave('percentOfGenesVsDistanceFromTTS_setFontSize.pdf', device = 'pdf', 
       height = 45, width = 45, units = 'mm')


View(filter(genePACplus1000MaxPacSourceInfo_tbl, source_index_short == 'alternative', cutIndex == '[100,1e+03]'))
View(filter(genePACplus1000MaxPacSourceInfo_tbl, source_index_short == 'longread', cutIndex == '[100,1e+03]'))

filter(genePACplus1000MaxPacSourceInfo_tbl, transcript_id_nn == 'CopciAB_493959')

# save environment ---------------------------------------------------------------------------------------------------

save.image(file='bind_PAS_to_PAC_230607_Environment.RData')

genePACplus1000final_tbl <- genePACplus1000_tbl
genePACplus1000final_tbl

genePACplus1000final_tbl %>% 
  mutate(transcript_id_nn = str_remove(transcript_id_nn, pattern = '\\.T0$')) %>% 
  left_join(CopciAB_id_dict, tID_T0, by=c('transcript_id_nn'='tID_T0')) -> genePACplus1000final_tbl

genePACplus1000final_tbl$maxPACindex <- genePACplus1000final_tbl$pac_index %in% genePACplus1000MaxPacSourceInfo_tbl$pac_index
genePACplus1000final_tbl$transcript_id_nn[genePACplus1000final_tbl$maxPACindex]

pac_all_df_p4_tbl[17:ncol(pac_all_df_p3_tbl)]

genePACplus1000final_tbl <- left_join(genePACplus1000final_tbl, pac_all_df_p3_tbl[17:ncol(pac_all_df_p3_tbl)], by='pac_index')

write_tsv(genePACplus1000final_tbl, 'pac_with_max_pas_readsupport_annot_pasreadsupport_noNA_annotEntry_FilteredRelevantPAC.tsv')


temp <- PACgeneAnnotLocExtend_tbl[PACgeneAnnotLocExtend_tbl$pac_index %in% genePAC_tbl$row_index,] # genePAC_tbl # pac helyek: 3UTR + kiterjesztett intergen
table(temp$annotEntry)


View(filter(genePACplus1000final_tbl, pac_index %in% filter(temp, str_detect(annotEntry, pattern= 'cdsPAC_sense'))$pac_index)[,1:18])
