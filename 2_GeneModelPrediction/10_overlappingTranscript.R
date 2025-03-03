# load libraries

suppressPackageStartupMessages({
  library(zoo)
  library(Biostrings)
  library(GenomicFeatures)
  library(rtracklayer)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(BSgenome)
})


# load data ------------------------------------------------------

iFormSource_path <- '/data/bhegedus/Projects/pj8_3genseq/Transcript_sequencing/Coprinopsis_cinerea/results/refannotation_summary/CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered_T0isoforme_source.tsv'
iFormSource_tbl <- read_tsv(iFormSource_path)

table(iFormSource_tbl$source_index_short)
iFormSourceLongreadOnly_tbl <- filter(iFormSource_tbl, source_index_short=="longread")


annot_path <- "../remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/clean_annotation_20220427_dir/CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered.gtf"
annot_gr <- import(annot_path)
annot_t0_gr <- annot_gr[str_detect(annot_gr$gene_id, pattern = "T0$")]
annot_t0_df <- as.data.frame(annot_t0_gr)

length(unique(annot_t0_df$gene_id))


annot_gr$gene_id[str_detect(annot_gr$gene_id, pattern = 'T0$')] %>% unique() %>% length()

annot_t0_df %>% filter(type %in% c('5UTR', '3UTR')) %>% select(gene_id, type) -> annot_t0_UTR_df
annot_t0_UTR_lt <- split(annot_t0_UTR_df$gene_id, as.character(annot_t0_UTR_df$type))
annot_t0_UTR_lt <- lapply(annot_t0_UTR_lt, unique)
sapply(annot_t0_UTR_lt, length)
annot_t0_bothUTR <- intersect(annot_t0_UTR_lt$`3UTR`, annot_t0_UTR_lt$`5UTR`)
length(annot_t0_bothUTR)
annot_t0_bothUTR_longreadOnly <- annot_t0_bothUTR[str_remove(annot_t0_bothUTR, pattern = '\\.T0$') %in% iFormSourceLongreadOnly_tbl$tID_T0]
length(annot_t0_bothUTR_longreadOnly)

annot_t0_lt <- split(annot_t0_df,annot_t0_df$seqnames)

intergen_annot_t0_lt <- mclapply(annot_t0_lt, function(x) {
  if(length(unique(x$gene_id))>1) {
    agr <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
    annot_TxDb <- makeTxDbFromGRanges(agr)
    genic <- genes(annot_TxDb)
    genic <- genic[order(start(ranges(genic)))]
    
    result_fct <- function(x) {
      gene_pair_id <- str_c(genic[x,]$gene_id, collapse = '_')
      gene_pair_strand <- str_c(as.character(strand(genic[x,])), collapse = '_')
      
      genic_reduced <- reduce(genic[x,], ignore.strand=T)
      intergenic <- gaps(genic_reduced)
      if(length(intergenic) == 2) {
        gene_pair_intergenStart <- start(intergenic[2])
        gene_pair_intergenEnd <- end(intergenic[2])
        gene_pair_intergenWidth <- width(intergenic[2])
        
        result_tbl = tibble(idPair = gene_pair_id,
                            strandPair = gene_pair_strand,
                            Start = gene_pair_intergenStart,
                            End = gene_pair_intergenEnd,
                            Width = gene_pair_intergenWidth,
                            overlap=0)
        
        return(result_tbl)
      } else {
        
        tscript_hits <- findOverlaps(genic[x,], ignore.strand=TRUE, drop.self=TRUE)
        tscript_olap <- pintersect(genic[x,][queryHits(tscript_hits)], genic[x,][subjectHits(tscript_hits)], ignore.strand=TRUE)
        
        tscript_olap_width <- width(tscript_olap)[1]
        
        result_tbl = tibble(idPair = gene_pair_id,
                            strandPair = gene_pair_strand,
                            Start = NA,
                            End = NA,
                            Width = 0,
                            overlap=tscript_olap_width)
        return(result_tbl)
      }
    }
    
    result_temp <- rollapply(seq_len(length(genic)), FUN=result_fct, by=1, width=2, align='left')
    result_temp <- as_tibble(result_temp)
    return(result_temp)
  } else {
    return(NA)
  }
}, mc.cores = 4)
intergen_annot_t0_lt <- intergen_annot_t0_lt[!is.na(intergen_annot_t0_lt)]

intergen_annot_t0_tbl <- bind_rows(intergen_annot_t0_lt, .id='sequence')
intergen_annot_t0_tbl$Start <- as.integer(intergen_annot_t0_tbl$Start)
intergen_annot_t0_tbl$End <- as.integer(intergen_annot_t0_tbl$End)
intergen_annot_t0_tbl$Width <- as.integer(intergen_annot_t0_tbl$Width)
intergen_annot_t0_tbl$overlap <- as.integer(intergen_annot_t0_tbl$overlap)
intergen_annot_t0_tbl$overlap[is.na(intergen_annot_t0_tbl$overlap)] <- 0 

table(intergen_annot_t0_tbl$Width == 0)

table(intergen_annot_t0_tbl$overlap == 0)
table(is.na(intergen_annot_t0_tbl$overlap))

table(is.na(intergen_annot_t0_tbl$Start))
filter(intergen_annot_t0_tbl, is.na(overlap))
filter(intergen_annot_t0_tbl, is.na(Start), overlap==0)
filter(intergen_annot_t0_tbl, Width==0, overlap==0)

# stat -------------------------------------------------


idPair_lt <- str_split(intergen_annot_t0_tbl$idPair, pattern = '(?<=T0)_')
idPair_index <- sapply(idPair_lt, function(x) all(x %in% annot_t0_bothUTR_longreadOnly))
intergen_annot_t0_longreadBothUTR_tbl <- intergen_annot_t0_tbl[idPair_index,]


write_tsv(intergen_annot_t0_tbl, 'intergenGenePairsT0_20230329.tsv')
write_tsv(intergen_annot_t0_longreadBothUTR_tbl, 'intergenGenePairsLongreadBothUTR_20230329.tsv')


intergenToCheck <- intergen_annot_t0_longreadBothUTR_tbl

summary(intergenToCheck$Width)
boxplot(intergenToCheck$Width)

summary(filter(intergenToCheck, overlap == 0)$Width)
boxplot(filter(intergenToCheck, overlap == 0)$Width)

# plot --------------------------------------------------

intergenToCheck %>% 
  group_by(strandPair) %>% 
  summarise(median_width = median(Width),
            mean_width = mean(Width),
            sd_width = sd(Width),
            group_n = n()) %>% 
  ungroup() -> intergen_strand_pairs_All_summary

ggplot(intergen_strand_pairs_All_summary, aes(x=strandPair, y=median_width)) +
  geom_col() +
  geom_text(
    aes(label = median_width, y = median_width + 10),
    position = position_dodge(0.9),
    vjust = 0
  )

ggplot(intergen_strand_pairs_All_summary, aes(x=strandPair, y=group_n)) +
  geom_col() +
  geom_text(
    aes(label = group_n, y = group_n + 80),
    position = position_dodge(0.9),
    vjust = 0
  )

intergenToCheck %>% 
  filter(overlap ==0) %>% 
  group_by(strandPair) %>% 
  summarise(median_width = median(Width),
            mean_width = mean(Width),
            sd_width = sd(Width),
            group_n = n()) %>% 
  ungroup() -> intergen_strand_pairs_NotOverlapping_summary

ggplot(intergen_strand_pairs_NotOverlapping_summary, aes(x=strandPair, y=median_width)) +
  geom_col() +
  geom_text(
    aes(label = median_width, y = median_width + 10),
    position = position_dodge(0.9),
    vjust = 0
  )

ggplot(intergen_strand_pairs_NotOverlapping_summary, aes(x=strandPair, y=group_n)) +
  geom_col() +
  geom_text(
    aes(label = group_n, y = group_n + 80),
    position = position_dodge(0.9),
    vjust = 0
  )

intergenToCheck %>% 
  filter(overlap >0) %>% 
  nrow()

p1 <- intergenToCheck %>% 
  filter(overlap >0) %>% 
  count(strandPair) %>% 
  ggplot(aes(x=strandPair, y=n)) +
  geom_col() +
  geom_text(
    aes(label = n, y = n + 80),
    position = position_dodge(0.9),
    vjust = 0
  ) +
  xlab('Overlapping gene pairs') +
  ylab('Number') +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

ggsave(plot = p1, filename = 'overlappingGenePairsDistribution_v01.pdf', device = 'pdf', width = 4, height = 4)

intergenToCheck %>% 
  filter(overlap >0) %>% 
  summarise(median_overlap = median(overlap), 
            mean_overlap = mean(overlap),
            sd_overlap = sd(overlap))

intergenToCheck %>% 
  filter(overlap >0) %>% 
  arrange(desc(overlap))
