# load packages ------------------------------------------------------------------------

suppressPackageStartupMessages({
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

SD <- format(Sys.Date(), '%Y%m%d')

# load data ---------------------------------------------------------------------

## ref. genome -------------------------------------------------
ref_genome <- readDNAStringSet("../ref_genome/fullgenome_jgi/merged_genome/CopciAB_new_jgi_20220113.fasta")

## ref. annot -------------------------------------------------------------
annot_path <- "../remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/clean_annotation_20220427_dir/CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered.gtf"
annot_gr <- import(annot_path)
annot_t0_gr <- annot_gr[str_detect(annot_gr$gene_id, pattern = "T0$")]
annot_t0_df <- as.data.frame(annot_t0_gr)

length(unique(annot_t0_df$gene_id))


annot_gr$gene_id[str_detect(annot_gr$gene_id, pattern = 'T0$')] %>% unique() %>% length()


length(unique(annot_gr$gene_id[str_detect(annot_gr$gene_id, pattern = "T0$", negate = T)]))

length(unique(str_remove(annot_gr$gene_id[str_detect(annot_gr$gene_id, pattern = "T0$", negate = T)], pattern = "\\.T[0-9]+$"))) # 1053

annot_gr$gene_id %>% 
  unique() %>% 
  str_split(pattern = '\\.', simplify = TRUE) %>% 
  as.data.frame() %>% setNames(c('geneId', 'isoformN')) -> isoCount_tbl

isoCount_tbl %>% count(geneId) %>% count(n)

isoCount_tbl %>% 
  group_by(geneId) %>% 
  summarise(maxIsoN = max(isoformN)) %>% 
  ungroup() %>% 
  count(maxIsoN)

iFormSource_tbl <- read_tsv('CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered_T0isoforme_source.tsv')

table(iFormSource_tbl$source_index_short)


annot_t0LongReadOnly_gr <- annot_t0_gr[str_remove(annot_t0_gr$gene_id, pattern = '\\.T0$') %in% filter(iFormSource_tbl, source_index_short == 'longread')$tID_T0]
length(unique(annot_t0LongReadOnly_gr$gene_id))

# motif identification ------------------------------------------------------------------------------
annot_t0LongReadOnly_df <- as.data.frame(annot_t0LongReadOnly_gr)
annot_t0LongReadOnly_TxDb <- makeTxDbFromGRanges(annot_t0LongReadOnly_gr)

# ref_genome

annot_genes_gr <- genes(annot_t0LongReadOnly_TxDb)
annot_exon_gr <- exonsBy(annot_t0LongReadOnly_TxDb, by = "tx", use.names=TRUE)
annot_cds_gr <- cdsBy(annot_t0LongReadOnly_TxDb, by = "tx", use.names=TRUE)
annot_introns_gr <- intronsByTranscript(annot_t0LongReadOnly_TxDb, use.names=TRUE)
annot_5UTR_gr <- fiveUTRsByTranscript(annot_t0LongReadOnly_TxDb, use.names=TRUE)
annot_3UTR_gr <- threeUTRsByTranscript(annot_t0LongReadOnly_TxDb, use.names=TRUE)

## TSS ----------------------------------------------------------------------------------------------------------

length(annot_5UTR_gr)

annot_genes_df <- as.data.frame(annot_genes_gr)
annot_genes_width5UTR_df <- filter(annot_genes_df, gene_id %in% names(annot_5UTR_gr))

annot_genes_width5UTR_lt <- split(annot_genes_width5UTR_df, annot_genes_width5UTR_df$gene_id)
tss_range_lt <- lapply(annot_genes_width5UTR_lt, function(x) {
  if(x$strand == "-") {
    result <- tibble(seqnames = x$seqnames,
                     start = x$end-4,
                     end = x$end +5,
                     strand = x$strand,
                     gene_id = x$gene_id)
  } else {
    result <- tibble(seqnames = x$seqnames,
                     start = x$start-5,
                     end = x$start +4,
                     strand = x$strand,
                     gene_id = x$gene_id)
  }
  return(result)
})
tss_range_tbl <- bind_rows(tss_range_lt)

tss_range_gr <- makeGRangesFromDataFrame(tss_range_tbl, keep.extra.columns = TRUE)
tss_range_seq <- BSgenome::getSeq(ref_genome, tss_range_gr)
names(tss_range_seq) <- tss_range_gr$gene_id

length(tss_range_seq)
table(width(tss_range_seq))

dnaPositionInfo <- as.data.frame(str_split(as.character(tss_range_seq), pattern = "", simplify = T))
dnaPositionInfo <- apply(dnaPositionInfo,2,table)
dnaPositionInfo <- dnaPositionInfo/length(tss_range_seq)

library("seqLogo")
library("ggseqlogo")

pwm <- makePWM(dnaPositionInfo)
summary(pwm)
seqLogo(pwm, ic.scale=TRUE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15)

ggplot() + geom_logo(as.character(tss_range_seq),method = 'bits') + theme_logo()


## TSS 50 ----------------------------------------------------------------------------------------------------------

annot_genes_df <- as.data.frame(annot_genes_gr)
annot_genes_width5UTR_df <- filter(annot_genes_df, gene_id %in% names(annot_5UTR_gr))

annot_genes_width5UTR_lt <- split(annot_genes_width5UTR_df, annot_genes_width5UTR_df$gene_id)
tss50_range_lt <- lapply(annot_genes_width5UTR_lt, function(x) {
  if(x$strand == "-") {
    result <- tibble(seqnames = x$seqnames,
                     start = x$end-49,
                     end = x$end +50,
                     strand = x$strand,
                     gene_id = x$gene_id)
  } else {
    result <- tibble(seqnames = x$seqnames,
                     start = x$start-50,
                     end = x$start +49,
                     strand = x$strand,
                     gene_id = x$gene_id)
  }
  return(result)
})
tss50_range_tbl <- bind_rows(tss50_range_lt)

tss50_range_gr <- makeGRangesFromDataFrame(tss50_range_tbl, keep.extra.columns = TRUE)
tss50_range_seq <- BSgenome::getSeq(ref_genome, tss50_range_gr)
names(tss50_range_seq) <- tss50_range_gr$gene_id

subseq(tss50_range_seq[1], 51, 100)

(table(as.character(subseq(tss50_range_seq, 51, 51))) / length(tss50_range_seq))*100
(table(as.character(subseq(tss50_range_seq, 50, 50))) / length(tss50_range_seq))*100

# tata box
table(as.character(subseq(tss50_range_seq, 17, 23))) %>% 
  as.data.frame() %>% 
  mutate(Freq_p = Freq/length(tss50_range_seq)) %>%
  arrange(desc(Freq)) %>% 
  head()

as.character(subseq(tss50_range_seq, 17, 23)) %>% str_detect(pattern = 'TATATA') %>% table() 
as.character(subseq(tss50_range_seq, 17, 23)) %>% str_detect(pattern = 'TATAAA') %>% table() 

as.character(subseq(tss50_range_seq, 17, 23)) %>% str_detect(pattern = 'TATA') %>% table()
as.character(subseq(tss50_range_seq, 17, 23)) %>% str_detect(pattern = 'TATA.{2}') %>% table()


dnaPositionInfo <- as.data.frame(str_split(as.character(tss50_range_seq), pattern = "", simplify = T))
dnaPositionInfo <- apply(dnaPositionInfo,2,table)
dnaPositionInfo <- dnaPositionInfo/length(tss50_range_seq)

library("seqLogo")
library("ggseqlogo")

pwm <- makePWM(dnaPositionInfo)
summary(pwm)
seqLogo(pwm, ic.scale=TRUE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15)

ggseqlogo(as.character(tss50_range_seq))

ggplot() + 
  geom_logo(as.character(tss50_range_seq),method = 'bits') + 
  theme_logo() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_x_discrete(limits=as.character(c(-50:-1, 1:50)))

ggsave('TSS_230308.pdf', device = 'pdf', width = 12, height = 4)

c(-50:-1, 1:50)[51]

### filter by gene expression --------------------------------------------------

iFormSource_tbl %>% filter(tID_T0 %in% str_remove(names(tss_range_seq), pattern = '\\.T0$')) %>% filter(read_n >= quantile(read_n, p=0.90, na.rm=TRUE)) -> iFormSource_topXp_tbl

tss_rangetopXp_seq <- tss50_range_seq[str_remove(names(tss50_range_seq), pattern = '\\.T0$') %in% iFormSource_topXp_tbl$tID_T0]


dnaPositionInfoXp <- as.data.frame(str_split(as.character(tss_rangetopXp_seq), pattern = "", simplify = T))
dnaPositionInfoXp <- apply(dnaPositionInfoXp,2,table)
dnaPositionInfoXp <- dnaPositionInfoXp/length(tss_rangetopXp_seq)

(table(as.character(subseq(tss_rangetopXp_seq, 51, 51))) / length(tss_rangetopXp_seq))*100
(table(as.character(subseq(tss_rangetopXp_seq, 50, 50))) / length(tss_rangetopXp_seq))*100

# tata box
table(as.character(subseq(tss_rangetopXp_seq, 17, 23))) %>%
  as.data.frame() %>%
  mutate(Freq_p = Freq/length(tss_rangetopXp_seq)) %>%
  arrange(desc(Freq)) %>%
  head()


pwm <- makePWM(dnaPositionInfoXp)
summary(pwm)
seqLogo(pwm, ic.scale=TRUE, xaxis=TRUE, yaxis=TRUE, xfontsize=15, yfontsize=15)

ggseqlogo(as.character(tss_rangetopXp_seq))

ggplot() + 
  geom_logo(as.character(tss_rangetopXp_seq),method = 'bits') + 
  theme_logo() + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_x_discrete(limits=as.character(c(-50:-1, 1:50)))

ggsave('TSS_top10p_230308.pdf', device = 'pdf', width = 12, height = 4)

### more precise identification of the TATA box --------------------------------------------------

inputseq <- subseq(tss_rangetopXp_seq,1,50)

ggplot() + 
  geom_logo(as.character(inputseq),method = 'bits') + 
  theme_logo() + 
  theme(axis.text.x = element_text(angle = 90)) 

table(str_count(as.character(inputseq), pattern = 'TATA.{2}'))
sum(table(str_count(as.character(inputseq), pattern = 'TATA.{2}'))[-1])/length(inputseq)

inputseq_tata_seq <- inputseq[str_detect(as.character(inputseq), pattern = 'TATA.{2}')]

tata_posSeq_lt <- lapply(inputseq_tata_seq, function(x) {
  tata_pos <- str_locate_all(x, pattern = 'TATA.{2}')
  tata_pos <- as_tibble(tata_pos[[1]])
  tata_pos$tata_seq <- apply(tata_pos, 1, function(y) {return(as.character(subseq(x, start = y['start'], end = y['end'])))})
  return(tata_pos)
})

tata_posSeq_tbl <- bind_rows(tata_posSeq_lt, .id='geneID')


barplot(table(tata_posSeq_tbl$start))
tata_posSeq_tbl %>% count(start) %>% arrange(desc(n))
tata_posSeq_tbl %>% count(start) %>% ggplot(aes(x=start, y=n)) + geom_bar(stat = 'identity')

tata_posSeq_tbl %>% filter(tata_seq %in% c('TATATA', 'TATAAA')) %>% summarise(unique_geneID = length(unique(geneID)))
tata_posSeq_tbl %>% filter(tata_seq %in% c('TATATA', 'TATAAA')) %>% count(start) %>% arrange(desc(n))
tata_posSeq_tbl %>% filter(tata_seq %in% c('TATATA', 'TATAAA')) %>%count(start) %>% ggplot(aes(x=start, y=n)) + geom_bar(stat = 'identity')


tata_posSeq_tbl %>% filter(!tata_seq %in% c('TATATA', 'TATAAA')) %>% count(start) %>% arrange(desc(n))
tata_posSeq_tbl %>% filter(!tata_seq %in% c('TATATA', 'TATAAA')) %>% count(start) %>% ggplot(aes(x=start, y=n)) + geom_bar(stat = 'identity')

tata_posSeq_tbl %>% 
  mutate(tata_index = ifelse(tata_seq %in% c('TATATA', 'TATAAA'), 'con', 'NOTcon')) %>% 
  count(tata_index, start) %>% 
  ggplot(aes(x=start, y=n, color=tata_index, group=tata_index)) + 
  geom_point(stat = 'identity') +
  geom_line() +
  theme_bw()

# -38:-31 = 13:20
tata_posSeq_tbl %>% mutate(tata_index = ifelse(start %in% 13:20, TRUE, FALSE)) %>% count(tata_index)
tata_posSeq_tbl %>% mutate(tata_index = ifelse(start %in% 13:20, TRUE, FALSE)) %>% group_by(tata_index) %>% summarise(unique_geneID = length(unique(geneID)))

tata_posSeq_tbl %>% 
  mutate(tata_index_conserv = ifelse(tata_seq %in% c('TATATA', 'TATAAA'), TRUE, FALSE),
         tata_index_position = ifelse(start %in% 13:20, TRUE, FALSE)) -> tata_posSeq_plus_tbl

tata_posSeq_plus_tbl %>% 
  filter(tata_index_conserv) %>% 
  group_by(tata_index_position) %>% 
  summarise(pos_sum = n()) %>% 
  ungroup() %>% 
  mutate(pos_p = pos_sum / sum(pos_sum))

tata_posSeq_plus_tbl %>% 
  filter(tata_index_conserv) %>% 
  group_by(tata_index_position) %>%
  summarise(uni_geneID=length(unique(geneID)))
             
length(unique(filter(tata_posSeq_plus_tbl, tata_index_conserv, tata_index_position)$geneID))

tata_posSeq_tbl %>% count(tata_seq) %>% arrange(desc(n)) # az összes
tata_posSeq_tbl %>% filter(tata_seq %in% c('TATATA', 'TATAAA')) %>% count(tata_seq) %>% arrange(desc(n)) 
tata_posSeq_tbl %>% filter(! tata_seq %in% c('TATATA', 'TATAAA')) %>% count(tata_seq) %>% arrange(desc(n)) 


## TIS / START -----------------------------------------------------------------------------------------------
# -12/+12

as.data.frame(sum(width(annot_5UTR_gr))) %>% 
  setNames('UTR5Width') %>% 
  rownames_to_column(var='seqnames') %>% 
  filter(UTR5Width > 12) %>% 
  mutate(start = UTR5Width-12+1,
         end = UTR5Width + 12,
         gene_id = seqnames) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) -> TISregion_gr

exon_seq <- extractTranscriptSeqs(ref_genome, annot_exon_gr)

setdiff(TISregion_gr$gene_id, names(exon_seq))

TISregion_seq <- BSgenome::getSeq(exon_seq, TISregion_gr)
names(TISregion_seq) <- TISregion_gr$gene_id

ggplot() + 
  geom_logo(as.character(TISregion_seq),method = 'bits') + 
  theme_logo() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(limits=as.character(c(-12:-1, 1:12)))
  
ggsave('TIS_motif_230310.pdf', device = 'pdf', width = 6, height = 5)

ggplot() + geom_logo(as.character(TISregion_seq),method = 'bits') + theme_logo()


TISregion_rangetopXp_seq <- TISregion_seq[str_remove(names(TISregion_seq), pattern = '\\.T0$') %in% iFormSource_topXp_tbl$tID_T0]

ggplot() + 
  geom_logo(as.character(TISregion_rangetopXp_seq),method = 'bits') + 
  theme_logo() + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_x_discrete(limits=as.character(c(-12:-1, 1:12)))

ggsave('TIS_motif_top10p_230310.pdf', device = 'pdf', width = 6, height = 5)

## TT / STOP -----------------------------------------------------------------------------------------------

as.data.frame(sum(width(annot_3UTR_gr))) %>% 
  setNames('UTR3Width') %>% 
  rownames_to_column(var='seqnames') -> UTR3Width_tbl

as.data.frame(sum(width(annot_exon_gr))) %>% 
  setNames('exonWidth') %>% 
  rownames_to_column(var='seqnames') -> exonWidth_tbl

exonUTR3Width_tbl <- left_join(exonWidth_tbl, UTR3Width_tbl, by='seqnames')
exonUTR3Width_tbl %>% 
  mutate(UTR3Width = ifelse(is.na(UTR3Width), 0, UTR3Width)) %>% 
  mutate(start=exonWidth-UTR3Width-2,
         end=exonWidth-UTR3Width,
         gene_id = seqnames) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) -> TTregion_gr

exon_seq <- extractTranscriptSeqs(ref_genome, annot_exon_gr)

setdiff(TTregion_gr$gene_id, names(exon_seq))

TTregion_seq <- BSgenome::getSeq(exon_seq, TTregion_gr)
names(TTregion_seq) <- TTregion_gr$gene_id

table(TTregion_seq)
length(TTregion_seq)
round((table(TTregion_seq)/length(TTregion_seq))*100, 1)


ggplot() + geom_logo(as.character(TTregion_seq),method = 'bits') + theme_logo()

ggsave('TT_motif_230314.pdf', device = 'pdf', width = 4, height = 4)

## Splice site ----------------------------------------------------------------------

annot_exon_length <- sapply(width(annot_exon_gr), length)
annot_multyexon_gr <- annot_exon_gr[annot_exon_length > 1]

length(annot_multyexon_gr)

### SS1 -----------------------------------------------------------------------------

splices5_range_lt <- lapply(annot_multyexon_gr, function(x) {
  if(any(strand(x) == "+")) {
    result <- tibble(start = end(x)[-length(x)]+1,
                     end = end(x)[-length(x)]+6,
                     strand = as.character(unique(strand(x))),
                     seqnames = as.character(unique(seqnames(x))))
  } else {
    result <- tibble(start = start(x)[-length(x)]-6,
                     end = start(x)[-length(x)]-1,
                     strand = as.character(unique(strand(x))),
                     seqnames = as.character(unique(seqnames(x))))
  }
  return(result)
})

splices5_range_tbl <- bind_rows(splices5_range_lt, .id = "gene_id")

splices5_range_gr <- makeGRangesFromDataFrame(splices5_range_tbl, keep.extra.columns = TRUE)
splices5_range_seq <- BSgenome::getSeq(ref_genome, splices5_range_gr)
names(splices5_range_seq) <- splices5_range_gr$gene_id

length(splices5_range_seq) 
length(unique(names(splices5_range_seq)))

p_splice5 <- ggplot() + 
  geom_logo(as.character(splices5_range_seq),method = 'bits') + 
  theme_logo()

ggsave(plot = p_splice5, filename = str_c('splice5_motif_',SD,'.pdf'), device = 'pdf', width = 4, height = 4)

sj5_mot <- sapply(str_split(as.character(splices5_range_seq), pattern = ''), function(x) str_c(x[1:2], collapse = '')) 
table(sj5_mot)
splices5_range_seq[!sj5_mot %in% c('GT', 'GC', 'AT')]

### SS2 -----------------------------------------------------------------------------

splices3_range_lt <- lapply(annot_multyexon_gr, function(x) {
    if(any(strand(x) == "+")) {
    result <- tibble(start = start(x)[-1]-3,
                     end = start(x)[-1]-1,
                     strand = as.character(unique(strand(x))),
                     seqnames = as.character(unique(seqnames(x))))
  } else {
    result <- tibble(start = end(x)[-1]+1,
                     end = end(x)[-1]+3,
                     strand = as.character(unique(strand(x))),
                     seqnames = as.character(unique(seqnames(x))))
  }
  return(result)
})

splices3_range_tbl <- bind_rows(splices3_range_lt, .id = "gene_id")

splices3_range_gr <- makeGRangesFromDataFrame(splices3_range_tbl, keep.extra.columns = TRUE)
splices3_range_seq <- BSgenome::getSeq(ref_genome, splices3_range_gr)
names(splices3_range_seq) <- splices3_range_gr$gene_id

length(splices3_range_seq) 

p_splice3 <- ggplot() + 
  geom_logo(as.character(splices3_range_seq),method = 'bits') + 
  theme_logo()

ggsave(plot = p_splice3, filename = str_c('splice3_motif_',SD,'.pdf'), device = 'pdf', width = 4, height = 4)

sj3_mot <- sapply(str_split(as.character(splices3_range_seq), pattern = ''), function(x) str_c(x[2:3], collapse = '')) 
table(sj3_mot)
splices3_range_seq[!sj3_mot %in% c('AC', 'AG')]

sj53_mot <- str_c(sj5_mot, sj3_mot, sep = '_')
table(sj53_mot)
splices3_range_gr[sj53_mot %in% c("AT_AC")]

table(sj53_mot %in% c("GT_AG", "GC_AG","AT_AC")) / length(sj53_mot)
sj53_mot_conserve <- sj53_mot[sj53_mot %in% c("GT_AG", "GC_AG","AT_AC")]
round((table(sj53_mot_conserve) / length(sj53_mot_conserve))*100,2)

## TTS / PAS -------------------------------------------------------------------------------------

annot_genes_lt <- split(annot_genes_df, annot_genes_df$gene_id)
annot_genes_lt <- annot_genes_lt[names(annot_genes_lt) %in%  names(annot_3UTR_gr)] 

TTS_lt <- lapply(annot_genes_lt, function(x) {
  if(x$strand == '-') {
    x$width <- NULL
    x$end <- x$start + 100
    x$start <- x$start - 100
  } else {
    x$width <- NULL
    x$start <- x$end - 100
    x$end <- x$end + 100
  }
  return(x)
})

TTS_tbl <- bind_rows(TTS_lt) 
TTS_gr <- makeGRangesFromDataFrame(TTS_tbl, keep.extra.columns = TRUE)

TTS_seq <- BSgenome::getSeq(ref_genome, TTS_gr)
names(TTS_seq) <- TTS_gr$gene_id

TTS_nucMatrix <- as.data.frame(str_split(as.character(TTS_seq), pattern = '', simplify = TRUE))
TTS_frecnucMatrix <- apply(TTS_nucMatrix, 2, table)
TTS_frecnucMatrix_norm <- TTS_frecnucMatrix/length(TTS_seq)

TTS_frecnucMatrix_norm %>% 
  as.data.frame() %>% 
  rownames_to_column(var='nuc') %>% 
  gather(key='position', value = 'nuc_p', -nuc) %>% 
  mutate(position = as.integer(str_remove(position, pattern = '^V'))) -> TTS_frecnucMatrix_norm_short

TTS_frecnucMatrix_norm_short %>% 
  mutate(position = position-101) %>%
  filter(position <= 50 & position >= -50) %>% 
  ggplot(aes(x=position, y = nuc_p, group=nuc, color=nuc)) +
  geom_line() +
  scale_color_manual(breaks = c('A', 'C', 'G', 'T', 'U'), 
                    values=c('#109648', '#255C99', '#F7B32B', '#D62839', '#D62839')) +
  geom_vline(aes(xintercept = 0), linetype=2, size=0.5) +
  theme_bw()

ggsave('TTS_PAS_nucfreq_230314.pdf', device = 'pdf', width = 14, height = 6)

TTS_frecnucMatrix_norm_center <- as.data.frame(TTS_frecnucMatrix_norm)[,c('V99','V100','V101', 'V102')]

TTS_center_seq <- Biostrings::subseq(TTS_seq, start = 51, end = 151)

ggplot() + 
  geom_logo(as.character(TTS_center_seq),method = 'bits') + 
  theme_logo() +
  scale_x_discrete(limits=as.character(-50:50)) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_vline(aes(xintercept = 0), linetype=2, size=0.5)

ggsave('TTS_PAS_nucfreq_logo_230317.pdf', device = 'pdf', width = 14, height = 6)

writeXStringSet(TTS_center_seq, 'TTS_seq.fasta')

### 6 mer seq --------------------------------
library(zoo)

TTS_center_6mers_lt <- lapply(TTS_center_seq, function(y) {
  seqkmers_m <- rollapply(unlist(str_split(as.character(y), pattern = '')), width=6, FUN=function(x) return(x), by=1, align='right')
  seqkmers_tbl <- tibble(seq = apply(seqkmers_m, 1, str_c, collapse=''), start_pos = seq_len(nrow(seqkmers_m)))
  return(seqkmers_tbl)
})

TTS_center_6mers_tbl <- bind_rows(TTS_center_6mers_lt, .id = 'geneID')
count(TTS_center_6mers_tbl, seq) %>% arrange(desc(n))

filter(TTS_center_6mers_tbl, start_pos %in% 25:35) %>% 
  count(seq) %>% 
  arrange(desc(n)) %>% 
  View()

TTS_center_6mers_tbl %>% 
  filter(start_pos == 1) %>% 
  count(seq) %>% 
  filter(seq == 'TTTCTT')

### 6 mer seq plus --------------------------------

TTS_center_6mers_tbl %>% count(seq, name='sum') %>% arrange(desc(sum)) %>% mutate(order = seq_len(length(sum))) -> TTS_center_6mers_seqCount_tbl
TTS_center_6mers_tbl %>% count(start_pos, seq, name='pos_count') -> TTS_center_6mers_posCount_tbl

TTS_center_6mers_posseqCount_tbl <- left_join(TTS_center_6mers_posCount_tbl, TTS_center_6mers_seqCount_tbl, by='seq')


TTS_center_6mers_posseqCount_tbl %>% 
  filter(order <=100) %>% 
  mutate(start_pos = start_pos-50) %>% 
  ggplot(aes(x=start_pos, y=pos_count, color=seq, group=seq)) +
  geom_line() +
  geom_vline(xintercept = c(-28, -20, -14, 1), linetype="dotted")

TTS_center_6mers_posseqCount_tbl %>% 
  filter(order <=100) %>% 
  mutate(start_pos = start_pos-50) %>% 
  filter(start_pos %in% -28:-14) %>% 
  arrange(desc(pos_count)) %>% 
  group_by(seq) %>% 
  mutate(sum_Arich_loc = sum(pos_count)) %>% 
  arrange(desc(sum_Arich_loc), seq) %>%
  ungroup() %>% 
  mutate(order_2=as.integer(factor(seq, levels = unique(seq), ordered = TRUE))) -> top100motifInRrichRegion_tbl

top100motifInRrichRegion_tbl %>% 
  filter(order_2 %in% 1:10) %>%
  ggplot(aes(x=start_pos, y=pos_count, color=seq, group=seq)) +
  geom_line() +
  geom_vline(xintercept = c(-28, -20,-14), linetype="dotted")

length(unique(top100motifInRrichRegion_tbl$seq))

top100motifInRrichRegion_tbl %>% 
  group_by(seq) %>% 
  summarise(sum_in_Arich_region=sum(pos_count)) %>%
  arrange(desc(sum_in_Arich_region)) %>% 
  View()


