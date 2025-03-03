# load libraries ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(ORFik)
  library(Biostrings)
  library(GenomicFeatures)
  library(rtracklayer)
  library(stringr)
  library(readr)
  library(BSgenome)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

# load data ----------------------------------------------------------------------------------
# load annotation

annot_gr <- import("../remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/clean_annotation_20220427_dir/CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered.gtf")
annotT0_gr <- annot_gr[str_detect(annot_gr$transcript_id, pattern = "\\.T0$")]

## load isoform table

iFormSource_path <- '/data/bhegedus/Projects/pj8_3genseq/Transcript_sequencing/Coprinopsis_cinerea/results/refannotation_summary/CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered_T0isoforme_source.tsv'
iFormSource_tbl <- read_tsv(iFormSource_path)

table(iFormSource_tbl$source_index_short)

annot_t0LongReadOnly_gr <- annotT0_gr[str_remove(annotT0_gr$gene_id, pattern = '\\.T0$') %in% filter(iFormSource_tbl, source_index_short == 'longread')$tID_T0]
length(unique(annot_t0LongReadOnly_gr$gene_id))
annot_t0LongReadOnly_gr <- annot_t0LongReadOnly_gr[annot_t0LongReadOnly_gr$gene_id %in% unique(annot_t0LongReadOnly_gr$gene_id[annot_t0LongReadOnly_gr$type == '5UTR'])]
length(unique(annot_t0LongReadOnly_gr$gene_id)) 

annot_TxDb <- makeTxDbFromGRanges(annot_t0LongReadOnly_gr)

# load genome 
genome_seq <- readDNAStringSet("../ref_genome/fullgenome_jgi/merged_genome/CopciAB_new_jgi_20220113.fasta")

# load functional annotation
eggnog_annot_path <- "/data/bhegedus/Projects/pj8_3genseq/Transcript_sequencing/Coprinopsis_cinerea/results/refannotation_summary/eggnog_dir_2210/out.emapper.annotations_cleanHeader"
ipro_annot_path <- "/data/bhegedus/Projects/pj8_3genseq/Transcript_sequencing/Coprinopsis_cinerea/results/remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/clean_annotation_20220427_dir/interproscan88_out/CopciAB_new_annot_CDS_20220425m1_iproscan88_headered_simple.txt"

eggnog_annot <- read_tsv(eggnog_annot_path)
ipro_annot <- read_tsv(ipro_annot_path)


## get 5' UTR ----------------------------------- 

UTR5_grl <- fiveUTRsByTranscript(annot_TxDb, use.names=TRUE)
UTR5_seq <- GenomicFeatures::extractTranscriptSeqs(genome_seq, UTR5_grl)

length(UTR5_seq)
hist(width(UTR5_seq))
boxplot(width(UTR5_seq), ylab="5'UTR length (nt)")
summary(width(UTR5_seq))

fiveUTRs_gr <- loadRegion(annot_TxDb, "leaders")
cds_gr <- loadRegion(annot_TxDb, "cds")

uORF_grList <- findUORFs(fiveUTRs_gr, 
                         "../ref_genome/fullgenome_jgi/merged_genome/CopciAB_new_jgi_20220113.fasta", 
                         startCodon ="ATG",
                         stopCodon = stopDefinition(1),
                         minimumLength = 1, 
                         cds = cds_gr)

uORF_gr <- Reduce(c,uORF_grList)

uORF_gr$gene_id <- names(uORF_gr)
uORF_gr$transcript_id <- uORF_gr$names
names(uORF_gr) <- NULL
uORF_gr$type <- "exon"

export(uORF_gr, "CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered_T0UTR5_uATGuORF.gtf")
uORF_df <- as.data.frame(uORF_gr)

length(unique(annot_t0LongReadOnly_gr$gene_id))

length(unique(str_remove(names(uORF_grList), pattern = "_[0-9]+$")))
length(uORF_grList)

length(unique(uORF_df$gene_id))
length(unique(uORF_df$transcript_id))

uORF_df %>% 
  group_by(gene_id) %>% 
  summarise(uORF_n = length(unique(transcript_id))) %>% 
  ungroup %>% 
  arrange(desc(uORF_n)) -> uORF_n_summary

uORF_df %>% 
  group_by(transcript_id) %>% 
  summarise(sum_width = sum(width)) %>% 
  ungroup() %>% 
  arrange(desc(sum_width)) %>% 
  mutate(aa_width = (sum_width/3)-1) -> uORF_width_summary

head(uORF_width_summary)
tail(uORF_width_summary)

min(uORF_width_summary$sum_width) 
summary(uORF_width_summary$sum_width)
boxplot(uORF_width_summary$sum_width)

summary(uORF_width_summary$aa_width) 
boxplot(uORF_width_summary$aa_width)

uORF_width_summary %>% 
  mutate(transcript_id_source = str_remove(transcript_id, pattern = "_[0-9]+$")) %>% 
  group_by(transcript_id_source) %>% 
  filter(sum_width == max(sum_width)) %>% 
  filter(seq_len(n()) == 1) -> uORF_widthGroupMax_summary

hist(uORF_widthGroupMax_summary$sum_width)
sum(uORF_widthGroupMax_summary$sum_width >= 500)
filter(uORF_widthGroupMax_summary, sum_width >= 500)


length(unique(uORF_n_summary$gene_id))

sum(!uORF_n_summary$gene_id %in% ipro_annot$Protein_accession) 
sum(uORF_n_summary$gene_id %in% ipro_annot$Protein_accession) 

uORFs_ipro_annots <- ipro_annot[ipro_annot$Protein_accession %in% uORF_n_summary$gene_id,]
length(unique(uORFs_ipro_annots$Protein_accession))
filter(uORFs_ipro_annots, Protein_accession=="CopciAB_446268.T0") # arg-2


COG_tbl <- dplyr::select(eggnog_annot, query, COG_category)
COG_category_sep <- as.data.frame(str_split(COG_tbl$COG_category, pattern = "", simplify = T))
COG_long_tbl <- bind_cols(COG_tbl, COG_category_sep)
COG_long_tbl %>% 
  dplyr::select(-COG_category) %>% 
  gather(key=COG_category_multy, value = COG_category, -query) -> COG_short_tbl

COG_short_tbl %>% 
  filter(str_count(COG_category)>0,
         COG_category != "-") -> COG_short_clean_tbl

table(COG_short_clean_tbl$COG_category_multy) 

COG_short_clean_tbl %>% 
  count(COG_category) %>% 
  arrange(desc(n)) %>% 
  mutate(COG_category, COG_category = factor(COG_category, levels = COG_category, ordered = T)) %>% 
  ggplot(aes(x=COG_category, y= n)) +
  geom_bar(stat = "identity", position = "dodge")


sum(!uORF_n_summary$gene_id %in% COG_short_clean_tbl$query) 
sum(uORF_n_summary$gene_id %in% COG_short_clean_tbl$query) 

uORFs_eggnog_annots <- COG_short_clean_tbl[COG_short_clean_tbl$query %in% uORF_n_summary$gene_id,]
filter(uORFs_eggnog_annots, query=="CopciAB_446268.T0") # arg-2

sort(table(uORFs_eggnog_annots$COG_category), decreasing = TRUE)
barplot(sort(table(uORFs_eggnog_annots$COG_category), decreasing = TRUE))

# uORF summary ----------------------------------------------------------------------------------------------------

uORF_width_summary
uORF_n_summary

uORF_n_m <- uORF_n_summary$uORF_n
barplot(table(uORF_n_m), xlab = "uORF count", ylab = "transcript count")
uORF_nGroups <- cut(uORF_n_m,breaks = c(0,1,2,3,4,5,6,7,8,9,10,max(uORF_n_m)), include.lowest = T)
barplot(table(uORF_nGroups)) 

(table(uORF_nGroups)/length(uORF_nGroups))*100
table(uORF_n_m[uORF_nGroups == '[0,1]'])
(sum(uORF_n_m == 1)/length(uORF_n_m))*100 

all_uORF_length  <- uORF_width_summary$sum_width
summary(all_uORF_length)
boxplot(all_uORF_length)


uORF_lt <- split(uORF_df, uORF_df$gene_id)
uORF_diffFormTSS_lt <- lapply(uORF_lt, function(x) {
  inner_list <- split(x, x$transcript_id)
  inner_result_lt <- lapply(inner_list, function(y) {
    if(any(y$strand == "-")) {
      result <- tibble(transcript_id = unique(y$transcript_id),
                       diffFromTSS = sum(unlist(apply(as.data.frame(ranges(fiveUTRs_gr[[unique(y$gene_id)]])),1,function(x) seq(x["start"], x["end"])))>max(y$end)))
    } else {
      result <- tibble(transcript_id = unique(y$transcript_id),
                       diffFromTSS = sum(unlist(apply(as.data.frame(ranges(fiveUTRs_gr[[unique(y$gene_id)]])),1,function(x) seq(x["start"], x["end"])))<min(y$start)))
    }
    return(result)
  })
  inner_result_tbl <- bind_rows(inner_result_lt)
  return(inner_result_tbl)
})

uORF_diffFormTSS_tbl <- bind_rows(uORF_diffFormTSS_lt) %>% arrange(diffFromTSS)
hist(uORF_diffFormTSS_tbl$diffFromTSS)
boxplot(uORF_diffFormTSS_tbl$diffFromTSS)
summary(uORF_diffFormTSS_tbl$diffFromTSS)
sum(uORF_diffFormTSS_tbl$diffFromTSS<=20)
table(uORF_diffFormTSS_tbl$diffFromTSS<=20)

table(sapply(uORF_diffFormTSS_lt, function(x) {all(x$diffFromTSS <=20)}))

as.data.frame(width(UTR5_grl)) %>% 
  dplyr::select(-group) %>% 
  rename("gene_id" = "group_name", "width" = "value") %>% 
  group_by(gene_id) %>% 
  summarise(sum_width = sum(width)) -> UTR5_width_tbl


uORF_diffFormTSS_tbl %>% 
  mutate(gene_id = str_remove(transcript_id, pattern = "_[0-9]+$")) %>% 
  left_join(UTR5_width_tbl, by="gene_id") -> uORF_diffFormTSS_plus_tbl

uORF_diffFormTSS_plus_tbl %>% 
  mutate(diffFromTSS_p = diffFromTSS/sum_width) %>% 
  ggplot(aes(x=diffFromTSS_p)) +
  geom_density() +
  xlim(c(0,1))


utrAnnot_gr <- import('CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered_T0UTR5_uATGuORF.gtf')
utrAnnot_exon_df <- as.data.frame(utrAnnot_gr)
utrAnnot_cds_df <- mutate(utrAnnot_exon_df, type = 'CDS')

utrAnnot_df <- bind_rows(utrAnnot_exon_df, utrAnnot_cds_df)
utrAnnot_df <- arrange(utrAnnot_df, nchar(as.character(seqnames)), as.character(seqnames), start)

export(utrAnnot_df, 'CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered_T0UTR5_uATGuORFcds.gtf')