library(readr)
library(stringr)
library(dplyr)
library(rtracklayer)
library(Biostrings)

############################# BLAST ----

# Run blast
# makeblastdb -in CopciAB_JGI_ONTXie.fasta -input_type fasta -dbtype prot -parse_seqids -out CopciAB_JGI_ONTXie
# blastp -db CopciAB_JGI_ONTXie -query ../../../CopciAB_new_annot_CDS_20220127.prot.fasta -evalue 0.05 -outfmt "6 qseqid sseqid qlen qstart qend sstart send evalue bitscore pident nident ppos positive length gaps" -out CopciAB_20220127_JGI_ONTXie.txt -num_threads 20

sj_path <- "Krizsi_JGI_Novogene_SJ.out.tab"
sj_names <- c("seqnames", "start", "end", "strand_index", "intron_motif_index", "annot_index", "uniquely_mapped_reads", "multimapping_reads", "maximum_spliced_alignment_overhang")
sj_tbl <- read_tsv(sj_path, col_names = sj_names)
sj_tbl <- mutate(sj_tbl, index = str_c(seqnames, start, end, sep = "-"))

sj_tbl$strand_index <- c(NA, "+", "-")[sj_tbl$strand_index+1]
sj_tbl$intron_motif_index <- c(NA, "GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AT")[sj_tbl$intron_motif_index+1]

annot_sj <- read_tsv("../../../CopciAB_new_annot_CDS_20220127_splice_info.tsv")
annot_sj <- mutate(annot_sj, index = str_c(seqnames, start, end, sep = "-"))

annot_sj_c <- left_join(annot_sj, sj_tbl[c(-1,-2,-3,-6,-9)], by="index")
write_tsv(annot_sj_c, "CopciAB_20220127_SJ_coverage_info.tsv")

annot_sj_c %>% 
  group_by(transcript_id) %>% 
  summarise(splice_c_unique = str_c(unique(splice_c), collapse = "_"),
            splice_c_unique_index = all(splice_c %in% c("GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AG")),
            uniquely_mapped_reads_mean = mean(uniquely_mapped_reads, na.rm=TRUE),
            uniquely_mapped_reads_sd = sd(uniquely_mapped_reads, na.rm=TRUE),
            uniquely_mapped_reads_min = min(uniquely_mapped_reads, na.rm=TRUE)) -> annot_sj_c_summary

write_tsv(annot_sj_c_summary, "CopciAB_20220127_SJ_coverage_info_summary.tsv")

annot_gr <- import("../../../CopciAB_new_annot_CDS_20220127.gtf")

annot_gr %>% 
  as.data.frame() %>% 
  filter(type == "exon") %>% 
  group_by(gene_id) %>% 
  summarise(exon_n=n()) %>% 
  ungroup() -> exon_n_tbl

annot_dict_tbl <- tibble(gene_id = unique(annot_gr$gene_id))
cds_index <- unique(annot_gr$gene_id[annot_gr$type == "CDS"])
annot_dict_tbl$protein_id <- cds_index[match(annot_dict_tbl$gene_id, cds_index)]

annot_dict_tbl <- left_join(annot_dict_tbl, exon_n_tbl, by="gene_id")

perfect_splicesite_n <- read_tsv("../../../splice_perfectmatch_site_match_reads_20220127_summary.tsv")
annot_dict_tbl <- left_join(annot_dict_tbl, perfect_splicesite_n, by=c("gene_id" = "transcript_id"))

sum(is.na(annot_dict_tbl$read_n[str_detect(annot_dict_tbl$gene_id, pattern = "^G")]))
sum(str_detect(annot_dict_tbl$gene_id, pattern = "^G") & annot_dict_tbl$exon_n == 1)
sum(str_detect(annot_dict_tbl$gene_id, pattern = "^G") & annot_dict_tbl$exon_n != 1 & is.na(annot_dict_tbl$read_n))
annot_dict_tbl[str_detect(annot_dict_tbl$gene_id, pattern = "^G") & annot_dict_tbl$exon_n != 1 & is.na(annot_dict_tbl$read_n),] # van olyan ONT annotáció amely nem kap splice lefedettséget mert rosszul sikerült az illesztés

sum(is.na(annot_dict_tbl$read_n) & annot_dict_tbl$exon_n > 1)
annot_dict_tbl$gene_id[is.na(annot_dict_tbl$read_n) & annot_dict_tbl$exon_n > 1]
table(sapply(str_split(annot_dict_tbl$gene_id[is.na(annot_dict_tbl$read_n) & annot_dict_tbl$exon_n > 1], pattern = "\\.|_|[0-9]+"), "[[", 1))


ipro_header <- c("Protein_accession", "MD5_digest", "Sequence_length", "Analysis", "Signature_accession", "Signature_description", "Start_location", "Stop_location", "Score_evalue", "Status", "Date", "InterPro_annotations_accession", "InterPro_annotations_description", "GO_annotations", "Pathways_annotations")
write_tsv(as.data.frame(matrix(ipro_header, nrow = 1)), file = "../../../interproscan_dir/ipro_header.tsv" ,col_names = F)
ipro_path <- "../../../interproscan_dir/CopciAB_new_annot_CDS_20220127_iproscan82_header.txt"
ipro_tbl <- read_tsv(ipro_path)

ipro_tbl %>% 
  filter(!Analysis %in% c("Coils", "MobiDBLite")) -> ipro_f1_tbl
  
ipro_f1_lt <- split(ipro_f1_tbl, ipro_f1_tbl$Protein_accession)
ipro_f1_compact_lt <- lapply(ipro_f1_lt, function(x) {
  u_ipro <- x[!duplicated(x$Signature_accession),]
  c_ipro <- tibble(Analisys = str_c(x$Analysis, collapse = "__"),
                   Signature_accession = str_c(x$Signature_accession, collapse = "__"),
                   InterPro_annotations_accession = str_c(x$InterPro_annotations_accession, collapse = "__"),
                   InterPro_annotations_description = str_c(x$InterPro_annotations_description, collapse = "__"))
  return(c_ipro)
})

ipro_f1_compact_tbl <- bind_rows(ipro_f1_compact_lt, .id = "protein_id")

annot_dict_tbl <- left_join(annot_dict_tbl, ipro_f1_compact_tbl, by=c("gene_id" = "protein_id"))

filter(annot_dict_tbl, str_detect(gene_id, pattern =  "G4236"))

annot_dict_tbl <- left_join(annot_dict_tbl, annot_sj_c_summary, by=c("gene_id"="transcript_id"))
annot_dict_tbl <- select(annot_dict_tbl, 1:4, 9:13, everything())

filter(annot_dict_tbl, str_detect(gene_id, pattern = "^CC2G"))


prot_seq <- readAAStringSet("../../../CopciAB_new_annot_CDS_20220127.prot.fasta")
length(prot_seq)

prote_width <- tibble(gene_id = names(prot_seq),
                      width_aa = width(prot_seq))

length(unique(annot_gr$gene_id))
length(unique(annot_gr$gene_id[annot_gr$type == "CDS"]))

blast_cnames <- c("qseqid sseqid qlen qstart qend sstart send evalue bitscore pident nident ppos positive length gaps")
blast_cnames <- unlist(str_split(blast_cnames, pattern = " "))
blast_tbl <- read_tsv("CopciAB_20220127_JGI_ONTXie.txt", col_names=blast_cnames)

blast_tbl <- mutate(blast_tbl, strain_index = ifelse(str_detect(sseqid, pattern = "^CopciAB_"), "JGI", "ONT")) 
blast_tbl <- mutate(blast_tbl, qir = nident/qlen)

setdiff(names(prot_seq),blast_tbl$qseqid)
table(sapply(str_split(setdiff(names(prot_seq),blast_tbl$qseqid), pattern = "_|\\.|[0-9]{2,}"), "[[", 1))


group_by(blast_tbl, strain_index) %>% summarise(length(unique(sseqid)))

blast_strain_lt <- split(blast_tbl, blast_tbl$strain_index) 


blast_lt <- split(blast_tbl, blast_tbl$qseqid)
sum(sapply(blast_lt, function(x) return(any(x$pident == 100) & any(x$strain_index == "JGI")))) 
sum(sapply(blast_lt, function(x) return(any(x$pident == 100) & any(x$strain_index == "ONT")))) 

sum(sapply(blast_lt, function(x) return(all(x$pident != 100) & any(x$strain_index == "JGI")))) 
sum(sapply(blast_lt, function(x) return(all(x$pident != 100) & any(x$strain_index == "ONT")))) 

# JGI ---------------------------------------------

JGI_blast_tbl <- blast_strain_lt$JGI
JGI_blast_lt <- split(JGI_blast_tbl, JGI_blast_tbl$qseqid)

length(setdiff(names(prot_seq),JGI_blast_tbl$qseqid))
table(sapply(str_split(setdiff(names(prot_seq),JGI_blast_tbl$qseqid), pattern = "_|\\.|[0-9]+\\."), "[[", 1))


sum(sapply(JGI_blast_lt, function(x) any(x$qir == 1)))
sum(sapply(JGI_blast_lt, function(x) any(x$qir >= 0.9)))

JGI_blast_max_lt <- lapply(JGI_blast_lt, function(x) return(x[x$qir == max(x$qir),]))

table(sapply(JGI_blast_max_lt, nrow))
JGI_blast_max_lt[sapply(JGI_blast_max_lt, nrow) == 3]


for(i in seq_len(length(JGI_blast_max_lt))) {
  JGI_blast_max_lt[[i]]$index <- i
}

JGI_blast_max_tbl <- bind_rows(JGI_blast_max_lt)

View(filter(JGI_blast_max_tbl, str_detect(qseqid, pattern="^G")))

hist(JGI_blast_max_tbl$qir)


annot_dict_tbl

JGI_blast_max_tbl %>% 
  filter(qir == 1) %>% 
  summarise(qseqid_dupli_n = sum(duplicated(qseqid)),
            sseqid_dupli_n = sum(duplicated(sseqid)))

sum(sapply(JGI_blast_max_lt, function(x) sum(x$qir == 1)) > 1)
JGI_blast_max_lt[sapply(JGI_blast_max_lt, function(x) sum(x$qir == 1)) > 1]

JGI_blast_max_tbl %>% 
  filter(qir >= 0.9) %>% 
  summarise(qseqid_dupli_n = sum(duplicated(qseqid)),
            sseqid_dupli_n = sum(duplicated(sseqid)))

arrange(JGI_blast_max_tbl, qir)

table(sapply(JGI_blast_max_lt, nrow))
JGI_blast_max_lt[sapply(JGI_blast_max_lt, nrow) > 1]

JGI_blast_max_compact_lt <- lapply(JGI_blast_max_lt, function(x) {
  result <- tibble(qseqid = str_c(x$sseqid, collapse = " | "),
                   pident = str_c(x$pident, collapse = " | "),
                   qir = str_c(x$qir, collapse = " | "),
                   target_n = nrow(x))
  return(result)
})

JGI_blast_max_compact_tb <- bind_rows(JGI_blast_max_compact_lt, .id = "transcript_id")

annot_dict2_tbl <- left_join(annot_dict_tbl, JGI_blast_max_compact_tb, by=c("gene_id"="transcript_id"))

filter(annot_dict2_tbl, str_detect(gene_id, pattern = "^CC2G"))


annot_dict2_tbl <- left_join(annot_dict2_tbl, prote_width, by="gene_id")
annot_dict2_tbl <- select(annot_dict2_tbl, gene_id, protein_id, width_aa, everything())


JGI_blast_max_target_lt <- split(JGI_blast_max_tbl, JGI_blast_max_tbl$sseqid)
JGI_blast_max_target_lt[sapply(JGI_blast_max_target_lt <- split(JGI_blast_max_tbl, JGI_blast_max_tbl$sseqid), nrow) > 1]

count_rpk_path <- "/node9_data/bhegedus/pj8_3genseq/Transcript_sequencing/remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/illumina_splice_site/FCounts_dir/CopciAB_2022-02-07_rpk_summary.tsv"
count_tbl <- read_tsv(count_rpk_path)

annot_dict2_tbl <- read_tsv("CopciAB_ref_20210415_dict.tsv")
annot_dict2_tbl <- left_join(annot_dict2_tbl, count_tbl, by=c("gene_id"="transcript_id"))

cor(annot_dict2_tbl$read_n, annot_dict2_tbl$uniquely_mapped_reads_mean,use="complete.obs")
cor(annot_dict2_tbl$read_n, annot_dict2_tbl$mean_count,use="complete.obs")
cor(annot_dict2_tbl$uniquely_mapped_reads_mean, annot_dict2_tbl$mean_count,use="complete.obs")

write_tsv(annot_dict2_tbl, "CopciAB_ref_20210415_dict2.tsv")

annot_dict3_tbl <- annot_dict2_tbl
annot_dict3_tbl$annot_index <- !is.na(annot_dict3_tbl$Analisys)
annot_dict3_tbl <- dplyr::select(annot_dict3_tbl, -Analisys, -Signature_accession, -InterPro_annotations_accession, -InterPro_annotations_description)

filter(annot_dict3_tbl, str_detect(gene_id, pattern = "000132"))
annot_dict3_tbl <- dplyr::rename(annot_dict3_tbl, fullmatch_longred_n = read_n)

annot_dict3_tbl$ont_index <- sapply(str_split(annot_dict3_tbl$gene_id, pattern="_|\\.|[0-9]+"), "[[", 1) %in% c("G", "man", "PB")

table(annot_dict3_tbl$ont_index[is.na(annot_dict3_tbl$fullmatch_longred_n)])

hist(annot_dict3_tbl[annot_dict3_tbl$ont_index,]$zero_count_p)
hist(annot_dict3_tbl[!annot_dict3_tbl$ont_index,]$zero_count_p)

ggplot(annot_dict3_tbl, aes(x = zero_count_p)) + 
  geom_histogram() +
  facet_wrap(~ont_index)

ggplot(annot_dict3_tbl, aes(x = mean_count)) + 
  geom_histogram(bins = 1000) +
  coord_cartesian(xlim = c(0,10000)) +
  facet_wrap(~ont_index)

annot_dict3_tbl %>% 
  group_by(ont_index) %>% 
  summarise(fullmatch_longred_n_mean = mean(fullmatch_longred_n, na.rm=T),
            mean_count_mean = mean(mean_count, na.rm = T),
            max_count_mean = mean(max_count, na.rm = T),
            zero_count_p_mean = mean(zero_count_p, na.rm = T))


annot_dict3_tbl %>% 
  filter(exon_n > 1) %>% 
  dplyr::select(gene_id,ont_index, fullmatch_longred_n, uniquely_mapped_reads_min, mean_count, zero_count_p) %>% 
  mutate(fullmatch_longred_n = ifelse(is.na(fullmatch_longred_n), 0, fullmatch_longred_n),
         mean_count = ifelse(is.na(mean_count), 0, mean_count),
         zero_count_p = ifelse(is.na(zero_count_p), 1, zero_count_p),
         uniquely_mapped_reads_min = as.numeric(uniquely_mapped_reads_min),
         uniquely_mapped_reads_min = ifelse(is.infinite(uniquely_mapped_reads_min), 0, uniquely_mapped_reads_min)) -> annot_dict_cluster

library(FactoMineR)
library(factoextra)
library(corrplot)

annot_dict_cluster_scale <- scale(annot_dict_cluster[,c(-1, -2)])
annot_dict_cluster_scale_PCA <- PCA(annot_dict_cluster_scale, scale.unit=F, graph=T)

factoextra::fviz_pca_var(annot_dict_cluster_scale_PCA, axes=1:2)


annot_dict_cluster_scale_cor <- cor(annot_dict_cluster_scale) 
corrplot(annot_dict_cluster_scale_cor)

annot_dict_cluster_scale_kmean <- kmeans(annot_dict_cluster_scale, centers = 2)
annot_dict_cluster$kmean <- annot_dict_cluster_scale_kmean$cluster

count(annot_dict_cluster, ont_index, kmean)
filter(annot_dict_cluster, ont_index==F, kmean == 1)

filter(annot_dict_cluster, ont_index==T, kmean == 2)
filter(annot_dict_cluster, ont_index==T, kmean == 1)

factoextra::fviz_pca_ind(annot_dict_cluster_scale_PCA, axes=1:2, col.ind=annot_dict_cluster$kmean)
factoextra::fviz_pca_ind(annot_dict_cluster_scale_PCA, axes=1:2, col.ind=annot_dict_cluster$ont_index)


############################# mmseqs RBH ----

# mkdir query_db
# mmseqs createdb ../../CopciAB_new_annot_CDS_20220127.prot.fasta query_db/CopciAB_20220127

# mkdir target_db
# mmseqs createdb Copci_AmutBmut1_GeneModels_FrozenGeneCatalog_20160912_aa_cleanname.fasta target_db/JGI

# mmseqs rbh 
# mmseqs rbh query_db/CopciAB_20220127 target_db/JGI new_vs_jgi ./tmp --threads 50

# mmseqs convertalis --format-output query,target,pident,evalue,bits,qstart,qend,tstart,tend,qlen,tlen,nident query_db/CopciAB_20220127 target_db/JGI new_vs_jgi new_vs_jgi_RBH.tsv

# qir = nident/qlen

rbh_cnames <- c("query,target,pident,evalue,bits,qstart,qend,tstart,tend,qlen,tlen,nident")
rbh_cnames <- unlist(str_split(rbh_cnames, pattern = ","))
rbh_tbl <- read_tsv("../../../mmseqs_rbh_dir/annot_20220127_dir/new_vs_jgi_RBH.tsv", col_names=rbh_cnames)

sum(duplicated(rbh_tbl$query)) # 9
sum(duplicated(rbh_tbl$target)) # 213

filter(rbh_tbl, query %in% unique(rbh_tbl$query[duplicated(rbh_tbl$query)])) %>% arrange(query)
filter(rbh_tbl, target %in% unique(rbh_tbl$target[duplicated(rbh_tbl$target)])) %>% arrange(target)

rbh_short_tbl <- select(rbh_tbl, query, target, pident, bits, qlen, tlen)
rbh_short_tbl$source <- "mmseqs_rbh"
rbh_colnames <- "query, target, pident, bits, qlen, tlen"
rbh_colnames <- unlist(str_split(rbh_colnames, pattern=", "))
JGI_short_tbl <- select(JGI_blast_max_tbl, qseqid, sseqid, pident, bitscore, qlen, length)
colnames(JGI_short_tbl) <- rbh_colnames
JGI_short_tbl$source <- "blastp"

comper_data_tbl <- bind_rows(rbh_short_tbl, JGI_short_tbl)
comper_data_tbl <- arrange(comper_data_tbl, query)

comper_data_lt <- split(comper_data_tbl, comper_data_tbl$query)
table(sapply(comper_data_lt, nrow))

table(sapply(comper_data_lt, nrow))
comper_data_short_lt <- comper_data_lt[sapply(comper_data_lt, nrow)==1]



extractClustersFromSelfHits <- function(hits) {
  stopifnot(is(hits, "Hits"))
  N <- queryLength(hits)
  stopifnot(N == subjectLength(hits))
  h <- GenomicRanges::union(hits, t(hits))
  qh <- queryHits(h)
  sh <- subjectHits(h)
  cid <- cid0 <- seq_len(N) 
  while (TRUE) {
    cid2 <- pmin(cid, selectHits(h, "first"))
    if (identical(cid2, cid))
      break
    cid <- cid2
    h <- Hits(qh, cid[sh], N, N)
  }
  unname(splitAsList(cid0, cid))
}


mergeConnectedRanges <- function(x, hits) {
  stopifnot(is(x, "GenomicRanges"))
  stopifnot(is(hits, "Hits"))
  stopifnot(queryLength(hits) == subjectLength(hits))
  stopifnot(queryLength(hits) == length(x))
  clusters <- extractClustersFromSelfHits(hits)
  ans <- range(extractList(x, clusters))
  if (any(elementNROWS(ans) != 1L))
    stop(wmsg("some connected ranges are not on the same ", "chromosome and strand, and thus cannot be ", "merged"))
  ans <- unlist(ans)
  mcols(ans)$revmap <- clusters
  ans
}



annot_tscripts_gr <- annot_gr[annot_gr$type == "transcript"]

hits <- findOverlaps(annot_tscripts_gr, drop.self = T) 

x <- annot_tscripts_gr[queryHits(hits)]
y <- annot_tscripts_gr[subjectHits(hits)]
relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))
hits <- hits[relative_overlap >= 0.6] 
gr1 <- mergeConnectedRanges(annot_tscripts_gr, hits)
gr1_df <- as.data.frame(gr1)
gr1_df <- cbind(gr1_df, plyr::ldply(gr1$revmap, rbind))
gr1_df <- select(gr1_df, -strand, -revmap)

gr1_df <- arrange(gr1_df, 
                  as.integer(str_extract(as.character(gr1_df$seqnames), pattern="[0-9]+")),
                  start)

gr1_df$overlap_index <- seq_len(nrow(gr1_df))

gr2_df <- tidyr::gather(gr1_df, 5:38, key=overlap_n, value=annot_index)
gr2_df <- filter(gr2_df, !is.na(annot_index)) 
gr2_df <- arrange(gr2_df, overlap_index)
gr2_df$transcript_id <- annot_tscripts_gr$gene_id[gr2_df$annot_index]

gr2_df_extra <- left_join(gr2_df, annot_dict2_tbl, by=c("transcript_id"="gene_id"))

gr2_df_extra$uniquely_mapped_reads_mean[is.na(gr2_df_extra$uniquely_mapped_reads_mean)] <- 0
any(is.na(gr2_df_extra$uniquely_mapped_reads_mean))
any(is.infinite(gr2_df_extra$uniquely_mapped_reads_mean))

gr2_df_extra$uniquely_mapped_reads_sd[is.na(gr2_df_extra$uniquely_mapped_reads_sd)] <- 0
any(is.na(gr2_df_extra$uniquely_mapped_reads_sd))
any(is.infinite(gr2_df_extra$uniquely_mapped_reads_sd))

gr2_df_extra$uniquely_mapped_reads_min[is.na(gr2_df_extra$uniquely_mapped_reads_min)] <- 0
gr2_df_extra$uniquely_mapped_reads_min[is.infinite(gr2_df_extra$uniquely_mapped_reads_min)] <- 0
any(is.na(gr2_df_extra$uniquely_mapped_reads_sd))
any(is.infinite(gr2_df_extra$uniquely_mapped_reads_sd))

gr2_df_extra$read_n[is.na(gr2_df_extra$read_n)] <- 0

table(sapply(str_split(gr2_df_extra$transcript_id, pattern="_|\\.|[0-9]+"), "[[", 1))
gr2_df_extra$ont_index <- sapply(str_split(gr2_df_extra$transcript_id, pattern="_|\\.|[0-9]+"), "[[", 1) %in% c("G", "man", "PB")
  

table(sapply(str_split(filter(gr2_df_extra, is.na(protein_id))$transcript_id, pattern="_|\\.|[0-9]+"), "[[", 1))
filter(gr2_df_extra, is.na(protein_id))$transcript_id[sapply(str_split(filter(gr2_df_extra, is.na(protein_id))$transcript_id, pattern="_|\\.|[0-9]+"), "[[", 1)=="G"]


gr2_df_extra_coding <- filter(gr2_df_extra, !is.na(protein_id))
gr2_df_extra_coding_lt <- split(gr2_df_extra_coding, gr2_df_extra_coding$overlap_index) 

sum(sapply(gr2_df_extra_coding_lt, nrow) == 1)

sum(sapply(gr2_df_extra_coding_lt, function(x) all(x$exon_n == 1)))
sum(sapply(gr2_df_extra_coding_lt, function(x) all(x$exon_n == 1) & length(x$width_aa == max(x$width_aa)) > 1))
sum(sapply(gr2_df_extra_coding_lt, function(x) all(x$ont_index)))

gr2_df_extra_coding_lt[sapply(gr2_df_extra_coding_lt, function(x) all(x$ont_index))]

filtered_overlap_lt <- lapply(gr2_df_extra_coding_lt, function(x) {
  if(nrow(x) == 1) {
    result <- x
  } else if(all(x$ont_index)) {
    result_index <- x$width_aa == max(x$width_aa) & x$read_n == max(x$read_n) & x$uniquely_mapped_reads_mean == max(x$uniquely_mapped_reads_mean)
    if(all(result_index)) {
      result <- x[1,]
    } else if(!all(result_index)) {
      result <- arrange(x[x$read_n == max(x$read_n),], desc(width_aa))[1,]
    } else {
      result <- x[result_index,][1,]
    }
 } else if(any(x$ont_index)) {
   result <- arrange(x[x$ont_index,], desc(width_aa))[1,]
  } else if(all(x$exon_n == 1)) {
    result <- x[x$width_aa == max(x$width_aa), ][1,]
  } else if(any(x$exon_n == 1)) {
    result <- x[x$width_aa == max(x$width_aa), ][1,]
  } else {
    result_index2 <- x$width_aa == max(x$width_aa) & x$read_n == max(x$read_n) & x$uniquely_mapped_reads_mean == max(x$uniquely_mapped_reads_mean)
    if(all(result_index2)) {
      result <- x[result_index2,][1,]
    } else {
      result <- arrange(x, desc(width_aa))[1,]
    }
  }
  return(result)
})

table(sapply(filtered_overlap_lt, nrow))
sum(sapply(filtered_overlap_lt, function(x) any(is.na(x$ont_index))))

filtered_overlap_tbl <- bind_rows(filtered_overlap_lt)

removable_overlap_index <- gr2_df_extra$transcript_id[!gr2_df_extra$transcript_id %in% filtered_overlap_tbl$transcript_id]

annot_dict2_tbl$removable_overlap <- annot_dict2_tbl$gene_id %in% removable_overlap_index


annot_cds_gr <- annot_gr[annot_gr$type == "CDS"]
annot_cds_df <- as.data.frame(annot_cds_gr)
annot_cds_lt <- split(annot_cds_df, annot_cds_df$gene_id)
annot_cds_lt <- lapply(annot_cds_lt, function(x) {
  result <- x[1,]
  result$end <- x$end[x$end == max(x$end)]
  return(result)
})
annot_cds_df <- bind_rows(annot_cds_lt)
annot_cds_gr <- makeGRangesFromDataFrame(annot_cds_df, keep.extra.columns = T)

hits <- findOverlaps(annot_cds_gr, drop.self = T) 

x <- annot_cds_gr[queryHits(hits)]
y <- annot_cds_gr[subjectHits(hits)]
relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))

gr1_cds <- mergeConnectedRanges(annot_cds_gr, hits)
gr1_cds_df <- as.data.frame(gr1_cds)
gr1_cds_df <- cbind(gr1_cds_df, plyr::ldply(gr1_cds$revmap, rbind))
gr1_cds_df <- dplyr::select(gr1_cds_df, -strand, -revmap)

gr1_cds_df <- arrange(gr1_cds_df, 
                      as.integer(str_extract(as.character(gr1_cds_df$seqnames), pattern="[0-9]+")),
                      start)

gr1_cds_df$overlap_index <- seq_len(nrow(gr1_cds_df))

gr2_cds_df <- tidyr::gather(gr1_cds_df, 5:26, key=overlap_n, value=annot_index)
gr2_cds_df <- filter(gr2_cds_df, !is.na(annot_index))
gr2_cds_df <- arrange(gr2_cds_df, overlap_index)
gr2_cds_df$transcript_id <- annot_cds_gr$gene_id[gr2_cds_df$annot_index]

gr2_cds_df_extra <- left_join(gr2_cds_df, annot_dict2_tbl, by=c("transcript_id"="gene_id"))

gr2_cds_df_extra$uniquely_mapped_reads_mean[is.na(gr2_cds_df_extra$uniquely_mapped_reads_mean)] <- 0
any(is.na(gr2_cds_df_extra$uniquely_mapped_reads_mean))
any(is.infinite(gr2_cds_df_extra$uniquely_mapped_reads_mean))

gr2_cds_df_extra$uniquely_mapped_reads_sd[is.na(gr2_cds_df_extra$uniquely_mapped_reads_sd)] <- 0
any(is.na(gr2_cds_df_extra$uniquely_mapped_reads_sd))
any(is.infinite(gr2_cds_df_extra$uniquely_mapped_reads_sd))

gr2_cds_df_extra$uniquely_mapped_reads_min[is.na(gr2_cds_df_extra$uniquely_mapped_reads_min)] <- 0
gr2_cds_df_extra$uniquely_mapped_reads_min[is.infinite(gr2_cds_df_extra$uniquely_mapped_reads_min)] <- 0
any(is.na(gr2_cds_df_extra$uniquely_mapped_reads_sd))
any(is.infinite(gr2_cds_df_extra$uniquely_mapped_reads_sd))

gr2_cds_df_extra$read_n[is.na(gr2_cds_df_extra$read_n)] <- 0

table(sapply(str_split(gr2_cds_df_extra$transcript_id, pattern="_|\\.|[0-9]+"), "[[", 1))
gr2_cds_df_extra$ont_index <- sapply(str_split(gr2_cds_df_extra$transcript_id, pattern="_|\\.|[0-9]+"), "[[", 1) %in% c("G", "man", "PB")

gr2_cds_df_extra_coding_lt <- split(gr2_cds_df_extra, gr2_cds_df_extra$overlap_index) 

sum(sapply(gr2_cds_df_extra_coding_lt, nrow) == 1)

sum(sapply(gr2_cds_df_extra_coding_lt, function(x) all(x$exon_n == 1)))
sum(sapply(gr2_cds_df_extra_coding_lt, function(x) all(x$exon_n == 1) & length(x$width_aa == max(x$width_aa)) > 1))
sum(sapply(gr2_cds_df_extra_coding_lt, function(x) all(x$ont_index)))

gr2_cds_df_extra_coding_lt[sapply(gr2_cds_df_extra_coding_lt, function(x) all(x$ont_index))]
gr2_cds_df_extra_coding_lt[sapply(gr2_cds_df_extra_coding_lt, function(x) any(x$protein_id == "G10501.1_t0_r1"))]

filtered_cds_overlap_lt <- lapply(gr2_cds_df_extra_coding_lt, function(x) {
  if(nrow(x) == 1) {
    result <- x
  } else if(all(x$ont_index)) {
    result_index <- x$width_aa == max(x$width_aa) & x$read_n == max(x$read_n) & x$uniquely_mapped_reads_mean == max(x$uniquely_mapped_reads_mean)
    if(all(result_index)) {
      result <- x[1,]
    } else if(!all(result_index)) {
      result <- arrange(x[x$read_n == max(x$read_n),], desc(width_aa))[1,]
    } else {
      result <- x[result_index,][1,]
    }
  } else if(any(x$ont_index)) {
    result <- arrange(x[x$ont_index,], desc(width_aa))[1,]
  } else if(all(x$exon_n == 1)) {
    result <- x[x$width_aa == max(x$width_aa), ][1,]
  } else if(any(x$exon_n == 1)) {
    result <- x[x$width_aa == max(x$width_aa), ][1,]
  } else {
    result_index2 <- x$width_aa == max(x$width_aa) & x$read_n == max(x$read_n) & x$uniquely_mapped_reads_mean == max(x$uniquely_mapped_reads_mean)
    if(all(result_index2)) {
      result <- x[result_index2,][1,]
    } else {
      result <- arrange(x, desc(width_aa))[1,]
    }
  }
  return(result)
})


filtered_cds_overlap_lt[sapply(filtered_cds_overlap_lt, function(x) any(x$protein_id == "G10501.1_t0_r1"))]


sum(sapply(filtered_cds_overlap_lt, function(x) any(is.na(x$ont_index))))

filtered_cds_overlap_tbl <- bind_rows(filtered_cds_overlap_lt)

removable_cds_overlap_index <- gr2_cds_df_extra$transcript_id[!gr2_cds_df_extra$transcript_id %in% filtered_cds_overlap_tbl$transcript_id]

annot_dict2_tbl$removable_cds_overlap <- annot_dict2_tbl$gene_id %in% removable_cds_overlap_index

sum(annot_dict2_tbl$removable_overlap & annot_dict2_tbl$removable_cds_overlap) 
sum(!annot_dict2_tbl$removable_overlap & annot_dict2_tbl$removable_cds_overlap) 
sum(annot_dict2_tbl$removable_overlap & !annot_dict2_tbl$removable_cds_overlap) 

annot_dict2_tbl$protein_id[annot_dict2_tbl$removable_overlap & !annot_dict2_tbl$removable_cds_overlap]
str_subset(na.omit(annot_dict2_tbl$protein_id[annot_dict2_tbl$removable_overlap & !annot_dict2_tbl$removable_cds_overlap]), pattern="^G|^PB|^man")

sum(annot_dict2_tbl$removable_overlap | annot_dict2_tbl$removable_cds_overlap & annot_dict2_tbl$gene_id != "G12835.5_t0_r1")
sum(annot_dict2_tbl$removable_overlap | annot_dict2_tbl$removable_cds_overlap)

annot_dict2_tbl$removable_clean_overlap <- annot_dict2_tbl$removable_overlap | annot_dict2_tbl$removable_cds_overlap
annot_dict2_tbl$removable_clean_overlap[annot_dict2_tbl$gene_id == "G12835.5_t0_r1"] <- FALSE


annot_dict2_tbl$ont_index <- sapply(str_split(annot_dict2_tbl$gene_id, pattern="_|\\.|[0-9]+"), "[[", 1) %in% c("G", "man", "PB")

write_tsv(annot_dict2_tbl, "CopciAB_ref_20210415_dict3.tsv")


annot_nooverlap_gr <- annot_gr[annot_gr$gene_id %in% annot_dict2_tbl$gene_id[!annot_dict2_tbl$removable_clean_overlap]]

export(annot_nooverlap_gr, "../../../CopciAB_new_annot_CDS_20220127_nooverlap2.gtf")


library(GenomicFeatures)
library(BSgenome)

ref_genome <- readDNAStringSet("../../../../../ref_genome/fullgenome_jgi/merged_genome/CopciAB_new_jgi_20220113.fasta")

new_annot_CDS_TxDb <- makeTxDbFromGRanges(annot_nooverlap_gr)

new_annot_CDS_TxDb_lt <- GenomicFeatures::as.list(new_annot_CDS_TxDb)
tx_CDS_name <- new_annot_CDS_TxDb_lt$transcripts$tx_name

cbt <- cdsBy(new_annot_CDS_TxDb, by = "tx")
names(cbt) <- tx_CDS_name[as.integer(names(cbt))]

cds_seq <- extractTranscriptSeqs(ref_genome, cbt)
prot_seq <- translate(cds_seq)

writeXStringSet(prot_seq, "CopciAB_new_annot_CDS_20220127_nooverlap2_prot.fasta")

prot_path <- "CopciAB_new_annot_CDS_20220127_nooverlap2_prot.fasta"
BUSCO_cmd <- str_c("run_BUSCO.py -f -i ", prot_path,
                   " -o ", str_replace(prot_path, pattern = ".prot.fasta", replacement = "_odb10"),
                   " -l /work/balintb/Busco/Databases/v4/basidiomycota_odb10/ -c 10 -m prot")
system(command = BUSCO_cmd)

annot_gr <- import("../../../CopciAB_new_annot_CDS_20220127_nooverlap2.gtf")


annot_cds_gr <- annot_gr[annot_gr$type == "CDS"]
annot_cds_df <- as.data.frame(annot_cds_gr)
annot_cds_lt <- split(annot_cds_df, annot_cds_df$gene_id)
annot_cds_lt <- lapply(annot_cds_lt, function(x) {
  result <- x[1,]
  result$end <- x$end[x$end == max(x$end)]
  return(result)
})
annot_cds_df <- bind_rows(annot_cds_lt)
annot_cds_gr <- makeGRangesFromDataFrame(annot_cds_df, keep.extra.columns = T)

hits <- findOverlaps(annot_cds_gr, drop.self = T, drop.redundant=T, ignore.strand=T) 

q_gene_id <- annot_cds_gr[queryHits(hits)]$gene_id
q_start <- start(annot_cds_gr[queryHits(hits)])
q_seqnames <- as.character(seqnames(annot_cds_gr[queryHits(hits)]))
h_gene_id <- annot_cds_gr[subjectHits(hits)]$gene_id

cds_olap_tbl <- tibble(querySeqnames = q_seqnames,
                       queryStart = q_start,
                       queryHits = q_gene_id,
                       subjectHits = h_gene_id)

cds_olap_lt <- split(cds_olap_tbl, cds_olap_tbl$queryHits)
cds_olap_spresd_lt <- lapply(cds_olap_lt, function(x) {
  x$index <- seq_len(nrow(x))
  x <- spread(x, key=index, value=subjectHits)
  return(x)
})
cds_olap_spresd_tbl <- bind_rows(cds_olap_spresd_lt)
cds_olap_spresd_tbl %>% 
  arrange(as.integer(str_extract(querySeqnames, pattern = "[0-9]+$")), queryStart) -> cds_olap_spresd_tbl

write_tsv(cds_olap_spresd_tbl, "CopciAB_ref_20210415_cds_olap_man.tsv")

library(readxl)
esheets <- excel_sheets("CopciAB_ref_20210415_cds_olap_man.xlsx")
cds_olap_spresd_man_tbl <- read_excel("CopciAB_ref_20210415_cds_olap_man.xlsx", esheets[[1]])
id_all <- unlist(cds_olap_spresd_man_tbl[,c(3:4)], use.names = F)
id_rem <- cds_olap_spresd_man_tbl$index_marad
id_rem <- str_trim(unlist(str_split(id_rem, pattern=";")))
szemet <- id_all[!id_all %in% id_rem]


annot_dict2_tbl <- read_tsv("CopciAB_ref_20210415_dict3.tsv")
annot_dict2_tbl$removable_cds_overlap2 <- annot_dict2_tbl$gene_id %in% szemet
annot_dict2_tbl$removable_clean_overlap2 <- annot_dict2_tbl$removable_clean_overlap | annot_dict2_tbl$removable_cds_overlap2
annot_dict2_tbl$cds_index <- is.na(annot_dict2_tbl$width_aa)
annot_dict2_tbl$removable_clean_overlap2_onlycoding <- annot_dict2_tbl$removable_clean_overlap2 | annot_dict2_tbl$cds_index

write_tsv(annot_dict2_tbl, "CopciAB_ref_20210415_dict4.tsv")


annot_nooverlap_onlycds_gr <- annot_gr[annot_gr$gene_id %in% annot_dict2_tbl$gene_id[!annot_dict2_tbl$removable_clean_overlap2_onlycoding]]

export(annot_nooverlap_onlycds_gr, "../../../CopciAB_new_annot_CDS_20220127_nooverlap3_onlycds.gtf")


library(GenomicFeatures)
library(BSgenome)

ref_genome <- readDNAStringSet("../../../../../ref_genome/fullgenome_jgi/merged_genome/CopciAB_new_jgi_20220113.fasta")

new_annot_CDS_TxDb <- makeTxDbFromGRanges(annot_nooverlap_onlycds_gr)

new_annot_CDS_TxDb_lt <- GenomicFeatures::as.list(new_annot_CDS_TxDb)
tx_CDS_name <- new_annot_CDS_TxDb_lt$transcripts$tx_name

cbt <- cdsBy(new_annot_CDS_TxDb, by = "tx")
names(cbt) <- tx_CDS_name[as.integer(names(cbt))]

cds_seq <- extractTranscriptSeqs(ref_genome, cbt)
prot_seq <- translate(cds_seq)

writeXStringSet(prot_seq, "CopciAB_new_annot_CDS_20220127_nooverlap3_onlycds_prot.fasta")

prot_path <- "CopciAB_new_annot_CDS_20220127_nooverlap3_onlycds_prot.fasta"
BUSCO_cmd <- str_c("run_BUSCO.py -f -i ", prot_path,
                   " -o ", str_replace(prot_path, pattern = ".prot.fasta", replacement = "_odb10"),
                   " -l /work/balintb/Busco/Databases/v4/basidiomycota_odb10/ -c 10 -m prot")
system(command = BUSCO_cmd)

library(ggplot2)

annot_UTR_tbl <- as.data.frame(annot_gr[annot_gr$type %in% c("3UTR", "5UTR")])
annot_UTR_tbl %>% 
  group_by(gene_id, type = as.character(type)) %>% 
  summarise(utr_length = sum(width)) %>% 
  ungroup() -> UTR_length_tbl

UTR_length_tbl %>% 
  group_by(type) %>% 
  summarise(mean_UTR_length = mean(utr_length),
            sd_UTR_length = sd(utr_length),
            UTR_n = n()) %>% 
ggplot(aes(type, mean_UTR_length)) +
  geom_col()

UTR_length_tbl %>% 
  ggplot(aes(type, utr_length, fill=type)) +
  geom_boxplot()


annot_dict3_tbl %>% 
  filter(!removable_overlap) %>% 
  dplyr::select(gene_id, qseqid, pident, qir) -> jgi_id_test_short

jgi_id_test_short %>% 
  filter(str_detect(gene_id, pattern = "^G"),
         !is.na(pident),
         pident == 100,
         !is.na(qir),
         sapply(str_split(qir, pattern = " \\| "), function(y) any(as.numeric(y) <= 0.8)))

jgi_id_test_short_qseqid_lt <- split(jgi_id_test_short, jgi_id_test_short$qseqid)
jgi_id_test_short_qseqid_lt[sapply(jgi_id_test_short_qseqid_lt, nrow) > 1]

