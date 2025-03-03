library(dplyr)
library(stringr)
library(readr)
library(tibble)
library(purrr)
library(tidyr)
library(ggplot2)
library(rtracklayer)

options(width=300, scipen = 999)

wd <- "/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/clean_annotation_20220427_dir/filter_by_RNAseq_readsupport"
setwd(wd)

group_lt <- list(
  qseq_gp = read_tsv("/node9_data/bhegedus/pjx_Copci_Quantseq_20220323/additionals/sample_path_and_grouping_info_plus_new.tsv", col_select = c("sample", "sample_group")),
  jgi_gp = read_tsv("/node9_data/bhegedus/pjx_Copci_JGI_20220401/additionals/sample_path_and_grouping_info_plus_new.tsv", col_select = c("Sample_Name", "group_name")),
  jgiKrizsi_gp = read_tsv("/node9_data/bhegedus/pjx_Copci_KrizsiJGI_20220405/additionals/sample_path_and_grouping_info_plus_new.tsv", col_select = c("sample_raw", "sample_group")),
  novo_gp = read_tsv("/node9_data/bhegedus/pjx_Copci_Novogene_20220329/additionals/sample_path_and_grouping_info_plus_new.tsv", col_select = c("sample_raw", "sample_group")),
  xie_gp = read_tsv("/node9_data/bhegedus/pjx_Copci_Xie_20220406/additionals/sample_path_and_grouping_info_plus_new.tsv", col_select = c("sample_name", "sample_group"))
)

group_lt <- lapply(group_lt, setNames, c("sample_name", "sample_group"))
group_tbl <- bind_rows(group_lt, .id="experiments")
group_tbl <- mutate(group_tbl, id_index = str_c(sample_name, str_remove(experiments, pattern = "_gp"), sep = "_"))
group_tbl <- mutate(group_tbl, group_index = str_c(sample_group, sapply(str_split(id_index, pattern = "_"), function(x) rev(x)[[1]]), sep = "_"))
group_tbl <- mutate(group_tbl, group_index2 = sapply(str_split(id_index, pattern = "_"), function(x) rev(x)[[1]]))
group_tbl <- mutate(group_tbl, group_index3 = ifelse(group_index2=="qseq", group_index2, "RNAseq"))

cpm_lt <- list(
  qseq_cpm = read_tsv("/node9_data/bhegedus/pjx_Copci_Quantseq_20220323/DEG_analysis/CopciAB_quantseq_CPM_20220427.tsv"),
  jgi_cpm = read_tsv("/node9_data/bhegedus/pjx_Copci_JGI_20220401/DEG_analysis/CopciAB_quantseq_CPM_20220427.tsv"),
  jgiKrizsi_cpm = read_tsv("/node9_data/bhegedus/pjx_Copci_KrizsiJGI_20220405/DEG_analysis/CopciAB_quantseq_CPM_20220427.tsv"),
  novo_cpm = read_tsv("/node9_data/bhegedus/pjx_Copci_Novogene_20220329/DEG_analysis/CopciAB_quantseq_CPM_20220427.tsv"),
  xie_cpm = read_tsv("/node9_data/bhegedus/pjx_Copci_Xie_20220406/DEG_analysis/CopciAB_quantseq_CPM_20220428.tsv")
)

cpm_lt <- lapply(cpm_lt, function(x) column_to_rownames(x, "gene_id"))

for(i in names(cpm_lt)) {
  cnames <- names(cpm_lt[[i]])
  names(cpm_lt[[i]]) <- str_c(cnames, str_remove(i, pattern = "_cpm"), sep = "_")
}

cpm_lt <- lapply(cpm_lt, function(x) rownames_to_column(x, "gene_id"))

cpm_tbl <- cpm_lt %>% purrr::reduce(full_join, by="gene_id")
cpm_tbl <- column_to_rownames(cpm_tbl, "gene_id")


dim(cpm_tbl)
setdiff(group_tbl$id_index, names(cpm_tbl))
setdiff(names(cpm_tbl), group_tbl$id_index)


cpm_t <- t(cpm_tbl)
cpm_t_tbl <- as.data.frame(cpm_t)

group_index <- group_tbl$group_index2[match(rownames(cpm_t_tbl), group_tbl$id_index)]

cpm_t_lt <- split(cpm_t_tbl, group_index)

sapply(cpm_t_lt, nrow)
rname_check_lt <- lapply(cpm_t_lt, function(x) {
  result <- tibble(rnames = rownames(x))
  return(result)
})
rname_check_tbl <- bind_rows(rname_check_lt, .id="orname") 

cpm_mean_lt <- lapply(cpm_t_lt, function(x) {
  return(as.data.frame(colMeans(x)))
})
cpm_mean_lt <- imap(cpm_mean_lt, function(x, y) setNames(x, y))
cpm_mean_lt <- lapply(cpm_mean_lt, function(x) return(rownames_to_column(x, var = "gene_id")))
cpm_mean_tbl <- cpm_mean_lt %>% purrr::reduce(full_join, by="gene_id")
cpm_mean_tbl <- column_to_rownames(cpm_mean_tbl, var = "gene_id")

cpm_mean_sort_tbl <- arrange(cpm_mean_tbl, jgiKrizsi)
cpm_mean_sort_tbl <- rownames_to_column(cpm_mean_sort_tbl, var = "gene_id")
cpm_mean_sort_g_tbl <- gather(cpm_mean_sort_tbl, key="seqType", value="CPM", -gene_id)


ggplot(cpm_mean_sort_g_tbl, aes(x=log2(CPM+0.25), colour=seqType)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="red")


cpm_mean_sort_g_lt <- split(cpm_mean_sort_g_tbl, cpm_mean_sort_g_tbl$seqType)
cpm_mean_sort_g_f1cpm_lt <- lapply(cpm_mean_sort_g_lt, function(x) return(filter(x, CPM<=1)))

barplot(sapply(cpm_mean_sort_g_f1cpm_lt, nrow), main = "mean gene count <= 1CPM")

library(UpSetR)
cpm_for_ipset_lt <- lapply(cpm_mean_sort_g_f1cpm_lt, "[[", "gene_id")

upset(fromList(cpm_for_ipset_lt), order.by = "freq")


group_index <- group_tbl$group_index[match(rownames(cpm_t_tbl), group_tbl$id_index)]

cpm_t_lt <- split(cpm_t_tbl, group_index)


table(sapply(cpm_t_lt, nrow))


cpm_mean_lt <- lapply(cpm_t_lt, function(x) {
  return(as.data.frame(colMeans(x)))
})
cpm_mean_lt <- imap(cpm_mean_lt, function(x, y) setNames(x, y))
cpm_mean_lt <- lapply(cpm_mean_lt, function(x) return(rownames_to_column(x, var = "gene_id")))
cpm_mean_tbl <- cpm_mean_lt %>% purrr::reduce(full_join, by="gene_id")
cpm_mean_tbl <- column_to_rownames(cpm_mean_tbl, var = "gene_id")


cpm_mean_t_tbl <- as.data.frame(t(cpm_mean_tbl))

group_index2 <- group_tbl$group_index2[match(rownames(cpm_mean_t_tbl), group_tbl$group_index)]

cpm_mean_t_lt <- split(cpm_mean_t_tbl, group_index2)


table(sapply(cpm_mean_t_lt, nrow))

sum(sapply(cpm_mean_t_lt, nrow))
sapply(cpm_mean_t_lt, nrow)


cpm_max_lt <- lapply(cpm_mean_t_lt, function(x) {
  return(as.data.frame(apply(x,2,max)))
})

cpm_max_lt <- imap(cpm_max_lt, function(x, y) setNames(x, y))
cpm_max_lt <- lapply(cpm_max_lt, function(x) return(rownames_to_column(x, var = "gene_id")))
cpm_max_tbl <- cpm_max_lt %>% purrr::reduce(full_join, by="gene_id")

cpm_max_g_tbl <- gather(cpm_max_tbl, key="seqType", value="CPM", -gene_id)

ggplot(cpm_max_g_tbl, aes(x=log2(CPM+0.25), colour=seqType)) +
  geom_density() +
  geom_vline(xintercept = 0, colour="red")


cpm_max_g_lt <- split(cpm_max_g_tbl, cpm_max_g_tbl$seqType)
cpm_max_g_f1cpm_lt <- lapply(cpm_max_g_lt, function(x) return(filter(x, CPM<=1)))

barplot(sapply(cpm_max_g_f1cpm_lt, nrow), main = "max mean gene count <= 1CPM")

library(UpSetR)
cpm_for_upset_lt <- lapply(cpm_max_g_f1cpm_lt, "[[", "gene_id")

upset(fromList(cpm_for_upset_lt), order.by = "freq")


cpm_max_g_lt <- split(cpm_max_g_tbl, cpm_max_g_tbl$seqType)

only_pcpm <- sapply(cpm_max_g_lt, function(x) {
  cpm_q <- quantile(x$CPM, probs=0.25)
  return(cpm_q)
})
only_pcpm_log2 <- log2(only_pcpm+0.25) 

ggplot(cpm_max_g_tbl, aes(x=log2(CPM+0.25), colour=seqType)) +
  geom_density() +
  geom_vline(xintercept = only_pcpm)


cpm_max_g_pcpm_lt <- lapply(cpm_max_g_lt, function(x) {
  cpm_q <- quantile(x$CPM, probs=0.25)
  return(filter(x, CPM<=cpm_q))
})

barplot(sapply(cpm_max_g_pcpm_lt, nrow), main = "max mean gene count <= quantile 25%")

library(UpSetR)
library(ComplexHeatmap)

cpm_for_upset_lt <- lapply(cpm_max_g_pcpm_lt, "[[", "gene_id")

upset(fromList(cpm_for_upset_lt), order.by = "freq")

cpm_cbat_matrix <- make_comb_mat(cpm_for_upset_lt)

comb_name(cpm_cbat_matrix, readable = TRUE)
comb_name(cpm_cbat_matrix)
comb_index <- extract_comb(cpm_cbat_matrix, "11111") 
length(comb_index)


library(Biostrings)

prot_seq <- readAAStringSet("../CopciAB_new_annot_CDS_20220425m_newIds_isoforms_mancur_prot.fasta")
sprot_seq <- prot_seq[names(prot_seq) %in% comb_index]

writeXStringSet(sprot_seq, "CopciAB_new_annot_CDS_20220425m_RNAseq_filtered.fasta")


ipro_tbl <- read_tsv("../interproscan88_out/CopciAB_new_annot_CDS_20220425m_iproscan88_headered.txt",na = "")
ipro_wiproid_tbl <- filter(ipro_tbl, str_detect(InterPro_annotations_accession, pattern="^IPR"))  
length(setdiff(comb_index, ipro_tbl$Protein_accession)) 
length(setdiff(comb_index, ipro_wiproid_tbl$Protein_accession))


refannot_gp <- read_tsv("../CopciAB_new_annot_CDS_20220425m_newIds_isoforms_mancur.genePred", col_names = F, col_types = cols(X9=col_character(), X10=col_character()))
refannot_gr <- import("../CopciAB_new_annot_CDS_20220425m_newIds_isoforms_mancur.gtf")


out_dir <- "mmseqs_out_20220425m_RNAseq_filtered"
dir.create(out_dir)

mmseqs_cdq_cmd <- str_c("mmseqs createdb ", "CopciAB_new_annot_CDS_20220425m_RNAseq_filtered.fasta ", out_dir, "/CopciAB_new_annot_CDS_20220425m_RNAseq_filtered_db")
system(command = mmseqs_cdq_cmd)

mmseqs_search_cmd <- str_c("mmseqs search ", "mmseqs_out_20220425m_RNAseq_filtered/CopciAB_new_annot_CDS_20220425m_RNAseq_filtered_db ", "/SSD5/Uniprot_v2022_01/uniref100_V2022_01_tax.db ", "lowRNAseq_vs_uniprot ", "/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/clean_annotation_20220427_dir/filter_by_RNAseq_readsupport/mmseqs_out_20220425m_RNAseq_filtered/temp ", "--threads 50")
system(command = mmseqs_search_cmd)

mmseqs_convert_cmd <- str_c("mmseqs convertalis --format-output query,target,pident,evalue,bits,qstart,qend,tstart,tend,qlen,tlen,qcov ", 
                            "mmseqs_out_20220425m_RNAseq_filtered/CopciAB_new_annot_CDS_20220425m_RNAseq_filtered_db ",
                            "/SSD5/Uniprot_v2022_01/uniref100_V2022_01_tax.db ",
                            "lowRNAseq_vs_uniprot ",
                            "lowRNAseq_vs_uniprot.tsv")
system(command = mmseqs_convert_cmd)


mmseqs_name <- unlist(str_split("query,target,pident,evalue,bits,qstart,qend,tstart,tend,qlen,tlen,qcov", pattern = ","))
mmseqs_uni <- read_tsv("lowRNAseq_vs_uniprot.tsv", col_names = mmseqs_name)
mmseqs_uni <- mutate(mmseqs_uni, TaxID = str_split(target, pattern = "\\|", simplify=T)[,2])
mmseqs_uni <- mutate(mmseqs_uni, TaxID = str_remove(TaxID, pattern = "TaxID:"))

sum(mmseqs_uni$TaxID %in% c("240176", "1132390"))

length(sprot_seq) 
length(unique(mmseqs_uni$query))

mmseqs_uni_f <- filter(mmseqs_uni, evalue <= 1e-10, qcov>=0.5, !TaxID %in% c("240176", "1132390"))
length(unique(mmseqs_uni_f$query))

mmseqs_uni_f_lt <- split(mmseqs_uni_f, mmseqs_uni_f$query)

NCBI_taxid <- readLines("NCBI_taxonomy/rankedlineage.dmp")
NCBI_taxid_lt <- lapply(NCBI_taxid, function(x) {
  return(str_split(x, pattern = "\t\\|\t|\t\\|$", simplify = T))
})
Copci_index <- str_detect(sapply(NCBI_taxid_lt, "[[", 2), pattern = "^Coprinopsis")
NCBI_taxid_Copci_lt <- NCBI_taxid_lt[Copci_index]
NCBI_taxid_Copci <- purrr::reduce(NCBI_taxid_Copci_lt, rbind) 
Copcis_Taxid <- NCBI_taxid_Copci[,1] 

sum(mmseqs_uni$TaxID %in% c("240176", "1132390"))
sum(mmseqs_uni$TaxID %in% Copcis_Taxid)

mmseqs_uni_f2 <- filter(mmseqs_uni, evalue <= 1e-10, qcov>=0.5, !TaxID %in% Copcis_Taxid)
length(unique(mmseqs_uni_f2$query))

mmseqs_uni_f2_lt <- split(mmseqs_uni_f2, mmseqs_uni_f2$query)

longread_support_tbl <- read_tsv("../splice_perfectmatch_monoexon_match_reads_20220425m_summary.tsv")


length(intersect(comb_index, longread_support_tbl$transcript_id))/length(comb_index) 
findex <- filter(longread_support_tbl,transcript_id %in% comb_index)$transcript_id 
filter(refannot_gp,X1 %in% findex)
table(filter(refannot_gp,X1 %in% findex)$X8)


length(intersect(names(mmseqs_uni_f2_lt), longread_support_tbl$transcript_id))/length(names(mmseqs_uni_f2_lt))


comb_index
longreadsupport_filter_tbl <- filter(longread_support_tbl,transcript_id %in% comb_index) 
longreadsupport_filter_tbl <- filter(longreadsupport_filter_tbl, read_n>=3)$transcript_id
uniprot_sim <- unique(mmseqs_uni_f2$query) 
interpro_acc <- intersect(comb_index, ipro_wiproid_tbl$Protein_accession)
stayin_index <- union(union(longreadsupport_filter_tbl, uniprot_sim), interpro_acc) 
comb_index2 <- setdiff(comb_index, stayin_index) 

upset(fromList(list(longread=longreadsupport_filter_tbl, uniprot=uniprot_sim, interpro=interpro_acc)), order.by = "freq")


sum(comb_index2 %in% ipro_tbl$Protein_accession) 
filter(ipro_tbl, Protein_accession %in% comb_index2)
table(filter(ipro_tbl, Protein_accession %in% comb_index2)$Signature_description) %>% knitr::kable()
filter(ipro_tbl, Protein_accession %in% comb_index2, 
       !Signature_description %in% c("consensus disorder prediction", "Coil"))


comb_index2_pat <- str_remove(comb_index2, pattern="[0-9]+$")
length(refannot_gp[str_remove(refannot_gp$X1, pattern="[0-9]+$") %in% comb_index2_pat,]$X1) == length(comb_index2)
refannot_filter_gp <- refannot_gp[!str_remove(refannot_gp$X1, pattern="[0-9]+$") %in% comb_index2_pat,]
nrow(refannot_filter_gp)
write_tsv(refannot_filter_gp, "../CopciAB_new_annot_CDS_20220425m_newIds_isoforms_mancur_filtered.genePred", col_names = F) 


refannot_gr
length(unique(refannot_gr[str_remove(refannot_gr$gene_id, pattern="[0-9]+$") %in% comb_index2_pat,]$gene_id)) == length(comb_index2)
refannot_filter_gr <- refannot_gr[!str_remove(refannot_gr$gene_id, pattern="[0-9]+$") %in% comb_index2_pat,]
length(unique(refannot_filter_gr$gene_id))
export(refannot_filter_gr, "../CopciAB_new_annot_CDS_20220425m_newIds_isoforms_mancur_filtered.gtf")


prot_filter_seq <- prot_seq[!str_remove(names(prot_seq), pattern="[0-9]+$") %in% comb_index2_pat,]
prot_filter_T0_seq <- prot_filter_seq[str_detect(names(prot_filter_seq), pattern = "T0$")]

writeXStringSet(prot_filter_seq, "../CopciAB_new_annot_CDS_20220425m_newIds_isoforms_mancur_filtered_prot.fasta")
writeXStringSet(prot_filter_T0_seq, "../CopciAB_new_annot_CDS_20220425m_newIds_isoforms_mancur_filtered_T0_prot.fasta")


busco_cmd <-  "run_BUSCO.py -f -i ../CopciAB_new_annot_CDS_20220425m_newIds_isoforms_mancur_filtered_T0_prot.fasta -o BUSCO_CopciAB_new_annot_CDS_20220425m_newIds_isoforms_mancur_filtered_T0_prot_odb10 -l /work/balintb/Busco/Databases/v4/basidiomycota_odb10/ -c 10 -m prot"
system(command = busco_cmd)