suppressPackageStartupMessages({
  library(stringr)
  library(parallel)
  library(stringr)
  library(tibble)
  library(readr)
  library(dplyr)
  library(rtracklayer)
  library(readxl)
})


# load ref data -------------------------------------------------------------------------------------
ref_annot <- read_tsv("../../CopciAB_new_annot_CDS_20220209_rev3.1_newIds.genePred", col_names = F, col_types = cols(X9= col_character(), X10=col_character()))
ref_annot_gp_se <- select(ref_annot, 1:5)
colnames(ref_annot_gp_se) <- str_c("ref_", colnames(ref_annot_gp_se))

ref_annot_df <- as.data.frame(import("../../CopciAB_new_annot_CDS_20220209_rev3.1_newIds.gtf"))
ref_annot_df %>% filter(type == "CDS") %>% group_by(gene_id) %>% summarise(ref_prot_length = sum(width)/3) -> ref_prot_tbl

# load isoforms -------------------------------------------------------------------------------------
isoform_annot_gp <- read_tsv("../tama_collapse/non_sj_correct_tssc_1000_filter10_cds.genePred", col_names = F, col_types = cols(X9= col_character(), X10=col_character()))
isoform_annot_gp_se <- select(isoform_annot_gp, 1:5)
colnames(isoform_annot_gp_se) <- str_c("isoform_", colnames(isoform_annot_gp_se))
isoform_annot_df <- as.data.frame(import("../tama_collapse/non_sj_correct_tssc_1000_filter10_cds.gtf"))

isoform_annot_df %>% filter(type == "CDS") %>% group_by(gene_id) %>% summarise(iso_prot_length = sum(width)/3) -> iso_prot_tbl

# count adatok betöltése ---------------------------------------------------------------------------------
isoform_count <- read_tsv("../tama_collapse_read_support_out/non_sj_correct_tssc_1000_collapse_read_support_filter10.txt")
pmatch_count <- read_tsv("../../splice_perfectmatch_site_match_reads_20220224_summary.tsv")
mmatch_count <- read_tsv("../../monoexon_match_reads_20220224_summary.tsv")
mmatch_count <- dplyr::rename(mmatch_count, transcript_id=ref_id)

jgi_dict <- read_tsv("../../prot_dict_dir/CopciAB_id_dict_20220224.tsv")

all_count <- bind_rows(pmatch_count, mmatch_count)
all_count <- left_join(all_count, jgi_dict[c("all_isoforms","new_name")], by=c("transcript_id"="all_isoforms"))


iso_class <- read_tsv("../sqanti3_out_1/tama_merge_classification.txt")

length(unique(isoform_annot_df$gene_id)) 
length(unique(iso_class$isoform)) 

iso_class <- filter(iso_class, isoform %in% unique(isoform_annot_df$gene_id))




iso_class_short <- select(iso_class, isoform, associated_gene)

filter(iso_class_short, isoform=="G13009.23") 

sum(str_detect(iso_class_short$associated_gene, pattern="novelGene"), na.rm=T) 
novelGene <- iso_class_short[str_detect(iso_class_short$associated_gene, pattern="novelGene") & !is.na(iso_class_short$associated_gene),]
isoform_annot_df[isoform_annot_df$gene_id == "G1053.3",] 
setdiff(novelGene$isoform,isoform_annot_df$gene_id) 
sum(novelGene$isoform %in% isoform_annot_df$gene_id)
novelGene_annot <- isoform_annot_df[isoform_annot_df$gene_id %in% novelGene$isoform,]
export(novelGene_annot, "novelGene_annot.gtf")
novelGene_annot_tscript <- filter(novelGene_annot, type == "transcript") %>% select(-exon_number, -exon_id, -source, -type, -score, -phase) %>% mutate(pas_index = str_c(seqnames, ":", start, "-", end))
novelGene_annot_tscript <- arrange(novelGene_annot_tscript, nchar(as.character(seqnames)), as.character(seqnames), start)
write_tsv(novelGene_annot_tscript, "novelGene_annot_simple.tsv")

iso_class_short_old <- iso_class_short 
iso_class_short <- iso_class_short[!is.na(iso_class_short$associated_gene),]
iso_class_short <- filter(iso_class_short, !isoform %in% novelGene_annot_tscript$gene_id)
sum(str_detect(iso_class_short$associated_gene, pattern="novel"))

iso_class_short_lt <- split(iso_class_short, iso_class_short$isoform)

iso_class_short_c1_lt <- lapply(iso_class_short_lt, function(x) {
  return(unlist(str_extract_all(x$associated_gene, pattern="CopciAB_[0-9]+")))
})

table(sapply(iso_class_short_c1_lt, length))

iso_class_short_c3_lt <- lapply(iso_class_short_c1_lt, function(x) {
  result <- data.frame(matrix(x, ncol = 1))
  colnames(result) <- "id"
  return(result)
})

length(iso_class_short_c3_lt[sapply(iso_class_short_c3_lt, nrow) > 1]) 
iso_class_short_c3_tbl <- bind_rows(iso_class_short_c3_lt, .id="isoform")

ref_annot_exon_df <- filter(ref_annot_df, type == "exon")
ref_annot_exon_lt <- split(ref_annot_exon_df, ref_annot_exon_df$gene_id)

get_splice_chain_fct <- function(x) {
  if(nrow(x)>1) {
    end <- 1:(nrow(x)-1)
    start <- 2:nrow(x)
    result <- tibble(splice_charin = str_c(str_c(x$end[end], x$start[start], sep = "-"), collapse = "_"))
    return(result)
  }
}

ref_splice_sites_lt <- lapply(ref_annot_exon_lt, get_splice_chain_fct)
ref_splice_sites_tbl <- bind_rows(ref_splice_sites_lt, .id = "gene_id")

isoform_annot_exon_df <- filter(isoform_annot_df, type=="exon")
isoform_annot_exon_lt <- split(isoform_annot_exon_df, isoform_annot_exon_df$gene_id)

isoform_splice_sites_lt <- lapply(isoform_annot_exon_lt, get_splice_chain_fct)
isoform_splice_sites_tbl <- bind_rows(isoform_splice_sites_lt, .id = "gene_id")

iso_class_short_c3_tbl_splicechain <- left_join(iso_class_short_c3_tbl, ref_splice_sites_tbl, by=c("id"="gene_id"))
iso_class_short_c3_tbl_splicechain <- dplyr::rename(iso_class_short_c3_tbl_splicechain, splice_charin_ref = splice_charin)
iso_class_short_c3_tbl_splicechain <- left_join(iso_class_short_c3_tbl_splicechain, isoform_splice_sites_tbl, by=c("isoform"="gene_id"))
iso_class_short_c3_tbl_splicechain <- dplyr::rename(iso_class_short_c3_tbl_splicechain, splice_charin_isoform = splice_charin)

iso_class_short_c3_tbl_splicechain_lt <- split(iso_class_short_c3_tbl_splicechain, iso_class_short_c3_tbl_splicechain$id)
splicechain_match_index <- sapply(iso_class_short_c3_tbl_splicechain_lt, function(x) {
  x <- iso_class_short_c3_tbl_splicechain_lt$CopciAB_161813
  return(any(str_detect(x$splice_charin_ref, pattern = x$splice_charin_isoform)))
})

table(splicechain_match_index, useNA="ifany") 
iso_class_short_c3_tbl$nongene_index <- !iso_class_short_c3_tbl$id %in% ref_annot$X1
iso_class_short_c3_tbl[iso_class_short_c3_tbl$nongene_index,]

iso_class_m <- left_join(iso_class, iso_class_short_c3_tbl, by="isoform") 
table(iso_class_m$nongene_index, useNA="ifany") 
dim(iso_class_short_c3_tbl)

iso_class_m <- filter(iso_class_m, !is.na(nongene_index))


length(unique(filter(iso_class_m, structural_category == "fusion")$isoform)) 

iso_class_fusio <- dplyr::select(iso_class_m, isoform, id, structural_category, subcategory)
iso_class_fusio_lt <- split(iso_class_fusio, iso_class_fusio$isoform)

sum(sapply(iso_class_fusio_lt, nrow) > 1)
multyrow_lt <- iso_class_fusio_lt[sapply(iso_class_fusio_lt, nrow) > 1]
sum(sapply(multyrow_lt, function(x) any(str_detect(x$structural_category, "fusion")))) 

multyrow_not_fusion_lt <- multyrow_lt[sapply(multyrow_lt, function(x) !all(str_detect(x$structural_category, "fusion")))]
length(multyrow_not_fusion_lt) 

table(sapply(multyrow_not_fusion_lt, function(x) unique(x$structural_category)))
table(sapply(multyrow_not_fusion_lt, function(x) unique(x$subcategory)))


multyrow_not_fusion_tbl <- bind_rows(multyrow_not_fusion_lt)
filter(ref_annot_df, gene_id=="CopciAB_501015")

multigene_no_fannot <- isoform_annot_df[isoform_annot_df$gene_id %in% multyrow_not_fusion_tbl$isoform,]
export(multigene_no_fannot, "multigene_no_fannot.gtf")
multigene_no_fannot_tscript <- filter(multigene_no_fannot, type == "transcript") %>% select(-exon_number, -exon_id, -source, -type, -score, -phase) %>% mutate(pas_index = str_c(seqnames, ":", start, "-", end)) %>% select(gene_id, pas_index)
multyrow_not_fusion_tbl <- left_join(multyrow_not_fusion_tbl, multigene_no_fannot_tscript, by=c("isoform"="gene_id"))
sum(is.na(multyrow_not_fusion_tbl$pas_index)) 
write_tsv(multyrow_not_fusion_tbl, "multigene_no_fannot.tsv")


multyrow_not_fusion_man_tbl <- read_excel("multigene_no_fannot.xlsx", excel_sheets("multigene_no_fannot.xlsx")[[1]])
multyrow_not_fusion_man_s_tbl <- select(multyrow_not_fusion_man_tbl, isoform, id, ok_index)


iso_class_fusio_f_lt <- iso_class_fusio_lt[sapply(iso_class_fusio_lt, function(x) any(str_detect(x$structural_category, "fusion")))]
length(iso_class_fusio_f_lt) 
iso_class_fusio_f_tbl <- bind_rows(iso_class_fusio_f_lt)

multigene_fannot <- isoform_annot_df[isoform_annot_df$gene_id %in% iso_class_fusio_f_tbl$isoform,]
export(multigene_fannot, "multigene_fannot.gtf")
multigene_fannot_tscript <- filter(multigene_fannot, type == "transcript") %>% select(-exon_number, -exon_id, -source, -type, -score, -phase) %>% mutate(pas_index = str_c(seqnames, ":", start, "-", end)) %>% select(gene_id, pas_index)
multyrow_fusion_tbl <- left_join(iso_class_fusio_f_tbl, multigene_fannot_tscript, by=c("isoform"="gene_id"))
sum(is.na(multyrow_fusion_tbl$pas_index)) 
write_tsv(multyrow_fusion_tbl, "multigene_fannot.tsv")

"multigene_fannot.xlsx" 
multyrow_fusion_man_tbl <- read_excel("multigene_fannot.xlsx", excel_sheets("multigene_fannot.xlsx")[[1]])
multyrow_fusion_man_s_tbl <- select(multyrow_fusion_man_tbl, isoform, id, ok_index) 


multyrow_all_man_s_tbl <- bind_rows(multyrow_not_fusion_man_s_tbl, multyrow_fusion_man_s_tbl)
table(multyrow_all_man_s_tbl$ok_index, useNA="ifany")
multyrow_forfiltering <- filter(multyrow_all_man_s_tbl, ok_index=="f")

iso_class_m2 <- iso_class_m[!(iso_class_m$isoform %in% multyrow_forfiltering$isoform & iso_class_m$id %in% multyrow_forfiltering$id),]
any(duplicated(iso_class_m2$isoform))


isoform_dict <- select(iso_class_m2, isoform, id)

isoform_ref_count <- left_join(isoform_dict, all_count, by=c("id"="new_name"))
isoform_ref_count <- left_join(isoform_ref_count, isoform_count[c("merge_trans_id", "gene_read_count","trans_read_count")], by=c("isoform"="merge_trans_id"))

apply(isoform_ref_count,2, function(x) sum(is.na(x)))
filter(isoform_ref_count, is.na(transcript_id))
length(unique(filter(isoform_ref_count, is.na(transcript_id))$id)) 
filter(ref_annot_df, gene_id == unique(filter(isoform_ref_count, is.na(transcript_id))$id)[[7]])

write_tsv(isoform_ref_count, "ref_iso_counts.tsv")

filter(isoform_ref_count, transcript_id=="G5678.3_t0_r1") 
filter(isoform_ref_count, id=="CopciAB_495339")

non_important_cs <- "gene_read_count;associated_transcript;diff_to_TSS;diff_to_TTS;diff_to_gene_TSS;diff_to_gene_TTS;RTS_stage;min_sample_cov;FL;n_indels;n_indels_junc;bite;iso_exp;gene_exp;ratio_exp;FSM_class;CDS_length;CDS_start;CDS_end;CDS_genomic_start;CDS_genomic_end;predicted_NMD;perc_A_downstream_TTS;seq_A_downstream_TTS;dist_to_cage_peak;within_cage_peak;dist_to_polya_site;within_polya_site;polyA_motif;polyA_dist"
non_important_cs <- unlist(str_split(non_important_cs, pattern = ";"))

iso_class_m_f <- iso_class_m2
iso_class_m_f <- iso_class_m_f[!colnames(iso_class_m_f) %in% non_important_cs]


count(iso_class_m_f, structural_category, subcategory) 

iso_class_m_f2 <- filter(iso_class_m_f, 
                         !structural_category %in% c("incomplete-splice_match", "antisense", "intergenic"),
                         min_cov >= 10 | is.na(min_cov),
                         all_canonical == "canonical"| is.na(all_canonical)) 
count(iso_class_m_f2, structural_category, subcategory)


not_in_filtered_class <- na.omit(isoform_ref_count$isoform)[!na.omit(isoform_ref_count$isoform) %in% iso_class_m_f2$isoform] 
isoform_ref_count_f <- isoform_ref_count[!isoform_ref_count$isoform %in% not_in_filtered_class,]

dim(isoform_ref_count_f)
dim(iso_class_m_f2)

iso_class_m_f2_c <- full_join(iso_class_m_f2, isoform_ref_count_f, by=c("isoform","id"))
dim(iso_class_m_f2_c)
length(unique(iso_class_m_f2_c[is.na(iso_class_m_f2_c$read_n),]$id)) 
iso_class_m_f2_c_f1 <- filter(iso_class_m_f2_c, !is.na(read_n))

iso_class_m_f2_c_f1_lt <- split(iso_class_m_f2_c_f1, iso_class_m_f2_c_f1$id) 

genes_with_isoform <- iso_class_m_f2_c_f1_lt

sum(sapply(genes_with_isoform, function(x) any(is.na(x$gene_read_count))))
sum(sapply(genes_with_isoform, function(x) any(is.na(x$read_n))))

genes_with_isoform_readp <- lapply(genes_with_isoform, function(x) {
  all_read_n <- sum(x$gene_read_count[[1]], x$read_n[[1]])
  x$read_p <- x$trans_read_count / all_read_n
  return(x)
})


genes_with_isoform_readp_10 <- lapply(genes_with_isoform_readp, function(x) {
  return(filter(x, read_p >= 0.1 & trans_read_count >= 10))
})

genes_with_isoform_readp_10 <- genes_with_isoform_readp_10[sapply(genes_with_isoform_readp_10, nrow) > 0]
length(genes_with_isoform_readp_10) 
genes_with_isoform_readp_10_tbl <- bind_rows(genes_with_isoform_readp_10)

count(genes_with_isoform_readp_10_tbl, structural_category, subcategory)
filter(genes_with_isoform_readp_10_tbl, subcategory == "at_least_one_novel_splicesite")

isoform_annot_readp_10_df <- filter(isoform_annot_df, gene_id %in% genes_with_isoform_readp_10_tbl$isoform)
length(unique(isoform_annot_readp_10_df$gene_id))

export(isoform_annot_readp_10_df, "read_supp_filter_isoform.gtf")

isoform_annot_readp_10_refid_df <- left_join(isoform_annot_readp_10_df, iso_class_m_f2_c_f1[c("isoform", "id")], by=c("gene_id"="isoform"))
any(duplicated(iso_class_m_f2_c_f1$isoform))

export(isoform_annot_readp_10_refid_df, "read_supp_filter_isoform_refid.gtf")

isoform_man_check_helper_tbl <- filter(isoform_annot_readp_10_refid_df, type=="transcript")

length(unique(isoform_man_check_helper_tbl$gene_id)) 
length(unique(isoform_man_check_helper_tbl$id)) 

isoform_man_check_helper_tbl <- isoform_man_check_helper_tbl %>% arrange(nchar(as.character(seqnames)), as.character(seqnames), start) %>% mutate(loc_index = str_c(seqnames, ":", start, "-", end)) %>% select(gene_id, id, loc_index)
write_tsv(isoform_man_check_helper_tbl, "read_supp_filter_isoform_refid_man_helper.tsv")

isoform_man_check_helper_man_tbl <- read_excel("read_supp_filter_isoform_refid_man_helper.xlsx", excel_sheets("read_supp_filter_isoform_refid_man_helper.xlsx")[[1]])
table(isoform_man_check_helper_man_tbl$ok_index, useNA="ifany")
filter(isoform_man_check_helper_man_tbl, ok_index=="ok") %>% summarise(unique_gene_n = length(unique(gene_id))) 

isoform_man_check_helper_man_okonly_tbl <- filter(isoform_man_check_helper_man_tbl, ok_index=="ok")
isoform_man_check_helper_man_okonly_tbl <- select(isoform_man_check_helper_man_okonly_tbl, gene_id, id)

setdiff(isoform_man_check_helper_man_okonly_tbl$gene_id, genes_with_isoform_readp_10_tbl$isoform)
genes_with_isoform_readp_10_tbl_manfilter <- filter(genes_with_isoform_readp_10_tbl, isoform %in% isoform_man_check_helper_man_okonly_tbl$gene_id) %>% select(isoform, id, read_p)

genes_with_isoform_readp_10_tbl_manfilter_lt <- split(genes_with_isoform_readp_10_tbl_manfilter, genes_with_isoform_readp_10_tbl_manfilter$id)

genes_with_isoform_readp_10_tbl_manfilter_lt <- lapply(genes_with_isoform_readp_10_tbl_manfilter_lt, function(x) {
  x$isoform_index <- str_c("T", seq_len(length(x$read_p)))[order(x$read_p, decreasing=T)]
  return(x)
})

genes_with_isoform_readp_10_tbl_manfilter_tbl <- bind_rows(genes_with_isoform_readp_10_tbl_manfilter_lt)
genes_with_isoform_readp_10_tbl_manfilter_tbl <- mutate(genes_with_isoform_readp_10_tbl_manfilter_tbl, id = str_c(id, isoform_index, sep = "."))


novelGene_path <- "novelGene_annot_simple.xlsx"
novelGene_man_tbl <- read_excel(novelGene_path, excel_sheets(novelGene_path)[[1]])


table(is.na(novelGene_man_tbl$remove)) 
novelGene_man_to_remove <- filter(novelGene_man_tbl, !is.na(remove))$remove
novelGene_man_to_remove <- str_remove(novelGene_man_to_remove, pattern="\\.[0-9]+$")


table(novelGene_man_tbl$ok_index, useNA="ifany")
novelGene_man_ok_genes <-  filter(novelGene_man_tbl, !is.na(ok_index)) %>% filter(!ok_index %in% c("f", "r")) %>% .$gene_id


novelGene_man_ok_genes_count <- filter(isoform_count[c("merge_trans_id", "gene_read_count","trans_read_count")], merge_trans_id %in% novelGene_man_ok_genes) # izoforma count
novelGene_man_ok_genes_count <- mutate(novelGene_man_ok_genes_count, gene_index = str_remove(merge_trans_id, pattern = "\\.[0-9]+"))
novelGene_man_ok_genes_count_lt <- split(novelGene_man_ok_genes_count, novelGene_man_ok_genes_count$gene_index)


novelGene_man_ok_genes_count_lt <- lapply(novelGene_man_ok_genes_count_lt, function(x) {
  x$count_p <- (x$trans_read_count / x$gene_read_count[[1]])
  return(x)
})


novelGene_man_ok_genes_count_lt <- lapply(novelGene_man_ok_genes_count_lt, function(x) {

  if(nrow(x) > 1) {
    if(any(x$count_p >= 0.1)) {
      x$filter_index <- x$count_p >= 0.1
    } else {
      x$filter_index <- x$count_p == max(x$count_p)
    }
    return(x)
  } else {
    x$filter_index <- TRUE
    return(x)
  }
})


novelGene_man_ok_genes_count_isoform_lt <- lapply(novelGene_man_ok_genes_count_lt, function(x) {
  x$isoform_index <- str_c("T", seq_len(length(x$trans_read_count))-1)[order(x$trans_read_count, decreasing=T)]
  return(x)
})

novelGene_man_ok_genes_count_isoform_tbl <- bind_rows(novelGene_man_ok_genes_count_isoform_lt)
novelGene_man_ok_genes_count_isoform_f_tbl <- filter(novelGene_man_ok_genes_count_isoform_tbl, filter_index) 
ref_annot 
isoform_annot_gp 

prot_dict_tbl <- read_tsv("../../prot_dict_dir/CopciAB_id_dict_20220224.tsv")
novelGene_man_to_remove <- filter(prot_dict_tbl, all_isoforms %in% novelGene_man_to_remove)$new_name 
novelGene_man_to_remove
novelGene_man_ok_genes_count_isoform_f_tbl

ref_annot_nGfilter <- filter(ref_annot, !X1 %in% novelGene_man_to_remove)
ref_annot_nGfilter <- mutate(ref_annot_nGfilter, X1 = str_c(X1, ".T0")) 

setdiff(novelGene_man_ok_genes_count_isoform_f_tbl$merge_trans_id, isoform_annot_gp$X1)
isoform_annot_gp_nGfilter <- filter(isoform_annot_gp, X1 %in% novelGene_man_ok_genes_count_isoform_f_tbl$merge_trans_id)
tindex <- novelGene_man_ok_genes_count_isoform_f_tbl$isoform_index[match(isoform_annot_gp_nGfilter$X1, novelGene_man_ok_genes_count_isoform_f_tbl$merge_trans_id)]
isoform_annot_gp_nGfilter$X1 <- str_c(isoform_annot_gp_nGfilter$X1, tindex, sep = ".")

new_refannot <- bind_rows(ref_annot_nGfilter, isoform_annot_gp_nGfilter)

setdiff(genes_with_isoform_readp_10_tbl_manfilter_tbl$isoform, isoform_annot_gp$X1)
isoform_annot_good_iso_gp <- filter(isoform_annot_gp, X1 %in% genes_with_isoform_readp_10_tbl_manfilter_tbl$isoform) 

isoform_annot_good_iso_gp$X1 <- genes_with_isoform_readp_10_tbl_manfilter_tbl$id[match(isoform_annot_good_iso_gp$X1, genes_with_isoform_readp_10_tbl_manfilter_tbl$isoform)]


new_refannot_withiso <- bind_rows(new_refannot, isoform_annot_good_iso_gp) 

write_tsv(new_refannot_withiso, "CopciAB_new_annot_CDS_20220419_newIds_isoforms.genePred", col_names = F)

togtf_cmd <- str_c('/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/tama_out_refannot/sep_joint_mapping_refannot_tama/genePredToGtf "file" -utr ',
                   " ", "CopciAB_new_annot_CDS_20220419_newIds_isoforms.genePred",
                   " ", str_replace("CopciAB_new_annot_CDS_20220419_newIds_isoforms.genePred", pattern = "genePred", "gtf"))
system(command = togtf_cmd)


                
iso_class_m_f2 <- filter(iso_class_m_f, 
                         !structural_category %in% c("incomplete-splice_match", "antisense", "intergenic"),
                         coding == "coding",
                         min_cov != 0 | is.na(min_cov)) 

iso_class_m_f2 <- left_join(iso_class_m_f2, ref_prot_tbl, by=c("id"="gene_id"))
iso_class_m_f2 <- left_join(iso_class_m_f2, iso_prot_tbl, by=c("isoform"="gene_id"))

iso_class_m_f2_lt <- split(iso_class_m_f2, iso_class_m_f2$id)

count(iso_class_m_f2, structural_category, subcategory)


sum(sapply(iso_class_m_f2_lt, function(x) {any(x$iso_prot_length > unique(x$ref_prot_length))}))

sum(sapply(iso_class_m_f2_lt, function(x) any(x$trans_read_count >= unique(x$read_n)) & 
             any(x$iso_prot_length > unique(x$ref_prot_length))), na.rm=T)

with_problems_lt <- iso_class_m_f2_lt[sapply(iso_class_m_f2_lt, function(x) any(x$trans_read_count >= unique(x$read_n)) & any(x$iso_prot_length > unique(x$ref_prot_length)))]
with_problems_lt <- with_problems_lt[!sapply(with_problems_lt, is.null)]
with_problems_f_lt <- lapply(with_problems_lt, function(x) x[which.max(x$trans_read_count),])
with_problems_f_tbl <- bind_rows(with_problems_f_lt)
with_problems_f_tbl <- left_join(with_problems_f_tbl, ref_annot[c(1,4)], by=c("id"="X1"))

write_tsv(with_problems_f_tbl, "iso_problems.tsv")


id_binder <- select(iso_class_m_f, isoform, id)
ref_iso_coord <- left_join(id_binder,isoform_annot_gp_se, by=c("isoform"="isoform_X1"))
ref_iso_coord <- left_join(ref_iso_coord, ref_annot_gp_se, by=c("id"="ref_X1"))
ref_iso_coord <- ref_iso_coord[!is.na(ref_iso_coord$isoform_X4),]

ref_iso_coord_diff <- mutate(ref_iso_coord, diff_X4 = ref_X4-isoform_X4, diff_X5 = ref_X5-isoform_X5) 


hist(ref_iso_coord_diff$diff_X4, main="", xlab="távolság a TSS-tól", breaks = 1000)
abline(v=c(-100, 100), col="red")
summary(ref_iso_coord_diff$diff_X4)
boxplot(ref_iso_coord_diff$diff_X4, ylab="távolság a TSS-tól")

hist(ref_iso_coord_diff$diff_X5, main="", xlab="távolság a TTS-tól", breaks = 1000)
abline(v=c(-100, 100), col="red")
summary(ref_iso_coord_diff$diff_X5)
boxplot(ref_iso_coord_diff$diff_X5)

ref_iso_coord_diff_index <- ref_iso_coord_diff$diff_X4 >= -100 & ref_iso_coord_diff$diff_X4 <= 100 & ref_iso_coord_diff$diff_X5 >= -100 & ref_iso_coord_diff$diff_X5 <= 100 
ref_iso_coord_diff_filter <- ref_iso_coord_diff[ref_iso_coord_diff_index,]

isoform_annot_se_df <- filter(isoform_annot_df, gene_id %in% ref_iso_coord_diff_filter$isoform)
export(isoform_annot_se_df, "isoform_filter_se.gtf")

iso_class_m_diff_f <- filter(iso_class_m_f, isoform %in% ref_iso_coord_diff_filter$isoform)
count(iso_class_m_diff_f, structural_category, subcategory)

filter(iso_class_m_diff_f, structural_category == "antisense")
filter(iso_class_m_diff_f, structural_category == "incomplete-splice_match")

length(unique(iso_class_m_diff_f$id)) 

wd <- "/data/bhegedus/Projects/pj8_3genseq/Transcript_sequencing/Coprinopsis_cinerea/results/remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/clean_annotation_20220222_dir/tama_out_220122/filtered_isoforms"

new_annot_CDS_gr <- import("../tama_merge_out/Copci_merged_annots_CDS.gtf")
new_annot_CDS_df <- as.data.frame(new_annot_CDS_gr)

ref_annot_gr <- import("../../CopciAB_new_annot_CDS_20220209_rev2.2.gtf")
ref_annot_df <- as.data.frame(ref_annot_gr)


library(readxl)
esheets <- excel_sheets("iso_problems.xlsx")
man_cur_1 <- read_excel("iso_problems.xlsx", esheets[1])

man_cur_1_onlyok <- filter(man_cur_1, !is.na(ok_index))

old_gene_ids <- str_trim(unlist(str_split(man_cur_1_onlyok$id, pattern = ";")))
old_gene_ids <- unique(old_gene_ids)
any(!old_gene_ids %in% ref_annot_df$gene_id)

new_gene_ids <- unique(man_cur_1_onlyok$isoform)
any(!new_gene_ids %in% new_annot_CDS_df$gene_id)

ref_annot_df_filtered <- ref_annot_df[!ref_annot_df$gene_id %in% old_gene_ids,]
new_annot_CDS_df_filtered <- new_annot_CDS_df[new_annot_CDS_df$gene_id %in% new_gene_ids,]

new_annot_CDS_df_filtered %>% mutate(gene_id = str_c(gene_id, "_rev4"),
                                     transcript_id = str_c(transcript_id, "_rev4"),
                                     exon_id = ifelse(!is.na(exon_number), str_c(transcript_id, ".", exon_number), NA)) -> new_annot_CDS_df_filtered_renamed

ref_annot_df_filtered_rev2 <- bind_rows(ref_annot_df_filtered, new_annot_CDS_df_filtered_renamed)
export(ref_annot_df_filtered_rev2, "../../CopciAB_new_annot_CDS_20220209_rev3.1.gtf")

length(unique(ref_annot_df_filtered_rev2$gene_id))
table(sapply(str_split(unique(ref_annot_df_filtered_rev2$gene_id), pattern = "_|\\.|[0-9]+"), "[[", 1))
table(sapply(str_split(unique(ref_annot_df_filtered_rev2$gene_id), pattern = "_|\\.|[0-9]+"), "[[", 1) %in% c("G", "PB", "man"))


