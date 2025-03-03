suppressPackageStartupMessages({
library(stringr)
library(parallel)
library(stringr)
library(tibble)
library(readr)
library(dplyr)
})

refgenome_path <- "/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/ref_genome/fullgenome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326.contigs.fsa"
aln_bam_path <- list.files(".", pattern="_subset_s.bam$", full.names=TRUE)


# tama collapse -------------------------------------------------

tama_collapse_runner_fct <- function(x) {
  out_tag <- str_extract(x, pattern = "NB[0-9]+")
  out_dir_path <- str_c(out_tag, "_tama_collapse_out")
  dir.create(out_dir_path)
  out_file_path <- str_c(out_dir_path, "/",out_tag)
  command <- str_c("tama_collapse.py -s ", x, " -b BAM -f ", refgenome_path, " -p ", out_file_path, " -x capped -icm ident_map -a 50 -z 50")
  system(command = command)
}

mclapply(aln_bam_path, FUN = tama_collapse_runner_fct, mc.cores = 2)


# clean reads header in the trans_read.bed file ---------------------------
wd <- "/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/tama_out"
setwd(wd)

trans_read_bed_path <- list.files("./tama_collapse_out/", patter="_trans_read.bed", full.names = TRUE)

read_header_cleaner_fct <- function(x) {
  out_tag <- str_extract(x, pattern = "NB[0-9]+_trans_read")
  trans_read_bed <- read_tsv(x, col_names=FALSE)
  trans_read_bed$X4 <- str_remove_all(trans_read_bed$X4, pattern="[0-9]+:[0-9]+\\|")
  write_tsv(trans_read_bed, str_c("./tama_collapse_out/", out_tag, "_clean.bed"), col_names=FALSE)
}

mclapply(trans_read_bed_path, FUN = read_header_cleaner_fct, mc.cores = length(trans_read_bed_path))

gc()


# tama read support -------------------------------------------------------

if(!dir.exists("./tama_collapse_read_support_out")) {
  dir.create("./tama_collapse_read_support_out")
}

# read cleaned trans_reads path 

trans_read_bed_path <- list.files("./tama_collapse_out/", pattern = "_trans_read_clean.bed", full.names=TRUE)

# create filelist.txt
filelist_creater_fct <- function(x) {
  out_tag <- str_extract(x, pattern = "NB[0-9]+")
  filelist_tbl <- tibble(X1 = out_tag, X2 = x, X3 = "trans_read")
  write_tsv(filelist_tbl, path = str_c("./tama_collapse_read_support_out/", out_tag, "_collapse_readsupport_filelist.txt"), col_names = FALSE)
}

lapply(trans_read_bed_path, FUN = filelist_creater_fct)

tama_read_sup_path <- "python /SSD/software/tama/tama_go/read_support/tama_read_support_levels.py"

filelist_path <- list.files("./tama_collapse_read_support_out/", pattern = "collapse_readsupport_filelist.txt", full.names = TRUE)

tama_read_sup_caller_fct <- function(x) {
  out_tag <- str_extract(x, pattern = "NB[0-9]+")
  command <- str_c(tama_read_sup_path, " -f ", x, " -m no_merge ", "-o ./tama_collapse_read_support_out/", out_tag, "_collapse")
  # print(command)
  system(command = command)
}

mclapply(filelist_path, FUN = tama_read_sup_caller_fct, mc.cores = length(filelist_path))


# tama filtering ----------------------------------------------------------

###### filter out single read models

if(!dir.exists("./tama_s_filtering_out")) {
  dir.create("./tama_s_filtering_out")
}

# annotation_bed
annot_bed_path <- list.files("./tama_collapse_out/", patter="NB[0-9]+.bed", full.names = TRUE)
# readsupport
readsupport_path <- list.files("./tama_collapse_read_support_out/", pattern = "read_support.txt$", full.names = TRUE)

tama_single_read_models_filter_path <- "python /SSD/software/tama/tama_go/filter_transcript_models/tama_remove_single_read_models_levels.py"

copci_samples <- str_extract(annot_bed_path, pattern = "NB[0-9]+")
tama_single_rmodel_filter_caller_fct <- function(x) {
  annot_bed <- str_subset(annot_bed_path, pattern = x)
  readsupport <- str_subset(readsupport_path, pattern = x)
  
  command <- str_c(tama_single_read_models_filter_path, " -b ", annot_bed, " -r ", readsupport, " -l transcript -k remove_multi", " -o ./tama_s_filtering_out/", x)
  # print(command)
  system(command = command)
}

mclapply(copci_samples, FUN = tama_single_rmodel_filter_caller_fct, mc.cores = length(copci_samples))

###### read support for filtered samples
trans_read_bed_path <- list.files("./tama_collapse_out/", pattern = "_trans_read_clean.bed", full.names=TRUE)
mergefile_path <- list.files("./tama_s_filtering_out/", pattern = "report", full.names = TRUE)

create_readsupport_filelist_fct <- function(x) {
  out_tag <- str_extract(x, pattern = "NB[0-9]+")
  readsupport_filelist_tbl <- tibble(X1 = out_tag,
                                     X2 = x,
                                     X3 = "trans_read")
  write_tsv(readsupport_filelist_tbl, str_c("./tama_s_filtering_out/" ,out_tag, "_readsupport_filelist.txt"), col_names = FALSE)
}

lapply(trans_read_bed_path, FUN = create_readsupport_filelist_fct)

readsupport_filelist_path <- list.files("./tama_s_filtering_out/", pattern = "readsupport_filelist.txt", full.names = TRUE)

tama_read_sup_path <- "python /SSD/software/tama/tama_go/read_support/tama_read_support_levels.py"

sample_index <- str_extract(readsupport_filelist_path, pattern = "NB[0-9]+")
tama_read_sup_caller2_fct <- function(x) {
  filelist <- str_subset(readsupport_filelist_path, pattern = x)
  mergefile <- str_subset(mergefile_path, pattern = x)
  command <- str_c(tama_read_sup_path, " -f ", filelist, " -m ", mergefile, " -o ./tama_s_filtering_out/", x, "_s_filter ", "-mt filter")
  system(command = command)
}

mclapply(sample_index, FUN = tama_read_sup_caller2_fct, mc.cores = length(sample_index))


# merge filtered filtered samples -----------------------------------------

if(!dir.exists("./tama_merge_out")) {
  dir.create("./tama_merge_out")
}


filtered_annot_bed_path <- list.files("./tama_s_filtering_out/", pattern = "NB[0-9]+.bed", full.name = TRUE)

# create file_list

annot_filelist_tbl <- tibble(X1 = filtered_annot_bed_path, 
                             X2 = "capped",
                             X3 = "1,1,1",
                             X4 = str_extract(filtered_annot_bed_path, pattern = "NB[0-9]+"))

annot_filelist_path <- "./tama_merge_out/annot_filelist.txt"
write_tsv(annot_filelist_tbl, annot_filelist_path, col_names = FALSE)

tama_merge_command <- str_c("tama_merge.py ", "-f ", annot_filelist_path, " -p ./tama_merge_out/Copci_merged_annots ", "-a 50 -z 50 ")
print(tama_merge_command)
system(command = tama_merge_command)

###### read support for merged samples

sfilter_readsupport_path <- list.files("./tama_s_filtering_out/", pattern="_filter_read_support.txt$", full.names = TRUE)

sfilter_readsupport_tbl <- tibble(X1 = str_extract(sfilter_readsupport_path, pattern = "NB[0-9]+"),
                                  X2 = sfilter_readsupport_path,
                                  X3 = "read_support")

sfilter_readsupport_filelist_path <- "./tama_merge_out/sfilter_readsupport_filelist.txt"
write_tsv(sfilter_readsupport_tbl, sfilter_readsupport_filelist_path, col_names = FALSE)

merge_file_path <- "./tama_merge_out/Copci_merged_annots_merge.txt"
tama_read_sup_path <- "python /SSD/software/tama/tama_go/read_support/tama_read_support_levels.py"

readsupport_command <- str_c(tama_read_sup_path, " -f ", sfilter_readsupport_filelist_path, " -o ./tama_merge_out/Copci_merged_annots ", "-m ", merge_file_path)
print(readsupport_command)

system(command = readsupport_command)


# convert merge bed annotation into gtf format ----------------------------

tama_bed_to_gtf_converter_path <- "python /SSD/software/tama/tama_go/format_converter/tama_convert_bed_gtf_ensembl_no_cds.py"
merged_bed_annot_path <- "./tama_merge_out/Copci_merged_annots.bed"
merged_gtf_annot_path <- "./tama_merge_out/Copci_merged_annots.gtf"

bed_to_gtf_command <- str_c(tama_bed_to_gtf_converter_path, " ", merged_bed_annot_path, " ", merged_gtf_annot_path)
print(bed_to_gtf_command)
system(command = bed_to_gtf_command)

#### convert the gtf annotation to into a sqanti2 comapatible format

suppressPackageStartupMessages({library(rtracklayer)})

merge_gtf <- import(merged_gtf_annot_path, format = "gtf")
merge_gtf_df <-  as.data.frame(merge_gtf)

merge_gtf_df %>%
  mutate(gene_id = str_replace(gene_id, pattern = "^G", replacement = "PB."),
         transcript_id = str_replace(transcript_id, pattern = "^G", replacement = "PB."),
         uniq_trans_id = str_replace(uniq_trans_id, pattern = "^G", replacement = "PB.")) %>% 
  filter(type != "gene") %>% 
  select(-uniq_trans_id, -exon_number) -> merge_gtf_mfor_sqanty2_df

export(merge_gtf_mfor_sqanty2_df, "./tama_merge_out/Copci_merged_annots_mfor_sqanti2.gtf", format = "gtf")

# quick summary -----------------------------------------------------


if(!dir.exists("./tama_summary_out")) {
  dir.create("./tama_summary_out")
}

merge_bed <- read_tsv("./tama_merge_out/Copci_merged_annots.bed", col_types = cols(.default = "c"), col_names = FALSE) # minden oszlop karaktert tartalmazzon!
merge_readsupport_tbl <- read_tsv("./tama_merge_out/Copci_merged_annots_read_support.txt", col_types = cols_only(merge_gene_id = "c",
                                                                                                                 merge_trans_id = "c",
                                                                                                                 gene_read_count = "i",
                                                                                                                 trans_read_count = "i",
                                                                                                                 source_line = "c"))

write_tsv(merge_readsupport_tbl, "./tama_merge_out/Copci_merged_annots_read_support_NoSourceLine.txt")

merge_readsupport_tbl %>% 
  group_by(merge_gene_id) %>% 
  filter(trans_read_count == max(trans_read_count)) %>% 
  mutate(bed_index = str_c(merge_gene_id, merge_trans_id, sep = ";"))-> top_gene_tbl


merge_bed_top_gene <- merge_bed[merge_bed$X4 %in% top_gene_tbl$bed_index,]
write_tsv(merge_bed_top_gene, "./tama_summary_out/Copci_merged_annots_top_gene.bed", col_names = FALSE)
