library(readr)
library(stringr)
library(dplyr)
library(knitr)
library(ggplot2)
library(tidyr)
library(tibble)

wd <- "/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/CopciAB_genome/refgenome_QC_dir_220805" # a szerveren

# load data
vcf_files_path <- list.files(".", pattern = "vcf$")

vcf_names <- "CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,S205.1"
vcf_file_lt <- lapply(vcf_files_path, function(x) read_tsv(x, comment = "#", col_names = unlist(str_split(vcf_names, pattern = ","))))
sapply(vcf_file_lt, nrow)

vcf_file <- read_tsv(vcf_files_path[4], comment = "#", col_names = unlist(str_split(vcf_names, pattern = ",")))

vcf_file_copy_lt <- vcf_file_lt
names(vcf_file_copy_lt) <- str_remove(vcf_files_path, pattern="\\.vcf")
as.data.frame(sapply(vcf_file_copy_lt, nrow)) %>% 
  rename('n'='sapply(vcf_file_copy_lt, nrow)') %>% 
  rownames_to_column(var='filename') %>% 
  mutate(filename = ifelse(str_detect(filename, pattern = "indel|SNP"), filename, str_c(filename, "_all"))) %>% 
  mutate(snv_type = sapply(str_split(filename, pattern = '_'), function(x) rev(x)[1])) %>% 
  mutate(strain=rep(c("CopciAB V1", "CopciAB V2", "CopciAB ONT"), each=3)) %>% 
  select(strain, snv_type, n) -> vcf_file_copy_tbl

write_tsv(vcf_file_copy_tbl, 'snv_type_summary_tbl.tsv')