# load packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(rtracklayer)
  library(seqinr)
  library(Biostrings)
  library(BSgenome)
  library(stringr)
  library(parallel)
  library(GenomicFeatures)
  library(readr)
})


# load data -------------------------------------------

ref_genome <- readDNAStringSet("/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/ref_genome/fullgenome_jgi/merged_genome/CopciAB_new_jgi_20220113.fasta")
ref_genome_tbl <- tibble(names = names(ref_genome),
                         width = width(ref_genome))

# DPAC paramters: --MinScore 75 --MinCount 20 --ReqAs 15
pas_f_df <- as.data.frame(import("All_PACSeq_combined-PASs.f.bedgraph", format = "bed"))
pas_r_df <- as.data.frame(import("All_PACSeq_combined-PASs.r.bedgraph", format = "bed"))

# check pas table
min(as.integer(pas_f_df$name)) 
min(pas_f_df$score) 

# process pas table --------------------------------
pas_all_df <- bind_rows(pas_f_df, pas_r_df)
pas_all_df <- arrange(pas_all_df, nchar(as.character(seqnames)), as.character(seqnames), start)
pas_all_df$pas_index <- str_c("pas_", seq_len(nrow(pas_all_df)))
pas_all_df$NA. <- NULL

## pas genome environment  -------------------------------------------

pas_all_m0_df <- mutate(pas_all_df, new_start = ifelse(strand == "-", start-10, start + 1),
                        new_end = ifelse(strand == "-", start-1, start + 10))
pas_all_m0_df <- dplyr::select(pas_all_m0_df, -end, -start, -width, -score)
pas_all_m0_df <- dplyr::rename(pas_all_m0_df, start = new_start, end = new_end)

pas_all_m0_lt <- split(pas_all_m0_df, as.character(pas_all_m0_df$strand))

pas_all_m0_clean_lt <- lapply(pas_all_m0_lt, function(x) {
    long_tbl <- left_join(x, ref_genome_tbl, by=c("seqnames" = "names"))
    result_index <- !(long_tbl$end > long_tbl$width | long_tbl$start < 1)
    result <- long_tbl[result_index,]
    result$width <- NULL
    return(result)
})

pas_up_seq_lt <- lapply(pas_all_m0_clean_lt, function(x) {
  pas_part_gr <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  pas_up_seq <- BSgenome::getSeq(ref_genome, pas_part_gr)
  names(pas_up_seq) <- pas_part_gr$pas_index
  return(pas_up_seq)
})

sapply(pas_up_seq_lt, function(x) all(width(x) == 10))
dim(pas_all_m0_df)
sum(sapply(pas_up_seq_lt, length))

names(pas_up_seq_lt) <- NULL
pas_up_seq_all <- do.call(c, pas_up_seq_lt)
pas_up_seq_all <- as.data.frame(pas_up_seq_all)
colnames(pas_up_seq_all) <- "up_seq"
pas_up_seq_all <- tibble::rownames_to_column(pas_up_seq_all, var = "pas_index")

pas_all_seq_df <- left_join(pas_all_df, pas_up_seq_all, by="pas_index") 
pas_all_seq_df$polyA_n <- nchar(str_extract(pas_all_seq_df$up_seq, pattern = "^A+"))
pas_all_seq_df$polyA_n[is.na(pas_all_seq_df$polyA_n)] <- 0
pas_all_seq_df$A_p <- sapply(str_extract_all(pas_all_seq_df$up_seq, pattern = "A"), length) / 10


## filter polyA pas --------------------------------------------
filter(pas_all_seq_df, polyA_n >= 4)
filter(pas_all_seq_df, A_p >= 0.7)

pas_all_seq_filtered_df <- pas_all_seq_df[!(pas_all_seq_df$A_p >= 0.7 | pas_all_seq_df$polyA_n >= 4),]
export(pas_all_seq_filtered_df, "All_PASs_rf.sorted_filtered.bed", format = "bed")

summary(pas_all_seq_filtered_df$polyA_n)
summary(pas_all_seq_filtered_df$A_p)

pas_all_seq_filtered_lt <- split(pas_all_seq_filtered_df, as.character(pas_all_seq_filtered_df$strand))

filtere_bed_files <- c("All_PASs_r_filtered.bed", "All_PASs_f_filtered.bed")
export(pas_all_seq_filtered_lt$`-`, filtere_bed_files[[1]])
export(pas_all_seq_filtered_lt$`+`, filtere_bed_files[[2]])

# clusater pas into pac (10 nt) -----------------------------------------------------------------------------

# sort
filtere_sorted_bed_files <- c("All_PASs_r.sorted_filtered.bed", "All_PASs_f.sorted_filtered.bed")
bed_sort_comm <- str_c("sort -k1,1 -k2,2n ", filtere_bed_files[[1]], " > ", filtere_sorted_bed_files[[1]])
system(command = bed_sort_comm)

bed_sort_comm <- str_c("sort -k1,1 -k2,2n ", filtere_bed_files[[2]], " > ", filtere_sorted_bed_files[[2]])
system(command = bed_sort_comm)

# cluster
merge_dist <- 10

filtere_sorted_merged_bed_files <- c("All_PASs_r.sorted_filtered_merged.bed", "All_PASs_f.sorted_filtered_merged.bed")
bed_merge_command <- str_c("bedtools merge -d ",merge_dist ," -s -c 4,5,6 -o sum,sum,distinct -i ", filtere_sorted_bed_files[[1]], ' | awk \'{OFS="\t"}{print $1, $2, $3, $4, $5, $6}\' > ', filtere_sorted_merged_bed_files[[1]]) 
system(command = bed_merge_command)

bed_merge_command <- str_c("bedtools merge -d ",merge_dist ," -s -c 4,5,6 -o sum,sum,distinct -i ", filtere_sorted_bed_files[[2]], ' | awk \'{OFS="\t"}{print $1, $2, $3, $4, $5, $6}\' > ', filtere_sorted_merged_bed_files[[2]]) 
system(command = bed_merge_command)


# Summary ---------------------------------------------------------------------------------------

readNumber_tsv <- read_tsv('rawReadNumberInMillion.tsv')

readNumber_tsv %>% 
  group_by(type) %>% 
  summarise(sum_readN_inM = sum(readN_inM)*10E6)

# quality filtered pas number with read support number
nrow(pas_all_df) 
sum(as.integer(pas_all_df$name))

mergedPAS_r <- read_tsv('All_PASs_r.sorted_filtered_merged.bed', col_names = FALSE)
mergedPAS_f <- read_tsv('All_PASs_f.sorted_filtered_merged.bed', col_names = FALSE)

mergedPAS_tbl <- bind_rows(mergedPAS_r, mergedPAS_f)

nrow(pas_all_seq_filtered_df) 
sum(as.integer(pas_all_seq_filtered_df$name)) 


nrow(mergedPAS_tbl) 
sum(mergedPAS_tbl$X4) 

sum(as.integer(pas_all_seq_filtered_df$name)) == sum(mergedPAS_tbl$X4)


