library(readr)
library(stringr)
library(rtracklayer)
library(dplyr)
library(Biostrings)
library(tidyr)
library(parallel)
library(BSgenome)



wd <- "/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/remap_annotatun_on_full_genome/jgi_PacBio_merged_genome"
setwd(wd)


# map the different annotation transcript on the reference genome an transform the resulting bam file into genePred format. 

# bamToBed -i CopciAB_jgiPacBio_nanopore_annot_s.bam -bed12 > CopciAB_jgiPacBio_nanopore_annot_s.bed
# bedToGenePred CopciAB_jgiPacBio_nanopore_annot_s.bed CopciAB_jgiPacBio_nanopore_annot_s.genePred


tscript_groups_lt <- list() 

# genome 
ref_genome <- readDNAStringSet("../../ref_genome/fullgenome_jgi/merged_genome/CopciAB_new_jgi_20220113.fasta")
ref_annot <- read.table("/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/isoforom_annotation_fullg_dir/new_ref_20210415/CopciAB_ref_20210415.genePred",
                        stringsAsFactors = FALSE)

# get coverage data -------------------------------------------------------------------------------

bam_path <- "Copci_jgiPacBio_all_ONT_PacBio_a6k22l14s2_s.bam"
cov_path <- str_replace(bam_path, pattern = ".bam", replacement = "_cover.bed")
read_cover_cmd <- str_c("bedtools genomecov -d -split -ibam ", bam_path, " > ", cov_path)
system(command = read_cover_cmd)

read_cover <- read_tsv(cov_path, col_names = F)
read_cover_lt <- split(read_cover, read_cover$X1)

genePred <- read.table("CopciAB_jgiPacBio_nanopore_annot_s.genePred",
                       stringsAsFactors = FALSE)

genePred_lt <- split(genePred, genePred$V1)
genePred_lt <- lapply(genePred_lt, function(x) return(mutate(x, V1 = str_c(V1, "_r", seq_len(length(V1))))))
genePred <- bind_rows(genePred_lt)
genePred_lt <- split(genePred, genePred$V1)


longergenom_index_0 <- !sapply(genePred_lt, function(x) {return(width(ref_genome[x$V2]) > x$V5)})

length(unique(str_remove(names(genePred_lt), pattern = "_r[0-9]+$")))
not_mapped <- setdiff(ref_annot$V1, unique(str_remove(names(genePred_lt), pattern = "_r[0-9]+$")))

tscript_groups_lt$not_mapped <- not_mapped


tscript_groups_lt$no_splice_site <- names(genePred_lt)[sapply(genePred_lt, "[[", "V8") == 1]
length(tscript_groups_lt$no_splice_site)

genePred_withsplice_lt <- genePred_lt[sapply(genePred_lt, "[[", "V8") > 1] 

genePred_cover_lt <- mclapply(genePred_withsplice_lt, function(x) {
  if(x$V3 == "-") {
    start_pos <- unlist(str_split(x$V9, pattern=","))
    start_pos <- as.integer(start_pos[str_detect(start_pos, pattern = "")]) + 1
    start_tbl <- tibble(start_pos = start_pos,
                        start_cov = read_cover_lt[[x$V2]][start_pos,"X3", drop=T])
    start_tbl <- start_tbl[-1,]
    
    end_pos <- unlist(str_split(x$V10, pattern=","))
    end_pos <- as.integer(end_pos[str_detect(end_pos, pattern = "")])
    end_tbl <- tibble(end_pos = end_pos,
                      end_cov = read_cover_lt[[x$V2]][end_pos,"X3", drop=T])
    end_tbl <- end_tbl[-nrow(end_tbl),]
    result <- cbind(start_tbl, end_tbl)
    result$strand <- x$V3
    return(result)
  } else {
    start_pos <- unlist(str_split(x$V9, pattern=","))
    start_pos <- as.integer(start_pos[str_detect(start_pos, pattern = "")]) + 1
    start_tbl <- tibble(start_pos = start_pos,
                        start_cov = read_cover_lt[[x$V2]][start_pos,"X3", drop=T])
    start_tbl <- start_tbl[-1,]
    
    end_pos <- unlist(str_split(x$V10, pattern=","))
    end_pos <- as.integer(end_pos[str_detect(end_pos, pattern = "")])
    end_tbl <- tibble(end_pos = end_pos,
                      end_cov = read_cover_lt[[x$V2]][end_pos,"X3", drop=T])
    end_tbl <- end_tbl[-nrow(end_tbl),]
    result <- cbind(start_tbl, end_tbl)
    result$strand <- x$V3
    return(result)
  }
}, mc.cores = 40)
genePred_cover_tbl <- bind_rows(genePred_cover_lt, .id = "transcript_id")

g_overhang_index <- names(genePred_cover_lt)[sapply(genePred_cover_lt, function(x) any(is.na(x)))] 
genePred_cover_m1_lt <- genePred_cover_lt[!sapply(genePred_cover_lt, function(x) any(is.na(x)))] 
genePred_cover_m1_tbl <- bind_rows(genePred_cover_m1_lt, .id = "transcript_id")

# genepred splice site ------------------------
genePred_withsplice_nogoverh_lt <- genePred_withsplice_lt[!names(genePred_withsplice_lt) %in% g_overhang_index]
splice_site_lt <- mclapply(genePred_withsplice_nogoverh_lt, function(x) {
  
  start_pos <- as.integer(unlist(str_extract_all(x$V9, pattern="[0-9]+")))
  start_tbl <- tibble(start = start_pos-1,
                      end = start_pos,
                      seqnames = x$V2,
                      strand = x$V3)
  start_tbl <- start_tbl[-1,]
  start_gr <- makeGRangesFromDataFrame(start_tbl, keep.extra.columns = TRUE)
  start_tbl$splice <- as.character(BSgenome::getSeq(ref_genome, start_gr))
  
  end_pos <- as.integer(unlist(str_extract_all(x$V10, pattern="[0-9]+"))) + 1
  end_tbl <- tibble(start = end_pos,
                    end = end_pos + 1,
                    seqnames = x$V2,
                    strand = x$V3)
  end_tbl <- end_tbl[-nrow(end_tbl),]
  end_gr <- makeGRangesFromDataFrame(end_tbl, keep.extra.columns = TRUE)
  end_tbl$splice <- as.character(BSgenome::getSeq(ref_genome, end_gr))
  if(x$V3 == "-") {
    start_tbl$index <- "start"
    end_tbl$index <- "end"
  } else {
    start_tbl$index <- "end"
    end_tbl$index <- "start"
  }
  result <- rbind(start_tbl, end_tbl)
  return(result)
}, mc.cores = 40)
sum(sapply(splice_site_lt, function(x) is.null(nrow(x))))
splice_site_tbl <- bind_rows(splice_site_lt, .id = "transcript_id")

splice_site_onlysites_lt <- lapply(splice_site_lt, function(x) {
  index_lt <- split(x, x$index)
  site_tbl <- bind_cols(index_lt$start["splice"], index_lt$end["splice"]) 
  site_tbl$splice_c <- apply(site_tbl, 1, function(x) return(str_c(x, collapse = "-")))
  site_tbl$strand <- unique(x$strand)
  return(site_tbl)
})
splice_site_onlysites_tbl <- bind_rows(splice_site_onlysites_lt, .id = "transcript_id")


# splice cover / splice site unification
genePred_cover_m1_lt$`CC2G_007305-T1_t0`
splice_site_onlysites_lt$`CC2G_007305-T1_t0`
identical(names(genePred_cover_m1_lt), names(splice_site_onlysites_lt))

genePred_cover_splice_site_tbl <- cbind(genePred_cover_m1_tbl, splice_site_onlysites_tbl)
genePred_cover_splice_site_tbl[c(6,7,8,9)] <- NULL
genePred_cover_splice_site_tbl <- genePred_cover_splice_site_tbl[c(1,2,4,3,5,6,7)]

exon_n <- sapply(genePred_withsplice_nogoverh_lt, "[[", "V8")
genePred_cover_splice_site_tbl$exon_n <- exon_n[match(genePred_cover_splice_site_tbl$transcript_id, names(exon_n))]

filter(genePred_cover_splice_site_tbl, transcript_id=="G343.19_t0_r1")


genePred_cover_splice_site_tbl %>% 
  group_by(transcript_id) %>% 
  summarise(u_splice_s = str_c(names(table(splice_c)), collapse = "_"),
            u_splice_n = str_c(table(splice_c), collapse = "_"),
            splice_comby_n = length(unique(splice_c))) -> genePred_cover_splice_site_summary_tbl


table(genePred_cover_splice_site_summary_tbl$splice_comby_n == 1)
barplot(table(genePred_cover_splice_site_summary_tbl$splice_comby_n == 1))

library(ggplot2)

filter(genePred_cover_splice_site_summary_tbl, splice_comby_n == 1) %>% 
  count(u_splice_s) %>% 
  arrange(desc(n)) %>% 
  mutate(u_splice_s = factor(u_splice_s, levels = u_splice_s, ordered = T)) %>% 
  ggplot(aes(x=u_splice_s, y=n)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90))

# Error detection ------------------------------

mapped_annot <- read.table("nanopore_annot_match.tsv", comment.char="@", stringsAsFactors = FALSE)
mapped_annot_lt <- split(mapped_annot, mapped_annot$V1)
mapped_annot_lt <- lapply(mapped_annot_lt, function(x) {
  x <- mutate(x, V1 = str_c(V1, "_r", seq_len(length(V1))))
  return(x)
})
mapped_annot <- bind_rows(mapped_annot_lt)


sum(str_detect(mapped_annot$V2, pattern="I"))

insertion_index <- str_detect(mapped_annot$V2, pattern="I")
insertion_tscript_tbl <- mapped_annot[insertion_index,]
str_extract_all(insertion_tscript_tbl$V2, pattern = "[0-9]+I")

length(unique(str_remove(insertion_tscript_tbl$V1, pattern = "_r[0-9]+$"))) 
table(sapply(str_split(unique(str_remove(insertion_tscript_tbl$V1, pattern = "_r[0-9]+$")), pattern="_|\\.|[0-9]+"), "[", 1))
barplot(table(sapply(str_split(unique(str_remove(insertion_tscript_tbl$V1, pattern = "_r[0-9]+$")), pattern="_|\\.|[0-9]+"), "[", 1)))

genePred_cover_splice_site_tbl$insert_index <- genePred_cover_splice_site_tbl$transcript_id %in% insertion_tscript_tbl$V1
length(unique(filter(genePred_cover_splice_site_tbl, insert_index)$transcript_id))

insertion_tscript_splice_lt <- str_extract_all(insertion_tscript_tbl$V2, pattern="[0-9]+N[0-9]+I[0-9]+M|[0-9]+M[0-9]+I[0-9]+N")
insertion_tscript_splice_index <- sapply(insertion_tscript_splice_lt, function(x) length(x) > 0)
sum(insertion_tscript_splice_index)
insertion_tscript_splice_tbl <- insertion_tscript_tbl[insertion_tscript_splice_index,]
length(unique(insertion_tscript_splice_tbl$V1))

genePred_cover_splice_site_tbl$insert_next_splice_index <- genePred_cover_splice_site_tbl$transcript_id %in% insertion_tscript_splice_tbl$V1
length(unique(filter(genePred_cover_splice_site_tbl, insert_next_splice_index)$transcript_id))


insertion_tscript_no_splice_lt <- str_extract_all(insertion_tscript_tbl$V2, pattern="[0-9]+M[0-9]+I[0-9]+M")
insertion_tscript_no_splice_index <- sapply(insertion_tscript_no_splice_lt, function(x) length(x) > 0)
sum(insertion_tscript_no_splice_index)
insertion_tscript_no_splice_tbl <- insertion_tscript_tbl[insertion_tscript_no_splice_index,]
length(unique(insertion_tscript_no_splice_tbl$V1))


sum(str_detect(mapped_annot$V2, pattern="D"))

deletion_index <- str_detect(mapped_annot$V2, pattern="D")
deletion_tscript_tbl <- mapped_annot[deletion_index,]
length(unique((deletion_tscript_tbl$V1)))

genePred_cover_splice_site_tbl$del_index <- genePred_cover_splice_site_tbl$transcript_id %in% deletion_tscript_tbl$V1


deletion_tscript_no_splice_lt <- str_extract_all(deletion_tscript_tbl$V2, pattern="[0-9]+M[0-9]+D[0-9]+M")
deletion_tscript_no_splice_index <- sapply(deletion_tscript_no_splice_lt, function(x) length(x) > 0)
sum(deletion_tscript_no_splice_index)
deletion_tscript_no_splice_tbl <- deletion_tscript_tbl[deletion_tscript_no_splice_index,]

setdiff(str_remove(genePred_cover_splice_site_tbl$transcript_id, pattern = "_r[0-9]+$"),ref_annot$V1)
genePred_cover_splice_site_tbl$ref_genome_exon_n <- ref_annot$V8[match(str_remove(genePred_cover_splice_site_tbl$transcript_id, pattern = "_r[0-9]+$"), ref_annot$V1)]

genePred_cover_splice_site_tbl <- left_join(genePred_cover_splice_site_tbl, genePred[c("V1", "V2")], by = c("transcript_id" = "V1"))

select(genePred_cover_splice_site_tbl, transcript_id, V2, everything()) %>% dplyr::rename("seqnames"="V2") -> genePred_cover_splice_site_tbl

write_tsv(genePred_cover_splice_site_tbl, "CopciAB_jgiPacBio_nanopore_annot_s_corrected_20220117_slice_info.tsv")

# summary ----------------------------------------------------------------------

annot_info_lt <- split(genePred_cover_splice_site_tbl, genePred_cover_splice_site_tbl$transcript_id)


annot_info_lt <- lapply(annot_info_lt, function(x) {
x <- annot_info_lt[["G5141.1_t0_r1"]]
  cover_all <- c(x$start_cov, x$end_cov)
  q1 <- quantile(cover_all, 0.25)
  q3 <- quantile(cover_all, 0.75)
  iqr_cover <- q3-q1
  lower <- q1-1.5*iqr_cover 
  upper <- q3+1.5*iqr_cover 
  if(any(cover_all<lower)) {
    result <- str_c(cover_all[cover_all<lower], collapse = "_")
  } else {
    result <- NA
  }
  x$cover_outlier <- result
  return(x)
})


annot_info_lt <- lapply(annot_info_lt, function(x) {
  cover_all <- c(x$start_cov, x$end_cov)
  if(all(cover_all==0)) {
    result <- NA
  } else {
    result <- (max(cover_all)-min(cover_all))/max(cover_all)
  }
  x$cover_outlier_2 <- result
  return(x)
})

hist(sapply(annot_info_lt, function(x) unique(x$cover_outlier_2)), xlab="(max(x) - min(x))/max(x)", main="")
abline(v=0.8, col="red")



length(annot_info_lt)

exon_n_problem_index <- sapply(annot_info_lt, function(x) return(all(x$exon_n != x$ref_genome_exon_n)))
sum(exon_n_problem_index) 
insert_problem_index <- sapply(annot_info_lt, function(x) return(all(x$insert_index)))
sum(insert_problem_index)

sum(exon_n_problem_index&insert_problem_index)
cor(exon_n_problem_index,insert_problem_index)
cor.test(as.integer(exon_n_problem_index),as.integer(insert_problem_index))

insert_at_splices_problem_index <- sapply(annot_info_lt, function(x) return(all(x$insert_next_splice_index)))
sum(insert_at_splices_problem_index)

sum(exon_n_problem_index&insert_at_splices_problem_index)
cor(exon_n_problem_index,insert_at_splices_problem_index)
cor.test(as.integer(exon_n_problem_index),as.integer(insert_at_splices_problem_index))

splice_site_problem_index <- !sapply(annot_info_lt, function(x) return(all(x$splice_c %in% c("GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AG"))))
sum(splice_site_problem_index)
sum(exon_n_problem_index&splice_site_problem_index)
cor(exon_n_problem_index,insert_at_splices_problem_index)
cor.test(as.integer(exon_n_problem_index),as.integer(insert_at_splices_problem_index))

cover_outlier_problem_index <- !sapply(annot_info_lt, function(x) return(all(is.na(x$cover_outlier))))
sum(cover_outlier_problem_index)
sum(exon_n_problem_index&cover_outlier_problem_index)

cover2_outlier_problem_index <- sapply(annot_info_lt, function(x) return(all(x$cover_outlier_2>0.8)))
sum(cover2_outlier_problem_index, na.rm = TRUE)
sum(exon_n_problem_index&cover2_outlier_problem_index, na.rm = TRUE)

# problem summary
library(corrplot)


posible_problems_tbl <- data.frame(exon_n_problem_index=exon_n_problem_index,
                                   insert_problem_index = insert_problem_index,
                                   insert_at_splices_problem_index=insert_at_splices_problem_index,
                                   conserved_splice_site_problem_index=splice_site_problem_index,
                                   cover_outlier_problem_index=cover_outlier_problem_index)

rownames(posible_problems_tbl) <- names(annot_info_lt)

cor(posible_problems_tbl[!apply(posible_problems_tbl,1,function(x) any(is.na(x))),])
corrplot(cor(posible_problems_tbl[!apply(posible_problems_tbl,1,function(x) any(is.na(x))),]))

library(ComplexHeatmap)
library(UpSetR)
m1 <- make_comb_mat(posible_problems_tbl)
UpSet(m1)


all_good_index <- apply(!posible_problems_tbl,1,all)
sum(all_good_index)

posible_problems_tbl_filter1 <- posible_problems_tbl[!all_good_index,]

m2 <- make_comb_mat(posible_problems_tbl_filter1)
UpSet(m2)

annot_info_lt[all_good_index]


# mark for correction  ---------

comb_name(m1, readable = TRUE)
comb_name(m1)
comb_size(m1)
comb_degree(m1)
length(extract_comb(m1, "01000"))
comb_index <- extract_comb(m1, "01000")
comb_index_tbl <- posible_problems_tbl[rownames(posible_problems_tbl)%in%comb_index,]


only_inser_problem_tbl <- comb_index_tbl


length(extract_comb(m1, "00001"))
only_cover_problem_index <- extract_comb(m1, "00001")
only_cover_problem_tbl <- posible_problems_tbl[rownames(posible_problems_tbl)%in%only_cover_problem_index,]


# filter good annotations ------
length(extract_comb(m1, "00000"))
good_annot_index <- extract_comb(m1, "00000")
good_annot_tbl <- posible_problems_tbl[rownames(posible_problems_tbl)%in%good_annot_index,]


all_good_annot_index <- c(comb_index, only_cover_problem_index, good_annot_index)
tscript_groups_lt$good_annot <- all_good_annot_index


genePred_all_bed_lt <- genePred_withsplice_nogoverh_lt[!names(genePred_withsplice_nogoverh_lt) %in% all_good_annot_index]


genePred_all_corrected_lt <- lapply(genePred_all_bed_lt, function(x) {
  transcript_id <- str_remove(x$V1, pattern = "_r[0-9]+$")
  genePred_ref <- filter(ref_annot, V1 == transcript_id)
  genePred_new <- x
  
  V4_diff <- genePred_new$V4 - genePred_ref$V4
  new_V9 <- as.integer(unlist(str_extract_all(genePred_ref$V9, pattern = "[0-9]+"))) + V4_diff
  new_V9_c <- str_c(str_c(new_V9, collapse = ","),",")
  new_V10 <- as.integer(unlist(str_extract_all(genePred_ref$V10, pattern = "[0-9]+"))) + V4_diff
  new_V10_c <- str_c(str_c(new_V10, collapse = ","),",")
  genePred_new$V9 <- new_V9_c
  genePred_new$V10 <- new_V10_c
  genePred_new$V8 <- length(new_V9)
  genePred_new["V5"] <- new_V10[length(new_V10)]
  genePred_new[c("V6", "V7")] <- 0
  return(genePred_new)
})

genePred_all_corrected_tbl <- bind_rows(genePred_all_corrected_lt)

out_name <- "Corrected_spliced_tscripts"
write_tsv(genePred_all_corrected_tbl, str_c(out_name, ".genePred"), col_names = F,)


togtf_cmd <- str_c('/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/tama_out_refannot/sep_joint_mapping_refannot_tama/genePredToGtf "file"  ',
                   " ", str_c(out_name, ".genePred"),
                   " ", str_c(out_name, ".gtf"))
system(command = togtf_cmd)


# test
bed_index <- names(genePred_all_bed_lt)[str_detect(names(genePred_all_bed_lt), pattern="^G")][[3]]
genePred_all_bed_lt[bed_index]

# test corrected annotatins ---------------------------------------------
sum(!sapply(genePred_all_corrected_lt, function(x) {return(width(ref_genome[x$V2]) > x$V5)}))
longergenom_index <- !sapply(genePred_all_corrected_lt, function(x) {return(width(ref_genome[x$V2]) > x$V5)})
genePred_all_corrected_lt <- genePred_all_corrected_lt[!longergenom_index]

tscript_groups_lt$corrected_overhang <- names(longergenom_index)[longergenom_index]
tscript_groups_lt$corrected <- names(longergenom_index)[!longergenom_index]

corrected_splice_site_lt <- mclapply(genePred_all_corrected_lt, function(x) {
  
  start_pos <- as.integer(unlist(str_extract_all(x$V9, pattern="[0-9]+")))
  start_tbl <- tibble(start = start_pos-1,
                      end = start_pos,
                      seqnames = x$V2,
                      strand = x$V3)
  start_tbl <- start_tbl[-1,]
  start_gr <- makeGRangesFromDataFrame(start_tbl, keep.extra.columns = TRUE)
  start_tbl$splice <- as.character(BSgenome::getSeq(ref_genome, start_gr))
  
  end_pos <- as.integer(unlist(str_extract_all(x$V10, pattern="[0-9]+"))) + 1
  end_tbl <- tibble(start = end_pos,
                    end = end_pos + 1,
                    seqnames = x$V2,
                    strand = x$V3)
  end_tbl <- end_tbl[-nrow(end_tbl),]
  end_gr <- makeGRangesFromDataFrame(end_tbl, keep.extra.columns = TRUE)
  end_tbl$splice <- as.character(BSgenome::getSeq(ref_genome, end_gr))
  if(x$V3 == "-") {
    start_tbl$index <- "start"
    end_tbl$index <- "end"
  } else {
    start_tbl$index <- "end"
    end_tbl$index <- "start"
  }
  result <- rbind(start_tbl, end_tbl)
  return(result)
}, mc.cores = 40)
sum(sapply(corrected_splice_site_lt, function(x) is.null(nrow(x))))
corrected_splice_site_tbl <- bind_rows(corrected_splice_site_lt, .id = "transcript_id")

corrected_splice_site_onlysites_lt <- lapply(corrected_splice_site_lt, function(x) {
  index_lt <- split(x, x$index)
  site_tbl <- bind_cols(index_lt$start["splice"], index_lt$end["splice"]) 
  site_tbl$splice_c <- apply(site_tbl, 1, function(x) return(str_c(x, collapse = "-")))
  site_tbl$strand <- unique(x$strand)
  return(site_tbl)
})
corrected_splice_site_onlysites_tbl <- bind_rows(corrected_splice_site_onlysites_lt, .id = "transcript_id")
corrected_splice_site_onlysites_lt[str_detect(names(corrected_splice_site_onlysites_lt), pattern="^G")]


sum(sapply(corrected_splice_site_onlysites_lt, function(x) return(all(x$splice_c %in% c("GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AG")))))

# save data ----

save(tscript_groups_lt, 
     genePred_all_corrected_lt,
     ref_annot, 
     genePred, 
     file = "Corrected_transcript_data_20200120.RData")


# CDS prediction ---------------

new_genePred_lt <- list(filter(genePred, V1 %in% tscript_groups_lt$good_annot),
                        filter(genePred, V1 %in% tscript_groups_lt$no_splice_site),
                        filter(genePred, V1 %in% tscript_groups_lt$corrected_overhang),
                        bind_rows(genePred_all_corrected_lt))

new_genePred <- bind_rows(new_genePred_lt)
new_genePred[c("V6", "V7")] <- 0

write_tsv(new_genePred, "CopciAB_new_annot_20220121.genePred", col_names = FALSE)

togtf_cmd <- str_c('/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/tama_out_refannot/sep_joint_mapping_refannot_tama/genePredToGtf "file"  ',
                   " ", "CopciAB_new_annot_20220121.genePred",
                   " ", str_replace("CopciAB_new_annot_20220121.genePred", pattern = "genePred", "gtf"))
system(command = togtf_cmd)

library(GenomicFeatures)

new_annot_gr <- import("CopciAB_new_annot_20220121.gtf")
new_annot_TxDb <- makeTxDbFromGRanges(new_annot_gr)

new_annot_TxDb_lt <- GenomicFeatures::as.list(new_annot_TxDb)
tx_name <- new_annot_TxDb_lt$transcripts$tx_name

ebt <- exonsBy(new_annot_TxDb, by = "tx")
names(ebt) <- tx_name

transcript_seq <- extractTranscriptSeqs(ref_genome, ebt)



library(ORFik)

man_ORFs <- findORFs(transcript_seq, 
                     startCodon = "ATG",
                     stopCodon = stopDefinition(1),
                     longestORF = TRUE, minimumLength = 30)

names(man_ORFs) <- names(transcript_seq)[as.integer(names(man_ORFs))]

man_ORF_prot_coord_lt <- mclapply(man_ORFs, function(x) {
  max_length_ORF <- x[width(x) == max(width(x))][1,]
  return(as.data.frame(max_length_ORF))
}, mc.cores = 10)

man_ORF_prot_coord_tbl <- bind_rows(man_ORF_prot_coord_lt, .id = "transcript_id")


new_genePred_lt <- split(new_genePred, new_genePred$V1)
new_genePred_CDS_lt <- mclapply(new_genePred_lt, function(x) {
  
  if(any(str_detect(man_ORF_prot_coord_tbl$transcript_id, x$V1))) {
    CDS_range <- filter(man_ORF_prot_coord_tbl, transcript_id == x$V1)
    if(x$V3 == "+") {
      exon_range <- t(rbind(as.integer(str_extract_all(x$V9, pattern = "[0-9]+", simplify = T)),
      as.integer(str_extract_all(x$V10, pattern = "[0-9]+", simplify = T))))
      exon_range_lt <- apply(exon_range,1, function(x) IRanges(x[[1]]+1, x[[2]]))
      tscript_range <- unlist(lapply(exon_range_lt, function(x) seq(start(x), end(x))))
      start_end <- tscript_range[c(CDS_range$start, CDS_range$end)]
      x$V6 <- start_end[[1]]-1
      x$V7 <- start_end[[2]]
    }  else {
      exon_range <- t(rbind(as.integer(str_extract_all(x$V9, pattern = "[0-9]+", simplify = T)),
                            as.integer(str_extract_all(x$V10, pattern = "[0-9]+", simplify = T))))
      exon_range_lt <- apply(exon_range,1, function(x) IRanges(x[[1]]+1, x[[2]]))
      tscript_range <- unlist(lapply(exon_range_lt, function(x) seq(start(x), end(x))))
      tscript_range_rev <- rev(tscript_range)
      start_end <- tscript_range_rev[c(CDS_range$start, CDS_range$end)]
      x$V6 <- start_end[[2]]-1
      x$V7 <- start_end[[1]]
      }
  }
  return(x)
}, mc.cores = 10)

new_genePred_CDS_tbl <- bind_rows(new_genePred_CDS_lt)

write_tsv(new_genePred_CDS_tbl, "CopciAB_new_annot_CDS_20220121.genePred", col_names = FALSE)

togtf_cmd <- str_c('/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/tama_out_refannot/sep_joint_mapping_refannot_tama/genePredToGtf "file" -utr ',
                   " ", "CopciAB_new_annot_CDS_20220121.genePred",
                   " ", str_replace("CopciAB_new_annot_CDS_20220121.genePred", pattern = "genePred", "gtf"))
system(command = togtf_cmd)


new_annot_CDS_gr <- import("CopciAB_new_annot_CDS_20220121.gtf")
new_annot_CDS_TxDb <- makeTxDbFromGRanges(new_annot_CDS_gr)

new_annot_CDS_TxDb_lt <- GenomicFeatures::as.list(new_annot_CDS_TxDb)
tx_CDS_name <- new_annot_CDS_TxDb_lt$transcripts$tx_name

cbt <- cdsBy(new_annot_CDS_TxDb, by = "tx")
names(cbt) <- tx_CDS_name[as.integer(names(cbt))]

cds_seq <- extractTranscriptSeqs(ref_genome, cbt)
prot_seq <- translate(cds_seq)
prot_seq["G7757.13_t0_r1"]

writeXStringSet(prot_seq, "CopciAB_new_annot_CDS_20220121.prot.fasta")


# BUSCO comparison ----
old_BUSCO <- read_tsv("/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/isoforom_annotation_fullg_dir/new_ref_20210415/run_run_BUSCO_CopciAB_ref_20210415_odb10/full_table_run_BUSCO_CopciAB_ref_20210415_odb10.tsv",
                      comment = "#",
                      col_names = FALSE)

new_BUSCO <- read_tsv("/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/run_BUSCO_CopciAB_new_20220121_odb10/full_table_BUSCO_CopciAB_new_20220121_odb10.tsv", 
                      comment = "#",
                      col_names = FALSE)

new_BUSCO <- mutate(new_BUSCO, index = str_remove(X3, pattern = "_r[0-9]+$"))

new_BUSCO_missing <- filter(new_BUSCO, X2=="Missing")

old_BUSCO_missing_short <- filter(old_BUSCO, X1 %in% new_BUSCO_missing$X1)
filter(old_BUSCO_missing_short, !is.na(X3)) %>% knitr::kable()

############################ 2 round ############################################

new_genePred_CDS_lt <- split(new_genePred_CDS_tbl, new_genePred_CDS_tbl$V1)
new_genePred_CDS_withsplice_lt <- new_genePred_CDS_lt[sapply(new_genePred_CDS_lt, "[[", "V8") > 1] # 15150

# genepred splice site ----

new_splice_site_lt <- mclapply(new_genePred_CDS_withsplice_lt, function(x) {
  start_pos <- as.integer(unlist(str_extract_all(x$V9, pattern="[0-9]+")))
  start_tbl <- tibble(start = start_pos-1,
                      end = start_pos,
                      seqnames = x$V2,
                      strand = x$V3)
  start_tbl <- start_tbl[-1,]
  start_gr <- makeGRangesFromDataFrame(start_tbl, keep.extra.columns = TRUE)
  start_tbl$splice <- as.character(BSgenome::getSeq(ref_genome, start_gr))
  
  end_pos <- as.integer(unlist(str_extract_all(x$V10, pattern="[0-9]+"))) + 1
  end_tbl <- tibble(start = end_pos,
                    end = end_pos + 1,
                    seqnames = x$V2,
                    strand = x$V3)
  end_tbl <- end_tbl[-nrow(end_tbl),]
  end_gr <- makeGRangesFromDataFrame(end_tbl, keep.extra.columns = TRUE)
  end_tbl$splice <- as.character(BSgenome::getSeq(ref_genome, end_gr))
  if(x$V3 == "-") {
    start_tbl$index <- "start"
    end_tbl$index <- "end"
  } else {
    start_tbl$index <- "end"
    end_tbl$index <- "start"
  }
  result <- rbind(start_tbl, end_tbl)
  return(result)
}, mc.cores = 40)

sum(sapply(new_splice_site_lt, function(x) is.null(nrow(x))))
new_splice_site_tbl <- bind_rows(new_splice_site_lt, .id = "transcript_id")

new_splice_site_onlysites_lt <- lapply(new_splice_site_lt, function(x) {
  index_lt <- split(x, x$index)
  site_tbl <- bind_cols(index_lt$start["splice"], index_lt$end["splice"]) 
  site_tbl$splice_c <- apply(site_tbl, 1, function(x) return(str_c(x, collapse = "-")))
  site_tbl$strand <- unique(x$strand)
  return(site_tbl)
})
new_splice_site_onlysites_tbl <- bind_rows(new_splice_site_onlysites_lt, .id = "transcript_id")

new_splice_site_problem_index <- !sapply(new_splice_site_onlysites_lt, function(x) return(all(x$splice_c %in% c("GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AG"))))
sum(new_splice_site_problem_index)

sum(splice_site_problem_index)
sum(new_splice_site_problem_index)

intersect(names(splice_site_problem_index)[splice_site_problem_index], names(new_splice_site_problem_index)[new_splice_site_problem_index])
length(intersect(names(splice_site_problem_index)[splice_site_problem_index], names(new_splice_site_problem_index)[new_splice_site_problem_index])) # 178 átfedés
length(setdiff(names(splice_site_problem_index)[splice_site_problem_index], names(new_splice_site_problem_index)[new_splice_site_problem_index])) # csak a javítatlanban: 266
length(setdiff(names(new_splice_site_problem_index)[new_splice_site_problem_index],
       names(splice_site_problem_index)[splice_site_problem_index]))

filter(old_BUSCO_missing_short, !is.na(X3)) %>% knitr::kable()
filter(old_BUSCO_missing_short, !is.na(X3)) -> old_BUSCO_missing_short2
old_BUSCO_missing_short2 <- mutate(old_BUSCO_missing_short2, index = str_c(X3, "_r1"))
sum(old_BUSCO_missing_short2$index %in% setdiff(names(new_splice_site_problem_index)[new_splice_site_problem_index],
        names(splice_site_problem_index)[splice_site_problem_index]))

filter(new_BUSCO, X2 == "Fragmented")
sum(filter(new_BUSCO, X2 == "Fragmented")$X3 %in% setdiff(names(new_splice_site_problem_index)[new_splice_site_problem_index],
                                                      names(splice_site_problem_index)[splice_site_problem_index]))

reannot_index <- setdiff(names(new_splice_site_problem_index)[new_splice_site_problem_index],
                         names(splice_site_problem_index)[splice_site_problem_index])



genePred_all_bed2_lt <- genePred_withsplice_nogoverh_lt[names(genePred_withsplice_nogoverh_lt) %in% reannot_index]


genePred_all_corrected2_lt <- lapply(genePred_all_bed2_lt, function(x) {
  transcript_id <- str_remove(x$V1, pattern = "_r[0-9]+$")
  genePred_ref <- filter(ref_annot, V1 == transcript_id) 
  genePred_new <- x
  
  V4_diff <- genePred_new$V7 - genePred_ref$V7
  new_V9 <- as.integer(unlist(str_extract_all(genePred_ref$V9, pattern = "[0-9]+"))) + V4_diff 
  new_V9_c <- str_c(str_c(new_V9, collapse = ","),",")
  new_V10 <- as.integer(unlist(str_extract_all(genePred_ref$V10, pattern = "[0-9]+"))) + V4_diff 
  new_V10_c <- str_c(str_c(new_V10, collapse = ","),",")
  genePred_new$V9 <- new_V9_c
  genePred_new$V10 <- new_V10_c
  genePred_new$V8 <- length(new_V9)
  genePred_new["V4"] <- new_V9[1]
  genePred_new["V5"] <- new_V10[length(new_V10)]
  genePred_new[c("V6", "V7")] <- 0
  return(genePred_new)
})
corrected2_overhang_index <- genePred_all_corrected2_tbl[str_detect(genePred_all_corrected2_tbl$V9, pattern="-"),]$V1
genePred_all_corrected2_lt <- genePred_all_corrected2_lt[!names(genePred_all_corrected2_lt) %in% corrected2_overhang_index]

genePred_all_corrected2_tbl <- bind_rows(genePred_all_corrected2_lt)

filter(genePred_all_corrected2_tbl, V1 == "G844.14_t0_r1")

out_name <- "Corrected2_spliced_tscripts"
write_tsv(genePred_all_corrected2_tbl, str_c(out_name, ".genePred"), col_names = F)


togtf_cmd <- str_c('/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/tama_out_refannot/sep_joint_mapping_refannot_tama/genePredToGtf "file"  ',
                   " ", str_c(out_name, ".genePred"),
                   " ", str_c(out_name, ".gtf"))
system(command = togtf_cmd)

# update annotation -------------------------------------------------------------


new_genePred_short <- new_genePred[!new_genePred$V1 %in% genePred_all_corrected2_tbl$V1,]
new2_genePred <- bind_rows(new_genePred_short, genePred_all_corrected2_tbl)


write_tsv(new2_genePred, "CopciAB_new_annot_20220124.genePred", col_names = FALSE)

togtf_cmd <- str_c('/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/tama_out_refannot/sep_joint_mapping_refannot_tama/genePredToGtf "file"  ',
                   " ", "CopciAB_new_annot_20220124.genePred",
                   " ", str_replace("CopciAB_new_annot_20220124.genePred", pattern = "genePred", "gtf"))
system(command = togtf_cmd)


library(GenomicFeatures)

new_annot_gr <- import("CopciAB_new_annot_20220124.gtf")
new_annot_TxDb <- makeTxDbFromGRanges(new_annot_gr)

new_annot_TxDb_lt <- GenomicFeatures::as.list(new_annot_TxDb)
tx_name <- new_annot_TxDb_lt$transcripts$tx_name

ebt <- exonsBy(new_annot_TxDb, by = "tx")
names(ebt) <- tx_name

transcript_seq <- extractTranscriptSeqs(ref_genome, ebt)


library(ORFik)

man_ORFs <- findORFs(transcript_seq, 
                     startCodon = "ATG",
                     stopCodon = stopDefinition(1),
                     longestORF = TRUE, minimumLength = 30)

names(man_ORFs) <- names(transcript_seq)[as.integer(names(man_ORFs))]


man_ORF_prot_coord_lt <- mclapply(man_ORFs, function(x) {
  max_length_ORF <- x[width(x) == max(width(x))][1,]
  return(as.data.frame(max_length_ORF))
}, mc.cores = 10)

man_ORF_prot_coord_tbl <- bind_rows(man_ORF_prot_coord_lt, .id = "transcript_id")


new2_genePred_lt <- split(new2_genePred, new2_genePred$V1)
new2_genePred_CDS_lt <- mclapply(new2_genePred_lt, function(x) {
  
  if(any(str_detect(man_ORF_prot_coord_tbl$transcript_id, x$V1))) {
    CDS_range <- filter(man_ORF_prot_coord_tbl, transcript_id == x$V1)
    if(x$V3 == "+") {

      exon_range <- t(rbind(as.integer(str_extract_all(x$V9, pattern = "[0-9]+", simplify = T)),
                            as.integer(str_extract_all(x$V10, pattern = "[0-9]+", simplify = T))))
      exon_range_lt <- apply(exon_range,1, function(x) IRanges(x[[1]]+1, x[[2]]))
      tscript_range <- unlist(lapply(exon_range_lt, function(x) seq(start(x), end(x))))
      start_end <- tscript_range[c(CDS_range$start, CDS_range$end)]
      x$V6 <- start_end[[1]]-1
      x$V7 <- start_end[[2]]
      
     }  else {

      exon_range <- t(rbind(as.integer(str_extract_all(x$V9, pattern = "[0-9]+", simplify = T)),
                            as.integer(str_extract_all(x$V10, pattern = "[0-9]+", simplify = T))))
      exon_range_lt <- apply(exon_range,1, function(x) IRanges(x[[1]]+1, x[[2]]))
      tscript_range <- unlist(lapply(exon_range_lt, function(x) seq(start(x), end(x))))
      tscript_range_rev <- rev(tscript_range)
      start_end <- tscript_range_rev[c(CDS_range$start, CDS_range$end)]
      x$V6 <- start_end[[2]]-1
      x$V7 <- start_end[[1]]
    }
  }
  return(x)
}, mc.cores = 10)

new2_genePred_CDS_tbl <- bind_rows(new2_genePred_CDS_lt)


write_tsv(new2_genePred_CDS_tbl, "CopciAB_new_annot_CDS_20220124.genePred", col_names = FALSE)

togtf_cmd <- str_c('/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/tama_out_refannot/sep_joint_mapping_refannot_tama/genePredToGtf "file" -utr ',
                   " ", "CopciAB_new_annot_CDS_20220124.genePred",
                   " ", str_replace("CopciAB_new_annot_CDS_20220124.genePred", pattern = "genePred", "gtf"))
system(command = togtf_cmd)


new_annot_CDS_gr <- import("CopciAB_new_annot_CDS_20220124.gtf")
new_annot_CDS_TxDb <- makeTxDbFromGRanges(new_annot_CDS_gr)

new_annot_CDS_TxDb_lt <- GenomicFeatures::as.list(new_annot_CDS_TxDb)
tx_CDS_name <- new_annot_CDS_TxDb_lt$transcripts$tx_name

cbt <- cdsBy(new_annot_CDS_TxDb, by = "tx")
names(cbt) <- tx_CDS_name[as.integer(names(cbt))]

cds_seq <- extractTranscriptSeqs(ref_genome, cbt)
prot_seq <- translate(cds_seq)
prot_seq["G6080.8_t0_r1"]

writeXStringSet(prot_seq, "CopciAB_new_annot_CDS_20220124.prot.fasta")


prot_path <- "CopciAB_new_annot_CDS_20220124.prot.fasta"
BUSCO_cmd <- str_c("run_BUSCO.py -f -i ", prot_path, 
                   " -o ", str_replace(prot_path, pattern = ".prot.fasta", replacement = "_odb10"),
                   " -l /work/balintb/Busco/Databases/v4/basidiomycota_odb10/ -c 10 -m prot")
system(command = BUSCO_cmd)

# BUSCO comparison ----
old_BUSCO <- read_tsv("/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/isoforom_annotation_fullg_dir/new_ref_20210415/run_run_BUSCO_CopciAB_ref_20210415_odb10/full_table_run_BUSCO_CopciAB_ref_20210415_odb10.tsv",
                      comment = "#",
                      col_names = FALSE)

new_BUSCO <- read_tsv("/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/remap_annotatun_on_full_genome/jgi_PacBio_merged_genome/run_CopciAB_new_annot_CDS_20220124_odb10/full_table_CopciAB_new_annot_CDS_20220124_odb10.tsv", 
                      comment = "#",
                      col_names = FALSE)

new_BUSCO <- mutate(new_BUSCO, index = str_remove(X3, pattern = "_r[0-9]+$"))

new_BUSCO_missing <- filter(new_BUSCO, X2=="Missing")

old_BUSCO_missing_short <- filter(old_BUSCO, X1 %in% new_BUSCO_missing$X1)
filter(old_BUSCO_missing_short, !is.na(X3)) %>% knitr::kable()

new_BUSCO_fragmented <- filter(new_BUSCO, X2=="Fragmented")


new_annot_CDS_gr[new_annot_CDS_gr$gene_id =="CC2G_009170-T1_t0_r1"]


############################ 3 round ############################################

new_genePred_CDS_lt <- split(new2_genePred_CDS_tbl, new2_genePred_CDS_tbl$V1)
new_genePred_CDS_withsplice_lt <- new_genePred_CDS_lt[sapply(new_genePred_CDS_lt, "[[", "V8") > 1] # 15150

# genepred splice site ----

new_splice_site_lt <- mclapply(new_genePred_CDS_withsplice_lt, function(x) {

  start_pos <- as.integer(unlist(str_extract_all(x$V9, pattern="[0-9]+")))
  start_tbl <- tibble(start = start_pos-1,
                      end = start_pos,
                      seqnames = x$V2,
                      strand = x$V3)
  start_tbl <- start_tbl[-1,]
  start_gr <- makeGRangesFromDataFrame(start_tbl, keep.extra.columns = TRUE)
  start_tbl$splice <- as.character(BSgenome::getSeq(ref_genome, start_gr))
  
  end_pos <- as.integer(unlist(str_extract_all(x$V10, pattern="[0-9]+"))) + 1
  end_tbl <- tibble(start = end_pos,
                    end = end_pos + 1,
                    seqnames = x$V2,
                    strand = x$V3)
  end_tbl <- end_tbl[-nrow(end_tbl),]
  end_gr <- makeGRangesFromDataFrame(end_tbl, keep.extra.columns = TRUE)
  end_tbl$splice <- as.character(BSgenome::getSeq(ref_genome, end_gr))
  if(x$V3 == "-") {
    start_tbl$index <- "start"
    end_tbl$index <- "end"
  } else {
    start_tbl$index <- "end"
    end_tbl$index <- "start"
  }
  result <- rbind(start_tbl, end_tbl)
  return(result)
}, mc.cores = 40)

sum(sapply(new_splice_site_lt, function(x) is.null(nrow(x))))
new_splice_site_tbl <- bind_rows(new_splice_site_lt, .id = "transcript_id")

new_splice_site_onlysites_lt <- lapply(new_splice_site_lt, function(x) {
  index_lt <- split(x, x$index)
  site_tbl <- bind_cols(index_lt$start["splice"], index_lt$end["splice"]) 
  site_tbl$splice_c <- apply(site_tbl, 1, function(x) return(str_c(x, collapse = "-")))
  site_tbl$strand <- unique(x$strand)
  return(site_tbl)
})
new_splice_site_onlysites_tbl <- bind_rows(new_splice_site_onlysites_lt, .id = "transcript_id")


new_splice_site_problem_index <- !sapply(new_splice_site_onlysites_lt, function(x) return(all(x$splice_c %in% c("GT-AG", "CT-AC", "GC-AG", "CT-GC", "AT-AC", "GT-AG"))))
sum(new_splice_site_problem_index) 

sum(splice_site_problem_index) 
sum(new_splice_site_problem_index) 

intersect(names(splice_site_problem_index)[splice_site_problem_index], names(new_splice_site_problem_index)[new_splice_site_problem_index])
length(intersect(names(splice_site_problem_index)[splice_site_problem_index], names(new_splice_site_problem_index)[new_splice_site_problem_index]))
length(setdiff(names(splice_site_problem_index)[splice_site_problem_index], names(new_splice_site_problem_index)[new_splice_site_problem_index]))
length(setdiff(names(new_splice_site_problem_index)[new_splice_site_problem_index],
               names(splice_site_problem_index)[splice_site_problem_index])) 

setdiff(names(new_splice_site_problem_index)[new_splice_site_problem_index],
        names(splice_site_problem_index)[splice_site_problem_index])[sapply(str_split(setdiff(names(new_splice_site_problem_index)[new_splice_site_problem_index],
               names(splice_site_problem_index)[splice_site_problem_index]),
          pattern = "[0-9]+|_"), "[[", 1) == "G"] %>% knitr::kable()

filter(old_BUSCO_missing_short, !is.na(X3)) %>% knitr::kable()
filter(old_BUSCO_missing_short, !is.na(X3)) -> old_BUSCO_missing_short2
old_BUSCO_missing_short2 <- mutate(old_BUSCO_missing_short2, index = str_c(X3, "_r1"))
sum(old_BUSCO_missing_short2$index %in% setdiff(names(new_splice_site_problem_index)[new_splice_site_problem_index],
                                                names(splice_site_problem_index)[splice_site_problem_index])) 

filter(new_BUSCO, X2 == "Fragmented")
sum(filter(new_BUSCO, X2 == "Fragmented")$X3 %in% setdiff(names(new_splice_site_problem_index)[new_splice_site_problem_index],
                                                          names(splice_site_problem_index)[splice_site_problem_index])) 

reannot_index <- setdiff(names(new_splice_site_problem_index)[new_splice_site_problem_index],
                         names(splice_site_problem_index)[splice_site_problem_index])




genees_splice_problem_lt <- new_genePred_CDS_withsplice_lt[new_splice_site_problem_index]
table(sapply(str_split(names(genees_splice_problem_lt), pattern="_|\\.|[0-9]+"), "[[", 1))


# quick check --------------

library(GenomicFeatures)

new_annot_gr <- import("CopciAB_new_annot_CDS_20220124_mHB1.gtf")
new_annot_gr$phase <- 0
new_annot_gr$phase <- as.integer(new_annot_gr$phase)
new_annot_TxDb <- makeTxDbFromGRanges(new_annot_gr)

new_annot_TxDb_lt <- GenomicFeatures::as.list(new_annot_TxDb)
tx_name <- new_annot_TxDb_lt$transcripts$tx_name

ebt <- exonsBy(new_annot_TxDb, by = "tx")
names(ebt) <- tx_name

transcript_seq <- extractTranscriptSeqs(ref_genome, ebt)



library(ORFik)

man_ORFs <- findORFs(transcript_seq, 
                     startCodon = "ATG",
                     stopCodon = stopDefinition(1),
                     longestORF = TRUE, minimumLength = 30)

names(man_ORFs) <- names(transcript_seq)[as.integer(names(man_ORFs))]


man_ORF_prot_coord_lt <- mclapply(man_ORFs, function(x) {
  max_length_ORF <- x[width(x) == max(width(x))][1,]
  return(as.data.frame(max_length_ORF))
}, mc.cores = 10)

man_ORF_prot_coord_tbl <- bind_rows(man_ORF_prot_coord_lt, .id = "transcript_id")

new_annot <- read.table("CopciAB_new_annot_CDS_20220124_mHB1.genePred", stringsAsFactors = FALSE)

new2_genePred_lt <- split(new_annot, new_annot$V1)
new2_genePred_CDS_lt <- mclapply(new2_genePred_lt, function(x) {
  
  if(any(str_detect(man_ORF_prot_coord_tbl$transcript_id, x$V1))) {
    CDS_range <- filter(man_ORF_prot_coord_tbl, transcript_id == x$V1)
    if(x$V3 == "+") {
      exon_range <- t(rbind(as.integer(str_extract_all(x$V9, pattern = "[0-9]+", simplify = T)),
                            as.integer(str_extract_all(x$V10, pattern = "[0-9]+", simplify = T))))
      exon_range_lt <- apply(exon_range,1, function(x) IRanges(x[[1]]+1, x[[2]]))
      tscript_range <- unlist(lapply(exon_range_lt, function(x) seq(start(x), end(x))))
      start_end <- tscript_range[c(CDS_range$start, CDS_range$end)]
      x$V6 <- start_end[[1]]-1
      x$V7 <- start_end[[2]]
    }  else {
      exon_range <- t(rbind(as.integer(str_extract_all(x$V9, pattern = "[0-9]+", simplify = T)),
                            as.integer(str_extract_all(x$V10, pattern = "[0-9]+", simplify = T))))
      exon_range_lt <- apply(exon_range,1, function(x) IRanges(x[[1]]+1, x[[2]]))
      tscript_range <- unlist(lapply(exon_range_lt, function(x) seq(start(x), end(x))))
      tscript_range_rev <- rev(tscript_range)
      start_end <- tscript_range_rev[c(CDS_range$start, CDS_range$end)]
      x$V6 <- start_end[[2]]-1
      x$V7 <- start_end[[1]]
    }
  }
  return(x)
}, mc.cores = 10)

new2_genePred_CDS_tbl <- bind_rows(new2_genePred_CDS_lt)

write_tsv(new2_genePred_CDS_tbl, "temp.genePred", col_names = FALSE)

togtf_cmd <- str_c('/SSD2/bhegedus/pj8_3genseq/Transcript_sequencing/Copci_ONT/all_good_pass_reads/tama_out_refannot/sep_joint_mapping_refannot_tama/genePredToGtf "file" -utr ',
                   " ", "temp.genePred",
                   " ", str_replace("temp.genePred", pattern = "genePred", "gtf"))
system(command = togtf_cmd)


new_annot_CDS_gr <- import("temp.gtf")
new_annot_CDS_TxDb <- makeTxDbFromGRanges(new_annot_CDS_gr)

new_annot_CDS_TxDb_lt <- GenomicFeatures::as.list(new_annot_CDS_TxDb)
tx_CDS_name <- new_annot_CDS_TxDb_lt$transcripts$tx_name

cbt <- cdsBy(new_annot_CDS_TxDb, by = "tx")
names(cbt) <- tx_CDS_name[as.integer(names(cbt))]

cds_seq <- extractTranscriptSeqs(ref_genome, cbt)
prot_seq <- translate(cds_seq)

writeXStringSet(prot_seq, "temp.prot.fasta")

file.rename("temp.prot.fasta", "CopciAB_new_annot_CDS_20220127.prot.fasta")
file.rename("temp.gtf", "CopciAB_new_annot_CDS_20220127.gtf")
file.rename("temp.genePred", "CopciAB_new_annot_CDS_20220127.genePred")


# Add splice junction coverage data ----


new2_genePred_CDS_tbl <- read.table("CopciAB_new_annot_CDS_20220127.genePred", stringsAsFactors = FALSE)
ref_genome <- readDNAStringSet("../../ref_genome/fullgenome_jgi/merged_genome/CopciAB_new_jgi_20220113.fasta")

new_genePred_CDS_lt <- split(new2_genePred_CDS_tbl, new2_genePred_CDS_tbl$V1)
new_genePred_CDS_withsplice_lt <- new_genePred_CDS_lt[sapply(new_genePred_CDS_lt, "[[", "V8") > 1]


new_splice_site_lt <- mclapply(new_genePred_CDS_withsplice_lt, function(x) {
  start_pos <- as.integer(unlist(str_extract_all(x$V9, pattern="[0-9]+")))
  start_tbl <- tibble(start = start_pos-1,
                      end = start_pos,
                      seqnames = x$V2,
                      strand = x$V3)
  start_tbl <- start_tbl[-1,]
  start_gr <- makeGRangesFromDataFrame(start_tbl, keep.extra.columns = TRUE)
  start_tbl$splice <- as.character(BSgenome::getSeq(ref_genome, start_gr))
  
  end_pos <- as.integer(unlist(str_extract_all(x$V10, pattern="[0-9]+"))) + 1
  end_tbl <- tibble(start = end_pos,
                    end = end_pos + 1,
                    seqnames = x$V2,
                    strand = x$V3)
  end_tbl <- end_tbl[-nrow(end_tbl),]
  end_gr <- makeGRangesFromDataFrame(end_tbl, keep.extra.columns = TRUE)
  end_tbl$splice <- as.character(BSgenome::getSeq(ref_genome, end_gr))
  if(x$V3 == "-") {
    start_tbl$index <- "start"
    end_tbl$index <- "end"
  } else {
    start_tbl$index <- "end"
    end_tbl$index <- "start"
  }
  result <- rbind(start_tbl, end_tbl)
  return(result)
}, mc.cores = 40)

sum(sapply(new_splice_site_lt, function(x) is.null(nrow(x))))
new_splice_site_tbl <- bind_rows(new_splice_site_lt, .id = "transcript_id")

new_splice_site_onlysites_lt <- mclapply(new_splice_site_lt, function(x) {

  index_lt <- split(x, x$index)
  if(all(x$strand == "+")) {
    site_tbl <- bind_cols(index_lt$start["splice"], index_lt$end["splice"]) 
    site_tbl$start <- index_lt$start$start
    site_tbl$end <- index_lt$end$end
    site_tbl <- tidyr::unite(site_tbl, splice_c, "splice...1", "splice...2", sep="-")
    site_tbl$strand <- unique(x$strand)
    site_tbl$seqnames <- unique(x$seqnames)
    return(site_tbl)
  } else {
    site_tbl <- bind_cols(index_lt$start["splice"], index_lt$end["splice"]) 
    site_tbl$start <- index_lt$end$start
    site_tbl$end <- index_lt$start$end
    site_tbl <- tidyr::unite(site_tbl, splice_c, "splice...1", "splice...2", sep="-")
    site_tbl$strand <- unique(x$strand)
    site_tbl$seqnames <- unique(x$seqnames)
    return(site_tbl)
  }
  
}, mc.cores = 20)
new_splice_site_onlysites_tbl <- bind_rows(new_splice_site_onlysites_lt, .id = "transcript_id")

write_tsv(new_splice_site_onlysites_tbl, "CopciAB_new_annot_CDS_20220127_splice_info.tsv")






