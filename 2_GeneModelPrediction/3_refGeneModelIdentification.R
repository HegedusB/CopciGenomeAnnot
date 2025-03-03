# Compare the newly created Tama transcripts with a reference genome and add information about the splice junction. -------
# python3 /installs/SQANTI3/sqanti3_qc.py --force_id_ignore --gtf ../refAnnot.gtf ../refGenome.fsa -o tama_full_genome -c /minimapSJ/Copci_all_SJ.out.tab

# Create a dataset which includes the TAMA transcripts, SQANTI3 result ------------------------------------------------------

# read support
read_support <- read_tsv("tama_merge_out/Copci_merged_annots_read_support_NoSourceLine.txt")
read_support <- mutate(read_support, merge_key = str_replace(merge_trans_id, pattern = "^G", replacement = "PB."))

# sqanti  
sqanti_junctions <- read_tsv("sqanti2_test_out/Copci_tama_merge_junctions.txt")
sqanti_classification <- read_tsv("sqanti2_test_out/Copci_tama_merge_classification.txt")

setdiff(sqanti_classification$isoform, read_support$merge_key)
setdiff(read_support$merge_key, sqanti_classification$isoform)

sqanti_classification_m <- left_join(sqanti_classification, read_support, by = c("isoform" = "merge_key"))
sqanti_classification_m <- sqanti_classification_m[str_order(sqanti_classification_m$merge_trans_id, numeric = TRUE),]
write_tsv(sqanti_classification_m, "sqanti2_test_out/Copci_tama_merge_classification_m.txt")

# load packages ----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(parallel)
  library(rtracklayer)
  library(GenomicRanges)
})

# identify the most relevant transcript per gene model
sqanti_annot <- import("Copci_tama_merge_corrected.gtf")

sqanti_class_tbl <- read_tsv("Copci_tama_merge_classification_m.txt")

sqanti_class_lt <- split(sqanti_class_tbl, sqanti_class_tbl$merge_gene_id)

trans_length_p_filter_fct <- function(x) {
  trans_length_summary_tbl <-  tibble(isoform = x$isoform,
                                      merge_gene_id = x$merge_gene_id,
                                      merge_trans_id = x$merge_trans_id,
                                      length = x$length,
                                      trans_mean_length = mean(x$length),
                                      trans_median_length = median(x$length),
                                      trans_geom_mean_length = exp(mean(log(x$length))), 
                                      trans_length_p = x$length / max(x$length),
                                      trans_read_count = x$trans_read_count,
                                      trans_read_sup_p = x$trans_read_count / max(x$trans_read_count))
  return(trans_length_summary_tbl)
}


trans_length_p_lt <- lapply(sqanti_class_lt, FUN = trans_length_p_filter_fct)


trans_legth_filter_fct <- function(x, med_length_th = 0.80) {
  trans_lfiltered <- x[x$length >= (unique(x$trans_geom_mean_length) * med_length_th),]
  
  longest <- trans_lfiltered$length == max(trans_lfiltered$length)
  longest <- trans_lfiltered[longest,]
  if(nrow(longest) > 1) {
    longest <- longest[longest$trans_read_count == max(longest$trans_read_count),]
  }
  
  largest <- trans_lfiltered$trans_read_count == max(trans_lfiltered$trans_read_count) 
  largest <- trans_lfiltered[largest,]
  if(nrow(largest) > 1) {
    largest <- largest[largest$length == max(largest$length),]
  }
  
  if((largest$trans_read_count %/% longest$trans_read_count) >=4) {
    largest$tag <- "largest"
    return(largest)
  } else {
    longest$tag <- "longest"
    return(longest)
  }
}

trans_length_filter_lt <- lapply(trans_length_p_lt, FUN = trans_legth_filter_fct)


trans_length_filter_tbl <- bind_rows(trans_length_filter_lt)
trans_length_filter_tbl <- trans_length_filter_tbl[str_order(trans_length_filter_tbl$isoform, numeric = TRUE),]

sqanti_annot_filter <- sqanti_annot[sqanti_annot$transcript_id %in% trans_length_filter_tbl$isoform]
export(sqanti_annot_filter, "Copci_tama_merge_corrected_filter1.gtf")

# Separate over clustered TAMA gene models ---------------------------------------------


annot <- import("Copci_tama_merge_corrected.gtf", format = "gtf")
annot_df <- as.data.frame(annot)

annot_tscript <- filter(annot_df, type == "transcript")
annot_tscript$gene_id <- str_extract(annot_tscript$transcript_id, pattern = "PB\\.[0-9]+")

annot_tscript_lt <- split(annot_tscript, annot_tscript$gene_id)

tscript_clusterer_fct <- function(z) {
  annot_tscript_gr <- makeGRangesFromDataFrame(z, keep.extra.columns = TRUE)
  
  hits <- findOverlaps(annot_tscript_gr, drop.self=TRUE, drop.redundant=TRUE)
  
  x <- annot_tscript_gr[queryHits(hits)]
  y <- annot_tscript_gr[subjectHits(hits)]
  relative_overlap <- width(pintersect(x, y)) / pmin(width(x), width(y))
  hits <- hits[relative_overlap >= 0.5]
  
  extractClustersFromSelfHits <- function(hits)
  {
    stopifnot(is(hits, "Hits"))
    N <- queryLength(hits)
    stopifnot(N == subjectLength(hits))
    h <- union(hits, t(hits))
    qh <- queryHits(h)
    sh <- subjectHits(h)
    cid <- cid0 <- seq_len(N)  # cluster ids
    while (TRUE) {
      cid2 <- pmin(cid, selectHits(h, "first"))
      if (identical(cid2, cid))
        break
      cid <- cid2
      h <- Hits(qh, cid[sh], N, N)
    }
    unname(splitAsList(cid0, cid))
  }
  
  clusters <- extractClustersFromSelfHits(hits)
  
  clusters_df <- data.frame(clusters)
  clusters_df$gene_id <- z[clusters_df$value, "gene_id"]
  clusters_df$transcript_id <- z[clusters_df$value, "transcript_id"]
  clusters_df <- clusters_df[,-c(2:3)]
  
  return(clusters_df)
}


transcript_cls_lt <- mclapply(annot_tscript_lt, FUN = tscript_clusterer_fct, mc.cores = 20)
transcript_cls_tbl <- bind_rows(transcript_cls_lt)

write_tsv(transcript_cls_tbl, "Copci_tama_merge_corrected_transcript_cls.tsv")

