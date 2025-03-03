library(Biostrings)
library(stringr)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

options(width=300)

# genome short summary ----------------------------------------------------

refgenome <- readDNAStringSet("CopciAB_new_jgi_20220113_C13.fasta")
length(refgenome)

par(mar = c(7, 4, 2, 2) + 0.2)
barplot(width(refgenome), names.arg = names(refgenome), las=2)

telomere_5_lt <- lapply(refgenome, function(x) {
  mpattern <- "CCCTAA"
  mtable <- NULL
  short_mtable <- NULL
  match_test <- 1 
  
  while(match_test != 0) {
    short_mpattern <- mpattern
    mpattern <- str_c(mpattern, "CCCTAA")
    
    short_mtable <- matchPattern(short_mpattern, x, max.mismatch = 1)
    mtable <- matchPattern(mpattern, x, max.mismatch = 1)
    match_test <- length(mtable)
  }
  return(short_mtable)
})

telomere_3_lt <- lapply(refgenome, function(x) {
  mpattern <- "TTAGGG"
  mtable <- NULL
  short_mtable <- NULL
  match_test <- 1 
  
  while(match_test != 0) {
    short_mpattern <- mpattern
    mpattern <- str_c(mpattern, "TTAGGG")
    
    short_mtable <- matchPattern(short_mpattern, x, max.mismatch = 1)
    mtable <- matchPattern(mpattern, x, max.mismatch = 1)
    match_test <- length(mtable)
  }
  return(short_mtable)
})

get_telo_fct <- function(x) {
  result <- data.frame(start = start(x),
                       end = end(x),
                       width = width(x),
                       contig_length = length(subject(x)))
  start(x)[start(x) < 1] <- 1 
  end(x)[end(x) > length(subject(x))] <- length(subject(x))
  result$seq <- as.character(x)
  return(result)
}

telomere_5_clean_lt <- lapply(telomere_5_lt, get_telo_fct)
telomere_3_clean_lt <- lapply(telomere_3_lt, get_telo_fct)


telomere_3_clean_df <- bind_rows(telomere_3_clean_lt, .id="seqname")
telomere_3_clean_f_df <- filter(telomere_3_clean_df, width >= 24) # 6 *4
colnames(telomere_3_clean_f_df) <- str_c(colnames(telomere_3_clean_f_df), "_telomere_3")

telomere_5_clean_df <- bind_rows(telomere_5_clean_lt, .id="seqname")
telomere_5_clean_f_df <- filter(telomere_5_clean_df, width >= 24) # 6 *4
colnames(telomere_5_clean_f_df) <- str_c(colnames(telomere_5_clean_f_df), "_telomere_5")

telomere_53_clean_f_df <- full_join(telomere_5_clean_f_df, telomere_3_clean_f_df, by=c("seqname_telomere_5"="seqname_telomere_3"))
telomere_53_clean_f_df <- select(telomere_53_clean_f_df, 1,5,10, 2:4, 7:9, 6,11)
telomere_53_clean_f_df <- arrange(telomere_53_clean_f_df, nchar(seqname_telomere_5), seqname_telomere_5)

p1 <- telomere_53_clean_f_df %>% 
  select(seqname_telomere_5, width_telomere_5, width_telomere_3) %>% 
  dplyr::rename(seqname = seqname_telomere_5, telomere_5prim=width_telomere_5, telomere_3prim=width_telomere_3) %>%
  mutate(seqname = factor(seqname, levels = seqname, ordered = T)) %>% 
  gather(key = "telomere_loc", value = "width", -seqname) %>% 
  ggplot(aes(x=seqname, y=width, color=telomere_loc)) +
  geom_point() +
  ylab('telomere length (nt)') +
  xlab('sequence name') +
  ylim(c(0,200)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, family = 'Helvetica'),
        legend.title = element_blank())

ggsave(plot=p1, filename = 'telomereLengthDistribution.pdf', width = 7, height = 4)


write_tsv(telomere_53_clean_f_df, "CopciAB_JGIfull_pacbio_telomere_summary_20220809.tsv")
