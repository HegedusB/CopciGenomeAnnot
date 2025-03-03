# load packages ------------------------------------------------------

suppressPackageStartupMessages({
  library("readr")
  library("dplyr")
  library("stringr")
})


reads_path_lt <- list(rpath1 = "../../Pecs_raw_reads/20200831_rawData_joint/",
                      rpath2 = "../../Pecs_raw_reads/20200914_rawData_joint/",
                      rpath3 = "../../Pecs_raw_reads/20211002_rawData_quantseq_extra/",
                      rpath4 ="../../Pecs_raw_reads/20211019_rawData_quantseq/")

reads_path <- unlist(lapply(reads_path_lt, function(x) list.files(x, pattern = ".fastq.gz$", full.names = TRUE)), use.names = FALSE)

reads_path_plus <- "../../Pecs_raw_reads/20210114_rawData/rawData/"
reads_path_plus <- list.files(reads_path_plus, pattern = ".fastq.gz$", full.names = TRUE)
reads_path_plus <- reads_path_plus[!str_detect(reads_path_plus, pattern = "CB[0-9]|CF[0-9]|MC[0-9]+")]

reads_path_all <- c(reads_path, reads_path_plus)

copci_meta_tbl <- tibble(sample = "",
                         treatment = "",
                         datafile = reads_path_all,
                         comment = "MiSeq")

copci_meta_tbl <- mutate(copci_meta_tbl, sample = sapply(str_split(reads_path_all, pattern="/"), function(x) rev(x)[[1]]))
copci_meta_tbl <- mutate(copci_meta_tbl, sample = str_remove(sample, pattern = ".fastq.gz$"))
copci_meta_tbl <- mutate(copci_meta_tbl, sample = str_remove(sample, pattern = "_R1$"))
copci_meta_tbl <- mutate(copci_meta_tbl, sample = str_remove(sample, pattern = "_S[0-9]+_R1_[0-9]+$"))

# M; MP : Zhihao protoplasted
# C43; C49 : Manis CRE

copci_meta_tbl <- filter(copci_meta_tbl, !str_detect(sample, pattern="^M[0-9]+|^MP-|^C[0-9]+"))
copci_meta_tbl <- mutate(copci_meta_tbl, sample = ifelse(str_detect(datafile, pattern="extra"), str_c(sample, "_extra"), sample))

copci_meta_tbl <- mutate(copci_meta_tbl, treatment = str_remove(sample, pattern = "[0-9]+$"))
copci_meta_tbl <- mutate(copci_meta_tbl, treatment = str_remove(treatment, pattern = "-$|_$"))

copci_meta_tbl[copci_meta_tbl$treatment == "OIDIA","treatment"] <- c("OIDIA_1", "OIDIA_1", "OIDIA_1", "OIDIA_2", "OIDIA_2", "OIDIA_2")

barplot(table(copci_meta_tbl$treatment), las=2)
filter(copci_meta_tbl, treatment %in% names(table(copci_meta_tbl$treatment))[table(copci_meta_tbl$treatment) > 3]) %>% arrange(sample) %>% knitr::kable()

# merge extras
copci_meta_tbl$treatment[str_detect(copci_meta_tbl$sample, pattern='LoKo-[0-9]$')] <- 'LoKo_1'
copci_meta_tbl$treatment[str_detect(copci_meta_tbl$sample, pattern='LoKo.$')] <- 'LoKo_2'
barplot(table(copci_meta_tbl$treatment), las=2)

all_extras <- str_subset(copci_meta_tbl$sample, pattern = 'extra$')
setdiff(str_remove(all_extras, pattern = '_extra'), copci_meta_tbl$sample)

outPath <- '/node9_data/bhegedus/Pecs_raw_reads/rawData_mergedWithExtra/'
all_extrasCat_lt <- lapply(all_extras, function(x) {
  sampleIndex <- str_remove(x, pattern = '_extra')
  catCommand <- str_c('cat ', copci_meta_tbl$datafile[copci_meta_tbl$sample == x], ' ', copci_meta_tbl$datafile[copci_meta_tbl$sample == sampleIndex],
                      ' > ', outPath, sampleIndex, '_merged.fastq.gz')
  result_tbl <- tibble(sample = sampleIndex, sample_extra = x, catCommand = catCommand)
  return(result_tbl)
})
all_extrasCat_tbl <- bind_rows(all_extrasCat_lt)

lapply(all_extrasCat_tbl$catCommand, function(x) system(command = x))

all_extrasCat_tbl$outPath <- sapply(str_split(all_extrasCat_tbl$catCommand, pattern=' '), '[[', 5)

copci_meta_NOextra_tbl <- filter(copci_meta_tbl, !sample %in% all_extrasCat_tbl$sample_extra)
copci_meta_NOextra_tbl$datafile[match(all_extrasCat_tbl$sample, copci_meta_NOextra_tbl$sample)] <- all_extrasCat_tbl$outPath

knitr::kable(copci_meta_NOextra_tbl[match(all_extrasCat_tbl$sample, copci_meta_NOextra_tbl$sample),])
barplot(table(copci_meta_NOextra_tbl$treatment), las=2)
table(sapply(copci_meta_NOextra_tbl$datafile, file.exists))

write_tsv(copci_meta_NOextra_tbl, "copci_meta_all_pecs.tsv")



