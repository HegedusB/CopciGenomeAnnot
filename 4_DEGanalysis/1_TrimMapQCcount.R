# load packages ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(parallel)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

# read raw data -----------------------------------------------------------

# sample description

sinfo_path <- "additionals/sample_info.xlsx"
sinfo_sheets <- excel_sheets(sinfo_path)
sinfo_tsv <- read_excel(sinfo_path, sinfo_sheets[[2]])

sinfo_tsv <- mutate(sinfo_tsv, id_index = str_remove(CPM_ID, pattern = "\\.*\\[.*\\]"))

sinfo_qseq_tsv <- filter(sinfo_tsv, `Sequencing method` == "Quantseq")
sinfo_qseq_tsv <- mutate(sinfo_qseq_tsv, id_index = str_remove_all(id_index, pattern = "\\.|\\_|-"))

# raw read path

reads_path_lt <- list(rpath1 = "../Pecs_raw_reads/20200831_rawData_joint/",
                      rpath2 = "../Pecs_raw_reads/20200914_rawData_joint/",
                      rpath3 = "../Pecs_raw_reads/20211002_rawData_quantseq_extra/",
                      rpath4 ="../Pecs_raw_reads/20211019_rawData_quantseq/")
reads_path <- unlist(lapply(reads_path_lt, function(x) list.files(x, pattern = ".fastq.gz$", full.names = TRUE)), use.names = FALSE)

reads_path_plus <- "../Pecs_raw_reads/20210114_rawData/rawData/"
reads_path_plus <- list.files(reads_path_plus, pattern = ".fastq.gz$", full.names = TRUE)
reads_path_plus <- reads_path_plus[!str_detect(reads_path_plus, pattern = "CB[0-9]|CF[0-9]|MC[0-9]+")]

reads_path_all <- c(reads_path, reads_path_plus)

reads_path_all_tbl <- tibble(path=reads_path_all)

reads_path_all_tbl <- mutate(reads_path_all_tbl, sample_raw = sapply(str_split(path, pattern="/"), function(x) rev(x)[[1]]))
reads_path_all_tbl <- mutate(reads_path_all_tbl, sample_raw = str_remove(sample_raw, pattern = ".fastq.gz$"))
reads_path_all_tbl <- mutate(reads_path_all_tbl, sample = str_remove(sample_raw, pattern = "_R1$"))
reads_path_all_tbl <- mutate(reads_path_all_tbl, sample = str_remove(sample, pattern = "_S[0-9]+_R1_[0-9]+$"))

reads_path_all_tbl <- mutate(reads_path_all_tbl, sample=ifelse(str_detect(path, pattern = "extra"), str_c(sample, "_extra"), sample))

reads_path_all_tbl <- mutate(reads_path_all_tbl, id_index = str_remove_all(sample, pattern = "\\.|\\_|-"))

reads_path_all_tbl <- mutate(reads_path_all_tbl, sample_group = str_remove(sample, pattern = "-[0-9]$|_[0-9]$|_[0-9]_extra|[0-9]$|[0-9]_extra$"))
reads_path_all_tbl <- arrange(reads_path_all_tbl, nchar(sample_group), sample_group)

reads_path_all_tbl$sample_group[str_detect(reads_path_all_tbl$sample, pattern = "OIDIA-[0-9]{2}$")] <- "OIDIA-0h"
reads_path_all_tbl$sample_group[str_detect(reads_path_all_tbl$sample, pattern = "OIDIA-[0-9]{4}$")] <- "OIDIA-18h"

reads_path_all_tbl$sample_group[str_detect(reads_path_all_tbl$sample, pattern = "LoKo[0-9]$")] <- "LoKo_Tier3"
reads_path_all_tbl$id_index[str_detect(reads_path_all_tbl$sample, pattern = "LoKo[0-9]$")] <- str_c(reads_path_all_tbl$id_index[str_detect(reads_path_all_tbl$sample, pattern = "LoKo[0-9]$")], "_Tier3")
reads_path_all_tbl$sample_group[str_detect(reads_path_all_tbl$sample, pattern = "LoKo-[0-9]$")] <- "LoKo_Tier2"
reads_path_all_tbl$id_index[str_detect(reads_path_all_tbl$sample, pattern = "LoKo-[0-9]$")] <- str_c(reads_path_all_tbl$id_index[str_detect(reads_path_all_tbl$sample, pattern = "LoKo-[0-9]$")], "_Tier2")

table(reads_path_all_tbl$sample_group) %>% knitr::kable()
table(table(reads_path_all_tbl$sample_group))

sinfo_qseq_tsv$id_index2 <- sapply(sinfo_qseq_tsv$id_index, function(x) {
  mindex <- reads_path_all_tbl$id_index[unique(which(str_detect(reads_path_all_tbl$id_index, pattern = str_c("^",x, "[0-9]+"))))]
  mindex <- str_c(mindex, collapse = ", ")
})

write_tsv(sinfo_qseq_tsv, "additionals/sample_info_m1.tsv")
write_tsv(reads_path_all_tbl, "additionals/sample_path_and_grouping_info.tsv")

sinfo_qseq_tsv_man <- read_tsv("additionals/sample_info_m1_man.tsv")

sinfo_qseq_tsv_man$id_index2
reads_path_all_tbl$id_index

k1 <- unlist(str_split(sinfo_qseq_tsv_man$id_index2, pattern=", "))

setdiff(k1, reads_path_all_tbl$id_index)
setdiff(reads_path_all_tbl$id_index, k1)

# raw read qc -------------------------------------------------------------

out_dir <- "raw_qc"
dir.create(out_dir)

raw_reads_path_concat <- str_c(reads_path_all_tbl$path, collapse=" ")

qc_command <- str_c("fastqc ", raw_reads_path_concat, " -t 20 -o ", out_dir)
system(command = qc_command)

# multiqc

system(command=str_c("multiqc ", out_dir, " -o ", out_dir, "/"))

# multiqc summary

qc_raw_gstat <- read_tsv("raw_qc/multiqc_data/multiqc_general_stats.txt")
qc_raw_totalrn <- select(qc_raw_gstat, Sample, `FastQC_mqc-generalstats-fastqc-total_sequences`)
colnames(qc_raw_totalrn)[[2]] <- "raw_total_sequences"

hist(qc_raw_totalrn$raw_total_sequences , main="", xlab = "raw_total_sequences")
qc_raw_totalrn_bp <- boxplot(qc_raw_totalrn$raw_total_sequences , main="", ylab = "raw_total_sequences")
qc_raw_totalrn[qc_raw_totalrn$raw_total_sequences %in% qc_raw_totalrn_bp$out,] %>% arrange(raw_total_sequences) %>% knitr::kable() # get outliers

setdiff(qc_raw_totalrn$Sample,reads_path_all_tbl$sample_raw)
setdiff(reads_path_all_tbl$sample_raw,qc_raw_totalrn$Sample)

reads_path_all_p_tbl <- left_join(reads_path_all_tbl, qc_raw_totalrn, by=c("sample_raw"="Sample")) 

# trimming ----------------------------------------------------------------

out_dir <- "trimmed_reads"

reads_path_all_p_tbl <- mutate(reads_path_all_p_tbl, trimmed_reads_path = str_c(out_dir,"/",sample,"_trimmed.fastq.gz"))
reads_path_all_p_lt <- split(reads_path_all_p_tbl, seq_len(nrow(reads_path_all_p_tbl)))

trimmer_fct <- function(x) {
  trimmer_command <- str_c("bbduk.sh in=", x[["path"]], " out=", x[["trimmed_reads_path"]], " ref=additionals/polyA.fa,additionals/TruSeq3-PE-2_plus.fa k=13 ktrim=r useshortkmers=t mink=5 qtrim=r trimq=10 minlength=20")
  system(command=trimmer_command)
}

lapply(reads_path_all_p_lt, FUN=trimmer_fct)


# trimmed reads qc --------------------------------------------------------

out_dir <- "trimmed_qc"
dir.create(out_dir)

trimmed_reads_path_concat <- str_c(reads_path_all_p_tbl$trimmed_reads_path, collapse=" ")

trimmed_qc_command <- str_c("fastqc ", trimmed_reads_path_concat, " -t 40 -o ", out_dir)

# print(trimmed_qc_command)
system(command = trimmed_qc_command)


# multiqc

system(command=str_c("multiqc ", out_dir, " -o ", out_dir, "/"))

# multiqc summary

qc_trimmed_gstat <- read_tsv("trimmed_qc/multiqc_data/multiqc_general_stats.txt")
qc_trimmed_totalrn <- select(qc_trimmed_gstat, Sample, `FastQC_mqc-generalstats-fastqc-total_sequences`)
colnames(qc_trimmed_totalrn)[[2]] <- "trimmed_total_sequences"

summary(qc_trimmed_totalrn$trimmed_total_sequences)

hist(qc_trimmed_totalrn$trimmed_total_sequences , main="", xlab = "trimmed_total_sequences")
qc_trimmed_totalrn_bp <- boxplot(qc_trimmed_totalrn$trimmed_total_sequences , main="", ylab = "trimmed_total_sequences")
qc_trimmed_totalrn[qc_trimmed_totalrn$trimmed_total_sequences %in% qc_trimmed_totalrn_bp$out,] %>% arrange(trimmed_total_sequences) %>% knitr::kable() # get outliers

setdiff(qc_trimmed_totalrn$Sample,reads_path_all_p_tbl$sample)
setdiff(reads_path_all_p_tbl$sample,qc_trimmed_totalrn$Sample)

reads_path_all_p_tbl <- left_join(reads_path_all_p_tbl, qc_trimmed_totalrn, by=c("sample"="Sample")) 
reads_path_all_p_tbl <- mutate(reads_path_all_p_tbl, trimmed_percent_sequences = trimmed_total_sequences/raw_total_sequences)

boxplot(reads_path_all_p_tbl$trimmed_percent_sequences, ylab="trimmed read %")

arrange(reads_path_all_p_tbl,trimmed_percent_sequences) %>% filter(trimmed_percent_sequences<0.95) %>% select(sample_raw,sample,id_index,sample_group,raw_total_sequences,trimmed_total_sequences,trimmed_percent_sequences)


# STAR aln --------------------------------------------------

ref_genome_path <- "ref_genome/"
ref_genome_fasta_path <- "CopciAB_new_jgi_20220113.fasta.gz"
ref_annot_path <- "CopciAB_V2_geneAnnotation.gtf.gz"


overhang <- 50
nCPU <- 30

genom_index.cmd <- str_c("STAR --runThreadN ",nCPU," --runMode genomeGenerate --genomeDir ",ref_genome_path," --genomeFastaFiles ",ref_genome_fasta_path," --sjdbGTFfile ",ref_annot_path," --sjdbOverhang ",overhang, " >ref_genome/STAR_index1.log",sep="")
system(command = genom_index.cmd)

star_outdir <- "STAR_align"
dir.create(star_outdir)

reads_path_all_p1_tbl <- mutate(reads_path_all_p_tbl, star_aln_out = str_replace(trimmed_reads_path, pattern = "^trimmed_reads", replacement = star_outdir))
reads_path_all_p1_tbl <- mutate(reads_path_all_p1_tbl, star_aln_out = str_remove(star_aln_out, pattern = "fastq.gz$"))
reads_path_all_p1_lt <- split(reads_path_all_p1_tbl, seq_len(nrow(reads_path_all_p1_tbl)))

ncpu <- 60

star_alligner_fct <- function(x) {
  star_command <- str_c("STAR --runThreadN ", ncpu, " --genomeDir ", ref_genome_path, " --readFilesIn ", x[["trimmed_reads_path"]],
                        " --readFilesCommand=zcat --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.1 --alignIntronMin 20 --alignIntronMax 4000 --alignMatesGapMax 4000 --outSAMattributes NH HI NM MD --outSAMtype BAM Unsorted --outSAMattrRGline ID:1 LB:1 SM:RNASeq PL:Illumina --outFileNamePrefix ", x[["star_aln_out"]])
  system(command = star_command)
}


lapply(reads_path_all_p1_lt, FUN = star_alligner_fct)

# sort alignment
aln_path <- list.files(star_outdir, pattern = "Aligned.out.bam", full.names = TRUE)

for(i in aln_path) {
  out_tag <- str_extract(i, pattern = ".*_trimmed.Aligned.out")
  command <- str_c("samtools sort -@", ncpu, " ", i, " -o ", out_tag, "_s.bam")
  system(command = command)
}

sort_path <- list.files(star_outdir, pattern = "trimmed.Aligned.out_s.bam", full.names = TRUE)

# index alignment 
lapply(sort_path, function(x) system(command = str_c("samtools index -@", ncpu, " ", x)))

file.remove(aln_path)

# mapping summary ---------------------------------------------------------

aln_sum_path <- list.files(star_outdir, pattern = "final.out", full.names = TRUE)
aln_name <- str_remove_all(aln_sum_path, pattern = "(STAR_align/|_trimmed.Log.final.out)")

aln_sum_lt <- lapply(aln_sum_path, readLines)
aln_sum_lt <- lapply(aln_sum_lt, function(x) return(str_subset(x, pattern = "\\|")))
aln_sum_lt <- lapply(aln_sum_lt, function(x) return(as.data.frame(str_split(x, pattern="\\|\t", simplify=T))))
names(aln_sum_lt) <- aln_name

aln_sum_df <- bind_rows(aln_sum_lt, .id="sample_name")
aln_sum_df <- apply(aln_sum_df,2,function(x) as.character(x))
colnames(aln_sum_df)[2:3] <- c("categ", "value")
aln_sum_df <- data.frame(apply(aln_sum_df,2,function(x) str_trim(x)), stringsAsFactors = FALSE)

aln_sum_df$value <- str_remove(aln_sum_df$value, pattern = "%$")
aln_sum_df$categ <- str_replace_all(aln_sum_df$categ, pattern = " ", replacement = "_")

write_tsv(aln_sum_df, str_c(star_outdir, "/STAR_all_summary.tsv")) 

aln_sum_cat_lt <- split(aln_sum_df, aln_sum_df$categ)

ump <- aln_sum_cat_lt[["Uniquely_mapped_reads_%"]]
ump$value <- as.numeric(ump$value)

omp_boxp <- boxplot(ump$value, main="Uniquely_mapped_reads_%")
ump[ump$value %in% omp_boxp$out,]

write_tsv(ump, str_c(star_outdir, "/STAR_all_summary_Uniquely_mapped_reads_%.tsv"))

umn <- aln_sum_cat_lt[["Uniquely_mapped_reads_number"]]
umn$value <- as.integer(umn$value)

umn_boxp <- boxplot(umn$value, main="Uniquely_mapped_reads_number")
umn[umn$value %in% umn_boxp$out,]

barplot(umn$value, names = umn$sample_name, las=2)
abline(h=1E+7, col = "green")
abline(h=1E+6, col = "red")

summary(umn$value)

write_tsv(umn, str_c(star_outdir, "/STAR_all_summary_Uniquely_mapped_reads_number.tsv"))

aln_sum_df %>% 
  filter(categ %in% c("Uniquely_mapped_reads_number", "Uniquely_mapped_reads_%")) %>% 
  mutate(value = as.numeric(value)) %>% 
  spread(key = categ, value=value) -> aln_sum_short1_df

aln_sum_short1_df %>% knitr::kable()
arrange(aln_sum_short1_df, `Uniquely_mapped_reads_%`) %>% head(n=10) %>% knitr::kable()
arrange(aln_sum_short1_df, `Uniquely_mapped_reads_number`) %>% head(n=10) %>% knitr::kable()

ggplot(aln_sum_short1_df, aes(x=sample_name, y=Uniquely_mapped_reads_number, fill=`Uniquely_mapped_reads_%`)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(yintercept = c(5E6, 1E7), colour=c("red", "green"))

filter(aln_sum_short1_df, Uniquely_mapped_reads_number < 5E6) %>% knitr::kable()

# collect data
setdiff(reads_path_all_p1_tbl$sample, aln_sum_short1_df$sample_name)
setdiff(aln_sum_short1_df$sample_name, reads_path_all_p1_tbl$sample)

reads_path_all_p2_tbl <- left_join(reads_path_all_p1_tbl, aln_sum_short1_df, by=c("sample"="sample_name"))

reads_path_all_p2_tbl %>% 
  mutate(sample_group = factor(sample_group, ordered = T, levels = reads_path_all_p2_tbl %>% group_by(sample_group) %>% summarise(median_umr=median(Uniquely_mapped_reads_number)) %>% arrange(desc(median_umr)) %>% .$sample_group)) %>% 
  ggplot(aes(y=Uniquely_mapped_reads_number, x=sample_group)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

reads_path_all_p2_tbl %>% 
  mutate(sample_group = factor(sample_group, ordered = T, levels = reads_path_all_p2_tbl %>% group_by(sample_group) %>% summarise(median_umr=median(`Uniquely_mapped_reads_%`)) %>% arrange(desc(median_umr)) %>% .$sample_group)) %>% 
  ggplot(aes(y=`Uniquely_mapped_reads_%`, x=sample_group)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))


# feature Counts  ---------------------------------------------------------

# create saf file

library(rtracklayer)
library(Rsubread)
library(tibble)

# exon
ref_annot_path <- "CopciAB_V2_geneAnnotation.gtf.gz"

annot <- import(ref_annot_path)
annot_df <- as.data.frame(annot)
annot_exon_df <- filter(annot_df, type == "exon")
annot_exon_df <- filter(annot_exon_df, str_detect(gene_id, pattern="T0$"))
length(unique(annot_exon_df$gene_id))

saf_exon_file <- dplyr::select(annot_exon_df, gene_id, seqnames, start, end, strand)
colnames(saf_exon_file) <- c("GeneID", "Chr", "Start", "End", "Strand")

saf_exon_path <- str_replace(ref_annot_path, pattern="\\.gtf.gz$", replacement="_exon_T0\\.saf")
write_tsv(saf_exon_file, saf_exon_path)


# get count values 

saf_path <- "CopciAB_V2_geneAnnotation_T0.saf"

sDate <- format(Sys.Date(), "%Y%m%d")
out_dir <- str_c("FCounts_dir_", sDate)
dir.create(out_dir)

strain_name <- "CopciAB"

bam_dir_path <- "STAR_align/"
bams_path <- list.files(bam_dir_path, pattern="_s.bam$", full.names = T)

genomef <- "CopciAB_new_jgi_20220113.fasta"

fcgene<-featureCounts(files=bams_path,genome=genomef,annot.ext=saf_path,isGTFAnnotationFile=F,
                      nthreads=20, isPairedEnd=F, minFragLength=40, maxFragLength=600,
                      countChimericFragments=F,countMultiMappingReads=T, fraction=T, byReadGroup=F, strandSpecific=1)

# save fc result

er <- as.data.frame(fcgene$counts)
er <- rownames_to_column(er, var="gene_id")
er_names <- str_remove(colnames(er), pattern="_trimmed.Aligned.out_s.bam")
colnames(er) <- er_names

write_tsv(er, str_c(out_dir, "/", strain_name, "_", sDate,"_rawcounts.tsv"))
save(er,fcgene,file=str_c(out_dir, "/", strain_name, "_", sDate,"_rawcounts.saved"))

# fc summary

total_alignments_sum <- apply(fcgene$stat[,-1], 2, sum)
names(total_alignments_sum) <- str_remove(names(total_alignments_sum), pattern = "_trimmed.Aligned.out_s.bam")
total_alignments_sum_df <- as.data.frame(t(total_alignments_sum))
total_alignments_sum_df <- mutate(total_alignments_sum_df, Status = "Total_read_n")

stat_plus <- as.data.frame((fcgene$stat[1,-1] / total_alignments_sum) *100)
names(stat_plus) <- str_remove(names(stat_plus), pattern = "_trimmed.Aligned.out_s.bam")
stat_plus <- mutate(stat_plus, Status = "Successfully_assigned_alignments_p")

stat_base <- fcgene$stat
names(stat_base) <- str_remove(names(stat_base), pattern = "_trimmed.Aligned.out_s.bam")
stat_plus <- bind_rows(stat_base, total_alignments_sum_df, stat_plus)

stat_plus <- mutate(stat_plus, Status = str_trim(Status))
write_tsv(stat_plus, str_c(out_dir, "/", strain_name, "_", sDate,"_countstat.tsv"))

saap <- unlist(filter(stat_plus, Status == "Successfully_assigned_alignments_p")[-1])
boxplot(saap, ylab="Successfully assigned alignments %")

summary(saap)

stat_plus %>% 
  gather(2:dim(.)[2], key="sample", value = "read_n") %>% 
  mutate(sample=str_remove(sample, pattern = ".Aligned.out")) %>% 
  filter(Status %in% c("Assigned", "Total_read_n")) %>% 
  mutate(Status = factor(Status, levels = c("Total_read_n", "Assigned"), ordered = TRUE)) %>% 
  ggplot(aes(x=sample, y=read_n, fill=Status)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 5E6, colour="red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

reads_path_all_p3_tbl <- read_tsv("additionals/sample_path_and_grouping_info_plus.tsv")

setdiff(reads_path_all_p3_tbl$sample,names(saap))
setdiff(names(saap), reads_path_all_p3_tbl$sample)

filter(stat_plus, Status %in% c("Total_read_n", "Successfully_assigned_alignments_p")) %>% 
  select(-Status) %>% 
  t(.) %>% 
  as.data.frame() %>% 
  setNames(., c("total_alligned_read_n_new", "saap_new")) %>% 
  tibble::rownames_to_column(var = "sample") -> stat_plus_short

reads_path_all_p3_tbl <- left_join(reads_path_all_p3_tbl, stat_plus_short, by="sample")


stat_plus_short %>% 
  mutate(sample = factor(sample, ordered = TRUE)) %>% 
  ggplot(aes(x=sample, y=total_alligned_read_n_new, fill=saap_new)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 5E6, colour="red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

reads_path_all_p3_tbl %>% 
  mutate(sample_group = factor(sample_group, ordered = T, levels = reads_path_all_p3_tbl %>% group_by(sample_group) %>% summarise(median_totread=median(total_alligned_read_n_new)) %>% arrange(desc(median_totread)) %>% .$sample_group)) %>% 
  ggplot(aes(y=total_alligned_read_n_new, x=sample_group)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

reads_path_all_p3_tbl %>% 
  mutate(sample_group = factor(sample_group, ordered = T, levels = reads_path_all_p3_tbl %>% group_by(sample_group) %>% summarise(median_saap=median(saap_new)) %>% arrange(desc(median_saap)) %>% .$sample_group)) %>% 
  ggplot(aes(y=saap_new, x=sample_group)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))


write_tsv(reads_path_all_p3_tbl, "additionals/sample_path_and_grouping_info_plus_new2.tsv")

reads_path_all_p3_tbl$sample[reads_path_all_p3_tbl$sample == "MP-1_5h1_extra"] <- "MP-1_5h_1_extra"
reads_path_all_p3_tbl$sample[reads_path_all_p3_tbl$sample == "MP-1_5h2_extra"] <- "MP-1_5h_2_extra"

extra_ids <- reads_path_all_p3_tbl$sample[str_detect(reads_path_all_p3_tbl$sample, pattern = "_extra")]
extra_plus_ids <- c(extra_ids, str_remove(extra_ids, pattern = "_extra"))

setdiff(extra_plus_ids, reads_path_all_p3_tbl$sample)
filter(reads_path_all_p3_tbl, sample %in% extra_plus_ids) %>% arrange(sample) %>% View()

er_new <- er
er_new <- dplyr::rename(er_new, "MP-1_5h_1_extra"="MP-1_5h1_extra", "MP-1_5h_2_extra"="MP-1_5h2_extra")
setdiff(extra_plus_ids, names(er_new))

names(er_new)[str_detect(names(er_new), pattern="^MP")]
extra_plus_ids[str_detect(extra_plus_ids, pattern="^MP")]

extra_plus_ids_index <- str_remove(extra_plus_ids, pattern = "_extra$")
extra_plus_ids_lt <- split(extra_plus_ids, extra_plus_ids_index)
for(i in names(extra_plus_ids_lt)) {
  internsum <- rowSums(er_new[extra_plus_ids_lt[[i]]])
  er_new[extra_plus_ids_lt[[i]]] <- NULL
  er_new[i] <- internsum
}

names(er_new)[str_detect(names(er_new), pattern="^M")]

setdiff(names(er_new), names(er))
setdiff(names(er), names(er_new))

write_tsv(er_new, str_c(out_dir, "/", strain_name, "_", sDate,"_rawcounts_NOextra.tsv"))
