# load libraries -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(parallel)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
  library(edgeR)
})


# load data ----------------------------------------------------


## load preprepared interpro annot table ----------------------------------------
ipro_simple_path <- "additionals/CopciABprotT0_iproscan93_header_clean.tsv"
ipro_simple <- read_tsv(ipro_simple_path)


## load sample description table ------------------------------------

sample_info_path <- 'additionals/sample_path_and_grouping_info_plus_new2_noExtra_simple.tsv'
sample_info <- read_tsv(sample_info_path)


## load TFs -------------------------------------------------------------
TFs_tbl <- read_tsv('additionals/CopciAB_TF_oversimply_20230705.tsv')
TFs <- unique(TFs_tbl$Protein_accession)

## load rawcount table ---------------------------------------------
rawReads_RNAseq <- read_tsv("/data/bhegedus/Projects/Copci_webpage/adatok_230703_final/expression_dir/Quantseq/CopciAB_20230724_rawcounts_NOextra.tsv")
rawReads_RNAseq <- column_to_rownames(rawReads_RNAseq, var = "gene_id")
rawReadsMatrix <- as.matrix(rawReads_RNAseq)
colnames(rawReadsMatrix) <- str_replace_all(colnames(rawReadsMatrix), pattern="-", replacement="_")

setdiff(colnames(rawReadsMatrix), sample_info$sample)
setdiff(sample_info$sample, colnames(rawReadsMatrix))

# prepare output dictionary and tags ---------------------------------------
SD<-format(Sys.time(),'%Y%m%d')
out_dir <- str_c("All_comparisons_", SD)
dir.create(out_dir)

fajnev<-"CopciAB"


## Development-related experiments -----------------------------------------------------------
### Basidiospore germination (BS) samples -----------------------------------------------------
BS_lt <- list(ref='BS_0h',
              samples = 'BS_4h,BS_8h,BS_12h')
BS_lt$samples <- unlist(str_split(BS_lt$samples, pattern = ','))

### Oidium germination samples ---------------------------------------------------------------
oidium_lt <- list(ref='OIDIA_0h',
                  samples = 'OIDIA_18h')

### Sclerotium ------------------------------------------------------------------------------
sclero_lt <- list(ref='DM_132h',
                  samples='Sclero')

### Dark Mycelium --------------------------------------------------------------------------
darkM_lt <- list(ref='BS_12h',
                 samples='DM_12h,DM_36h,DM_60h,DM_84h,DM_108h,DM_132h,DM_156h_AE,DM_156h_AM')
darkM_lt$samples <- unlist(str_split(darkM_lt$samples, pattern = ','))

### Early light induction and light induction samples ---------------------------------------

AM_lt <- list(ref="DM_156h_AM",
              samples="ELI1hAM,ELI2hAM,LI2hAM,LI6hAM,LI12hAM,LI18hAM,LI_24h_AM")
AM_lt$samples <- unlist(str_split(AM_lt$samples, pattern=","))

AE_lt <- list(ref="DM_156h_AE",
              samples="ELI1hAE,ELI2hAE,LI2hAE,LI6hAE,LI12hAE,LI18hAE,LI_24h_AE")
AE_lt$samples <- unlist(str_split(AE_lt$samples, pattern=","))

### Light induction 24h (LI 24h) samples --------------------------------------------------

LI24AM_lt <- list(ref='LI18hAM',
                  samples='LI_24h_AM')

LI24AE_lt <- list(ref='LI18hAE',
                  samples='LI_24h_AE')

### Early fruiting in light/dark conditions (LD) samples --------------------------------

EFruit_lt <- list(ref='LI_24h_HK',
                  samples='L_D_6h,L_D_12h,L_D_18h,L_D_24h')
EFruit_lt$samples <- unlist(str_split(EFruit_lt$samples, pattern=","))

### Hyphal knot formation in the dark --------------------------------------------------

HKdark_lt <- list(ref='DM_60h',
                  samples='a_30,a_22_5,a_15,a_7_5,HK,HK_plus_7_5')
HKdark_lt$samples <- unlist(str_split(HKdark_lt$samples, pattern=","))

### Mycelium formation/specification --------------------------------------------------

Dark_mycelium156h_AM_lt <- list(ref="DM_156h_AM",
                                samples="DM_156h_AE")

Early_light_1h_AM_lt <- list(ref="ELI1hAM",
                             samples="ELI1hAE")

Early_light_2h_AM_lt <- list(ref="ELI2hAM",
                             samples="ELI2hAE")

Late_light_2h_pLI_AM_lt <- list(ref="LI2hAM",
                                samples="LI2hAE")

Late_light_6h_pLI_AM_lt <- list(ref="LI6hAM",
                                samples="LI6hAE")

Late_light_12h_pLI_AM_lt <- list(ref="LI12hAM",
                                 samples="LI12hAE")

Late_light_18h_pLI_AM_lt <- list(ref="LI18hAM",
                                 samples="LI18hAE")

Late_light_24h_pLI_AM_lt <- list(ref="LI_24h_AM",
                                 samples="LI_24h_AE")

### Starvation response (SWAT samples) --------------------------------------------------------------

SWAT_lt <- list(ref="DM_60h",
                samples="SWAT2h,SWAT4h,SWAT8h,SWAT16h,SWAT24h")
SWAT_lt$samples <- unlist(str_split(SWAT_lt$samples, pattern=","))

### Starvation ------------------------------------------------------------------------------------

DM_starv_lt <- list(ref="DM_60h",
                    samples="DM-132h,DM-156h-AM")
DM_starv_lt$ref <- unlist(str_split(DM_starv_lt$ref, pattern=","))
DM_starv_lt$samples <- unlist(str_split(DM_starv_lt$samples, pattern=","))

LI_starv_lt <- list(ref="DM_60h",
                    samples="LI6hAM,LI12hAM,LI18hAM,LI-24h-AM")
LI_starv_lt$ref <- unlist(str_split(LI_starv_lt$ref, pattern=","))
LI_starv_lt$samples <- unlist(str_split(LI_starv_lt$samples, pattern=","))

### Cobalt chloride --------------------------------------------------------------------------------

CoCl2_lt <- list(ref="CoK",
                 samples="Cocl")

## Stress conditions, plated cultures --------------------------------------------------------------

Freezing_lt <- list(ref="DM_132h",
                    samples="min20")

Hshock_lt <- list(ref="DM_132h",
                  samples="HS47")

CO2_lt <- list(ref="DM_132h",
               samples="CO2")

Drought_lt <- list(ref="DM_84h",
                   samples="DR")

TrichoMP_lt <- list(ref="DM_60h",
                    samples="CTMP")

TrichoGT_lt <- list(ref="DM_108h",
                    samples="CTGT")

Scrached_lt <- list(ref='DM_84h',
                    samples='Scratchy')

Cold_lt <- list(ref='DM_84h',
                samples='Cold')

HypoxiaSress_lt <- list(ref="HypK",
                        samples="Hyp")


## Stress conditions, liquid cultures ------------------------------------------------------------------

OxiSress_lt <- list(ref="LoKo_Tier2",
                    samples="H2O2")

AcidSress_lt <- list(ref="LoKo_Tier2",
                     samples="pH_4")

AlkalineSress_lt <- list(ref="LoKo_Tier2",
                         samples="pH_9")

OsmoSress_lt <- list(ref="LoKo_Tier2",
                     samples="Sor")

CellWallSress_lt <- list(ref="LoKo_Tier2",
                         samples="CR")

Calci_lt <- list(ref="LoKo_Tier3",
                 samples="Calci")

Coper_lt <- list(ref="LoKo_Tier3",
                 samples="Cu")

Vorico_lt <- list(ref='LoKo_Tier3',
                  samples='Vor')

## deltaCRE (Manish) ------------------------------------------------------------

deltaCRE <- list(ref='DM_84h',
                 samples='C43,C91')
deltaCRE$samples <- unlist(str_split(deltaCRE$samples, pattern=","))

# use all samples ---------------------------------------------------------------------------

all_samples <- unlist(c(BS_lt,
                        oidium_lt,
                        sclero_lt,
                        darkM_lt,
                        AM_lt,AE_lt,
                        LI24AM_lt,LI24AE_lt,
                        EFruit_lt,
                        HKdark_lt,
                        Dark_mycelium156h_AM_lt,Early_light_1h_AM_lt,Early_light_2h_AM_lt,Late_light_2h_pLI_AM_lt,Late_light_6h_pLI_AM_lt,Late_light_12h_pLI_AM_lt,Late_light_18h_pLI_AM_lt,Late_light_24h_pLI_AM_lt,
                        SWAT_lt,
                        DM_starv_lt,LI_starv_lt,
                        CoCl2_lt,
                        Freezing_lt,Hshock_lt,CO2_lt,Drought_lt,TrichoMP_lt,TrichoGT_lt,Scrached_lt,Cold_lt,HypoxiaSress_lt,
                        OxiSress_lt,AcidSress_lt,AlkalineSress_lt,OsmoSress_lt,CellWallSress_lt,Calci_lt,Coper_lt,Vorico_lt,
                        deltaCRE), use.names = F)

all_samples <- str_replace_all(all_samples, pattern="-", replacement="_")
all_samples <- unique(all_samples)

setdiff(all_samples, sample_info$sample_group)
setdiff(sample_info$sample_group, all_samples)

# DGE anapizis -----------------------------------------------

groups <- sample_info$sample_group[match(colnames(rawReadsMatrix), sample_info$sample)]
table(colnames(rawReadsMatrix), groups)
groups <- as.factor(groups)

dge <- DGEList(counts=rawReadsMatrix)
dim(dge)

dge <- calcNormFactors(dge)

cutoff <- 1
drop <- which(apply(cpm(dge), 1, max) < cutoff)
length(drop)
dge <- dge[-drop,] 
dim(dge)

man_palette <- rainbow(max(as.numeric(groups)))
sidecols <- man_palette[as.numeric(groups)]

design <- model.matrix(~0+groups) 
colnames(design) <- gsub("groups","",colnames(design))

pdf(file = paste0(out_dir, "/", fajnev,"_Scale_Voom_Total_gene_reads",SD,".pdf"),width=16,height=12, pointsize = 18, family = "sans", bg = "white")
v=voom(dge,design,plot=TRUE,normalize.method="none") 
dev.off()

pdf(file = paste0(out_dir, "/", fajnev,"_MDS_Scale_Voom_Total_gene",SD,".pdf"),width=16,height=12, pointsize = 18, family = "sans", bg = "white")
plotMDS(v, labels = colnames(v),main=paste0(fajnev,"_Scale_Voom_Total_gene_reads_MDS_plot"),cex=0.75,col=sidecols)
dev.off()



fit <- lmFit(v, design)


comp_list <- list(BS_lt,
                  oidium_lt,
                  sclero_lt,
                  darkM_lt,
                  AM_lt,AE_lt,
                  LI24AM_lt,LI24AE_lt,
                  EFruit_lt,
                  HKdark_lt,
                  Dark_mycelium156h_AM_lt,Early_light_1h_AM_lt,Early_light_2h_AM_lt,Late_light_2h_pLI_AM_lt,Late_light_6h_pLI_AM_lt,Late_light_12h_pLI_AM_lt,Late_light_18h_pLI_AM_lt,Late_light_24h_pLI_AM_lt,
                  SWAT_lt,
                  DM_starv_lt,LI_starv_lt,
                  CoCl2_lt,
                  Freezing_lt,Hshock_lt,CO2_lt,Drought_lt,TrichoMP_lt,TrichoGT_lt,Scrached_lt,Cold_lt,HypoxiaSress_lt,
                  OxiSress_lt,AcidSress_lt,AlkalineSress_lt,OsmoSress_lt,CellWallSress_lt,Calci_lt,Coper_lt,Vorico_lt,
                  deltaCRE)


contrast_vector_lt <- lapply(comp_list, function(y) {
  result <- sapply(str_replace_all(y$samples, pattern = "-", replacement = "_"), function(x) str_c(x, str_replace_all(y$ref, pattern = "-", replacement = "_"), sep = "-"))
  return(result)
})

contrast_vector_c <- unlist(contrast_vector_lt)
contrast.matrix <- makeContrasts(contrasts=contrast_vector_c, levels=design)

fit2= contrasts.fit(fit, contrast.matrix)
fit2= eBayes(fit2)

comparisons <- vector(mode = "list", length = length(contrast_vector_c))
comparisons <- lapply(seq_len(length(contrast_vector_c)), function(x) return(contrast_vector_c[[x]] <- x))
comparisons_names <- sapply(str_split(contrast_vector_c, pattern="-"), function(x) return(str_c(x[[1]], "_vs_", x[[2]])))
comparisons <- setNames(comparisons, comparisons_names)


volcano.plot=function(x){
  cat("Working with ",x,".\n",sep="")
  tt=topTable(fit2, coef=comparisons[[x]], adjust="BH",number=nrow(v))
  basename=x
  write.csv(tt,file=paste(out_dir, "/", basename,".diff_exp.csv",sep=""))
  pdf(width=10,height=10,file=paste(out_dir, "/",basename,".volcano.pdf",sep=""))
  plot(tt[,"logFC"],-log10(tt[,"adj.P.Val"]),ylab="-log10(BH QVal)",xlab="logFC",pch=19,cex=0.30,main=basename)
  abline(h=-log10(0.05),lty=3,col="red")
  abline(v=c(-1,1),lty=3,col="red")
  dev.off()
  signif.upreg=tt[tt[,"logFC"]>=1&tt[,"adj.P.Val"]<=0.05,]
  signif.downreg=tt[tt[,"logFC"]<=-1&tt[,"adj.P.Val"]<=0.05,]
  write.csv(signif.upreg,file=paste(out_dir, "/",basename,".signif.upreg.csv",sep="")) 
  write.csv(signif.downreg,file=paste(out_dir, "/",basename,".signif.downreg.csv",sep="")) 
  signif.change=c(rownames(signif.downreg),rownames(signif.upreg))
  return(signif.change)
}

signif.genes <- unique(unlist(lapply(names(comparisons),FUN=volcano.plot))) 
length(signif.genes)

signif.data <- v[rownames(v)%in%signif.genes,]

# deg analysis summary ----------------------------------

signif_path <- list.files(out_dir, pattern="*.signif.*", full.names=T)
signif_names <- str_remove_all(sapply(str_split(signif_path, pattern="/"), function(x) rev(x)[[1]]), pattern="\\.signif|\\.csv")
signif_lt <- lapply(signif_path, function(x) return(read.csv(x, row.names = 1)))
names(signif_lt) <- signif_names
signif_tbl <- as.data.frame(sapply(signif_lt, nrow)) %>% setNames("sig_n") %>% arrange(sig_n) %>% rownames_to_column(var="comp")
signif_tbl %>%
  mutate(comp = factor(comp, levels = comp, ordered = T)) %>% 
  ggplot(aes(x=comp, y=sig_n)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme(axis.text.x = element_text(angle = 90))

signif_tbl <- mutate(signif_tbl, comp_index = str_remove(comp, pattern = "\\..*"))

comp_lev <- names(comparisons)
comp_lev_index <- c(rep('Basidiospore_germination', 3),
                    rep('Oidium_germination', 1),
                    rep('Sclerotium', 1),
                    rep('Dark_Mycelium', 8),
                    rep("Light_induction",14),
                    rep('Light_induction_24', 2),
                    rep("Fruiting",4),
                    rep('HKformationInDark', 6),
                    rep("Mycelium_specification", 8),
                    rep("Starvation",11), 
                    rep('Cobalt_chloride', 1),
                    rep('StressPlatedCulture', 9),
                    rep('StressLiquidCulture', 8),
                    rep('deltaCRE', 2)) 

comp_lev_tbl <- tibble(comp_lev = comp_lev,
                       comp_lev_index = comp_lev_index)
comp_lev_tbl <- mutate(comp_lev_tbl, index=seq_len(nrow(comp_lev_tbl)))

comp_lev_tbl <- left_join(signif_tbl, comp_lev_tbl, by=c("comp_index"="comp_lev"))
comp_lev_tbl <- arrange(comp_lev_tbl, index)
comp_lev_tbl <- mutate(comp_lev_tbl, sig_type = sapply(str_split(comp, pattern = "\\."), "[[", 2))
comp_lev_tbl <- arrange(comp_lev_tbl, index,desc(sig_type))

write_tsv(comp_lev_tbl, str_c(out_dir,"/Number_of_the_signifgenes.tsv"))

colourCount = length(unique(comp_lev_tbl$comp_lev_index))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


comp_lev_tbl %>%
  mutate(comp = factor(comp, levels = comp, ordered = T)) %>% 
  mutate(comp_lev_index = factor(comp_lev_index, levels = unique(comp_lev_tbl$comp_lev_index), ordered = TRUE)) %>% 
  ggplot(aes(x=comp, y=sig_n, fill=comp_lev_index)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Saple pairs compared") +
  ylab("Number of DEG genes") +
  theme_bw() +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.title = element_text(size = 15),
        legend.title = element_blank())


ggsave(str_c(out_dir,"/Number_of_the_signifgenes.pdf"), device = "pdf", width = 18, height = 9)

## plot only upreg ------------------------------------------

colourCount = length(unique(comp_lev_tbl$comp_lev_index))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

comp_lev_tbl %>%
  filter(sig_type == 'upreg') %>% 
  mutate(comp = str_remove(comp, pattern='\\.upreg')) %>% 
  mutate(comp = factor(comp, levels = comp, ordered = T)) %>% 
  mutate(comp_lev_index = factor(comp_lev_index, levels = unique(comp_lev_tbl$comp_lev_index), ordered = TRUE)) %>% 
  ggplot(aes(x=comp, y=sig_n, fill=comp_lev_index)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Saple pairs compared") +
  ylab("Number of upregulated genes") +
  theme_bw() +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        axis.title = element_text(size = 15),
        legend.title = element_blank())

ggsave(str_c(out_dir,"/Number_of_the_signifgenesOnlyUpreg.pdf"), device = "pdf", width = 18, height = 9)


library(UpSetR)
library(ComplexHeatmap)

signif_genes_lt <- lapply(signif_lt, rownames)
signif_genes_lt <- signif_genes_lt[match(comp_lev_tbl$comp, names(signif_genes_lt))]
comp_lev_lt <- split(signif_genes_lt, comp_lev_tbl$comp_lev_index)

names(comp_lev_lt)
comp <- "Starvation"
comp_lt <- comp_lev_lt[[comp]]
sigtype <- "upreg"
comp_lt <- comp_lt[str_detect(names(comp_lt), pattern=sigtype)]

upreg_genes_mat <- make_comb_mat(comp_lt)
UpSet(upreg_genes_mat)

sapply(comp_lt, function(x) return(any(x == "CopciAB_447627.T0")))


m1 <- upreg_genes_mat
comb_name(m1, readable = TRUE)
comb_name(m1)
comb_size(m1)
comb_degree(m1)
length(extract_comb(m1, "11111111111"))
comb_index <- extract_comb(m1, "11111111111")

# get log2CPM -----------------------------------------------------------------------------------
log2cpm <- as.data.frame(v$E)
normCPM <- 2^log2cpm %>% rownames_to_column(var = "gene_id")
write_tsv(normCPM, str_c(out_dir, "/CopciAB_Quantseq_cpm.tsv"))

normCPM %>% 
  gather(key="samples", value = "cpm", -gene_id) -> cpm_short

cpm_short$sample_group <- sample_info$sample_group[match(cpm_short$samples, sample_info$sample)]

cpm_short2 <- cpm_short %>% 
  group_by(gene_id, sample_group) %>% 
  summarise(mean_cpm = mean(cpm)) %>% 
  ungroup()

setdiff(unlist(str_split(contrast_vector_c, pattern = "-")),unique(cpm_short2$sample_group))

cpm_wide2 <- tidyr::spread(cpm_short2, key = sample_group, value=mean_cpm)
cpm_wide2 <- column_to_rownames(cpm_wide2, var = "gene_id")

write_tsv(rownames_to_column(cpm_wide2, var="gene_id"), str_c(out_dir, "/CopciAB_Quantseq_cpm_pooled.tsv"))

# max mean list -----------------------------------------------------------------------------

max_mean_lt <- lapply(contrast_vector_c, function(x) {
  compairs <- unlist(str_split(x, pattern = "-")) 
  max_mean <- apply(cpm_wide2[compairs],1,max) 
  max_mean <- as.data.frame(max_mean)
  max_mean <- rownames_to_column(max_mean, var = "gene_id")
  return(max_mean)
})
names(max_mean_lt) <- contrast_vector_c

max_mean_tbl <- bind_rows(max_mean_lt, .id="comparisons")
max_mean_tbl$max_mean <- round(max_mean_tbl$max_mean, 3)
max_mean_wide_tbl <- spread(max_mean_tbl, key = "comparisons", value = "max_mean")

signif_m_lt <- signif_lt
setdiff(str_remove(str_replace(names(signif_m_lt), pattern = "_vs_", replacement = "-"), pattern="\\..*$"), max_mean_tbl$comparisons)

for(i in names(signif_m_lt)) {
  deg_tbl <- signif_m_lt[[i]]
  deg_tbl <- rownames_to_column(deg_tbl, var = "gene_id")
  findex <- str_remove(str_replace(i, pattern = "_vs_", replacement = "-"), pattern="\\..*$")
  result <- left_join(deg_tbl, max_mean_wide_tbl[c("gene_id",findex)], by="gene_id")
  colnames(result)[length(colnames(result))] <- "max_meanCPM"
  signif_m_lt[[i]] <- result
}

for(i in names(signif_m_lt)) {
  write_csv(signif_m_lt[[i]], file = str_c(out_dir, "/", i, "_max_meanCPM.csv"))
}


# filter upreg genes (max mean 2 CPM) ------------------------------

signif_mf_lt <- lapply(signif_m_lt, function(x) {
  result <- filter(x, max_meanCPM >=2)
  return(result)
})

table(sapply(signif_mf_lt, nrow) > 0)

for(i in names(signif_mf_lt)) {
  write_csv(signif_mf_lt[[i]], file = str_c(out_dir, "/", i, "_max_meanCPM_th2cpm.csv"))
}

# upreg gene number after filtering ------------------------------------ 

signif_mf_upregn_tbl <- as.data.frame(sapply(signif_mf_lt, nrow)) %>% setNames("sig_n") %>% rownames_to_column(var = "comp_lev")
any(is.na(match(signif_mf_upregn_tbl$comp_lev,comp_lev_tbl$comp)))
signif_mf_upregn_tbl <- left_join(signif_mf_upregn_tbl, comp_lev_tbl[c("comp_index","comp","comp_lev_index","index", "sig_type")], by=c("comp_lev"="comp"))
signif_mf_upregn_tbl <- rename(signif_mf_upregn_tbl, "comp"="comp_lev")
signif_mf_upregn_tbl <- arrange(signif_mf_upregn_tbl, index, desc(sig_type))

write_tsv(signif_mf_upregn_tbl, str_c(out_dir,"/Number_of_the_signifgenes_filtered.tsv"))

# plot 

colourCount = length(unique(signif_mf_upregn_tbl$comp_lev_index))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

signif_mf_upregn_tbl %>%
  mutate(comp = factor(comp, levels = comp, ordered = T)) %>% 
  mutate(comp_lev_index = factor(comp_lev_index, levels = unique(comp_lev_tbl$comp_lev_index), ordered = TRUE)) %>% 
  ggplot(aes(x=comp, y=sig_n, fill=comp_lev_index)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Saple pairs compared") +
  ylab("Number of upregulated genes") +
  theme_bw() +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1), 
        axis.title = element_text(size = 15),
        legend.title = element_blank())


ggsave(str_c(out_dir,"/Number_of_the_signifgenes_filtered.pdf"), device = "pdf", width = 18, height = 9)


## plot only upreg ---------------------------------------------------------

colourCount = length(unique(signif_mf_upregn_tbl$comp_lev_index))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

signif_mf_upregn_tbl %>%
  mutate(comp = ifelse(comp_lev_index == 'Mycelium_specification' & sig_type == 'upreg', str_c(str_remove(comp, pattern='_vs_.*'), '_upreg'), comp)) %>%
  mutate(comp = ifelse(comp_lev_index == 'Mycelium_specification' & sig_type == 'downreg', str_c(str_extract(comp, pattern='(?<=_vs_)(.*)(?=\\.)'), '_upreg'), comp)) %>% 
  mutate(sig_type = ifelse(comp_lev_index == 'Mycelium_specification', 'upreg', sig_type)) %>% 
  filter(sig_type == 'upreg') %>% 
  mutate(comp = str_remove(comp, pattern='\\.upreg')) %>% 
  mutate(comp_lev_index = factor(comp_lev_index, levels = unique(comp_lev_tbl$comp_lev_index), ordered = TRUE)) %>% 
  arrange(comp_lev_index) %>% 
  mutate(comp = factor(comp, levels = comp, ordered = T)) %>% 
  ggplot(aes(x=comp, y=sig_n, fill=comp_lev_index)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Saple pairs compared") +
  ylab("Number of upregulated genes") +
  theme_bw() +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust = 1), 
        axis.title = element_text(size = 15),
        legend.title = element_blank())

ggsave(str_c(out_dir,"/Number_of_the_signifgenes_filteredOnlyUpreg.pdf"), device = "pdf", width = 18, height = 9)

## till plot all devreg --------------------------------------------------------

signif_mf_upregn_tbl %>% 
  filter(!comp_index %in% c('C91_vs_DM_84h', 'C43_vs_DM_84h')) %>% 
  filter(!comp_index %in% c('LI_24h_AM_vs_DM_156h_AM', 'LI_24h_AE_vs_DM_156h_AE')) %>% 
  arrange(desc(index)) %>% 
  mutate(comp_index = factor(comp_index, levels=unique(comp_index), ordered = TRUE)) %>% 
  mutate(sig_type = factor(sig_type, levels=c('upreg', 'downreg'), ordered = TRUE)) %>% 
  ggplot(aes(x=sig_type, y=comp_index, fill=sig_n)) +
  geom_tile(color = "gray50") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = median(signif_mf_upregn_tbl$sig_n), limit = c(min(signif_mf_upregn_tbl$sig_n),max(signif_mf_upregn_tbl$sig_n)), space = "Lab", 
                       name="DEG\nnumber") +
  geom_text(aes(label = sig_n), color = "black", size = 2) +
  xlab('') +
  ylab('Comparisons') +
  theme_minimal()

ggsave(str_c(out_dir,"/Heatmap_Number_of_the_signifgenes_filtered_noLI24h.pdf"), device = "pdf", width = 5, height = 9)

# DEG gene count summary -----------------------------------------------------------------

signif_mf_geneID_lt <- lapply(signif_mf_lt, function(x) {return(select(x, gene_id))}) 
signif_mf_geneID_tbl <- bind_rows(signif_mf_geneID_lt, .id='comps') 
signif_mf_geneID_tbl %>% 
  mutate(deg_type = sapply(str_split(comps, pattern='\\.'), '[[', 2)) %>% 
  mutate(comp_index = str_remove(comps, pattern='\\..*')) -> signif_mf_geneID_tbl

write_tsv(signif_mf_geneID_tbl, str_c(out_dir, '/AllCompSignifDEGWithGeneID.tsv'))


signif_mf_geneID_tbl %>% 
  count(comps) %>% 
  full_join(signif_mf_upregn_tbl, by=c('comps'='comp')) %>% 
  summarise(ident=n==sig_n) %>% 
  .$ident %>% 
  table()

length(unique(sample_info$sample_group))

table(colnames(rawReadsMatrix) %in% sample_info$sample)
table(sample_info$sample %in% colnames(rawReadsMatrix))
  

allSampleInComps <-  unique(unlist(str_split(signif_mf_geneID_tbl$comp_index, pattern = '_vs_')))
length(allSampleInComps) 

table(count(sample_info, sample_group)$n) 
count(sample_info, sample_group) %>% 
  arrange(n)

sample_info %>% 
  filter(!sample_group %in% c('C43', 'C91', 'M', 'MP_1_5h', 'MP_3_5h')) -> sample_info_NOothers
length(unique(sample_info_NOothers$sample_group))


allSampleInComps_NoManish <- allSampleInComps[!allSampleInComps %in% c('C43', 'C91')]
length(allSampleInComps_NoManish) 

sample_info_NOothers %>% 
  filter(!sample_group %in% allSampleInComps_NoManish)

setdiff(allSampleInComps, sample_info$sample_group)


length(unique(signif_mf_geneID_tbl$gene_id))
(length(unique(signif_mf_geneID_tbl$gene_id)) / 13617)*100

signif_mf_geneID_tbl %>% 
  group_by(deg_type) %>% 
  filter(!duplicated(gene_id)) %>% 
  summarise(n=n(),
            p=(n/13617)*100)


# Create Summary tables -----------------------------------------------------------------------

signif_genes_lt <- lapply(signif_mf_lt, "[[", "gene_id")
signif_genes_lt <- signif_genes_lt[match(comp_lev_tbl$comp, names(signif_genes_lt))]
comp_lev_lt <- split(signif_genes_lt, comp_lev_tbl$comp_lev_index)


# Starvation --------------------------------------------------------------------------------------
# SWAT_DEGs_union

names(comp_lev_lt)
comp <- "Starvation"
comp_lt <- comp_lev_lt[[comp]]
sigtype <- "upreg" # upreg, downreg
comp_lt <- comp_lt[str_detect(names(comp_lt), pattern=sigtype)]
comp_lt <- comp_lt[str_detect(names(comp_lt), pattern = "^SWAT")]
names(comp_lt)

knitr::kable(sapply(comp_lt, length))
SWAT_DEGs_union <- unique(unlist(comp_lt,use.names=FALSE))
SWAT_DEGs_union_tbl <- tibble(gene_id = SWAT_DEGs_union)
SWAT_DEGs_union_tbl <- left_join(SWAT_DEGs_union_tbl, ipro_simple, by=c("gene_id"="Protein_accession"))
dim(SWAT_DEGs_union_tbl)

write_tsv(SWAT_DEGs_union_tbl, file = str_c(out_dir, "/SWAT_DEGs_union.tsv"))

filter(SWAT_DEGs_union_tbl, gene_id=="CopciAB_447627.T0") # Roc1
filter(SWAT_DEGs_union_tbl, gene_id=="CopciAB_369621.T0") # WC-1
filter(SWAT_DEGs_union_tbl, gene_id=="CopciAB_8610.T0") # WC-2
filter(SWAT_DEGs_union_tbl, gene_id=="CopciAB_501445.T0") # early light


# DM/LI starvation

names(comp_lev_lt)
comp <- "Starvation"
comp_lt <- comp_lev_lt[[comp]]
sigtype <- "upreg" # upreg, downreg
comp_lt <- comp_lt[str_detect(names(comp_lt), pattern=sigtype)]
comp_lt <- comp_lt[str_detect(names(comp_lt), pattern = "^SWAT", negate=TRUE)]
names(comp_lt)

knitr::kable(sapply(comp_lt, length))

DMstarv <- unique(unlist(comp_lt, use.names = FALSE))

DMstarv_tbl <- tibble(gene_id = DMstarv)
dim(DMstarv_tbl)
DMstarv_tbl <- left_join(DMstarv_tbl, ipro_simple, by=c("gene_id"="Protein_accession"))

write_tsv(DMstarv_tbl, file = str_c(out_dir, "/DM_LI_starvation.tsv"))

filter(DMstarv_tbl, gene_id=="CopciAB_447627.T0") # Roc1
filter(DMstarv_tbl, gene_id=="CopciAB_369621.T0") # WC-1
filter(DMstarv_tbl, gene_id=="CopciAB_8610.T0") # WC-2
filter(DMstarv_tbl, gene_id=="CopciAB_501445.T0") # nob1 early light
dim(DMstarv_tbl)

# all genes starvation genes

all_starv <- unique(c(SWAT_DEGs_union_tbl$gene_id, DMstarv_tbl$gene_id))
length(all_starv)

# Core starvation response gene list

Core_starvation_response_gene_list <- DMstarv_tbl[DMstarv_tbl$gene_id %in% SWAT_DEGs_union_tbl$gene_id,]
setdiff(Core_starvation_response_gene_list$gene_id,intersect(DMstarv_tbl$gene_id,SWAT_DEGs_union_tbl$gene_id))
setdiff(intersect(DMstarv_tbl$gene_id,SWAT_DEGs_union_tbl$gene_id), Core_starvation_response_gene_list$gene_id)
dim(Core_starvation_response_gene_list) 

write_tsv(Core_starvation_response_gene_list, file = str_c(out_dir, "/Core_starvation_response_gene_list.tsv"))

filter(Core_starvation_response_gene_list, gene_id=="CopciAB_447627.T0") # Roc1
filter(Core_starvation_response_gene_list, gene_id=="CopciAB_369621.T0") # WC-1
filter(Core_starvation_response_gene_list, gene_id=="CopciAB_8610.T0") # WC-2
filter(Core_starvation_response_gene_list, gene_id=="CopciAB_501445.T0") # nob1 early light

# Core starvation and TFs

TFs_coreStarv_index <- TFs %in% Core_starvation_response_gene_list$gene_id
table(TFs_coreStarv_index)
TFs_coreStarv <- TFs[TFs_coreStarv_index]

### clean summary table -----------------------------------------------------
SWAT_DEGs_union_tbl$gene_id
DMstarv_tbl$gene_id
Core_starvation_response_gene_list

starvation_induction_summary_tbl <- tibble(comp=c('standard_starvation_induced', 
                                                  'SWAT_starvation_induced', 
                                                  'core_starvation_induced'),
                                      group_n =c(length(DMstarv_tbl$gene_id), 
                                                 length(SWAT_DEGs_union_tbl$gene_id), 
                                                 length(Core_starvation_response_gene_list$gene_id)),
                                      gene_id = c(str_remove_all(str_c(DMstarv_tbl$gene_id, collapse = ','), pattern = '\\.T0'),
                                                  str_remove_all(str_c(SWAT_DEGs_union_tbl$gene_id, collapse = ','), pattern = '\\.T0'),
                                                  str_remove_all(str_c(Core_starvation_response_gene_list$gene_id, collapse = ','), pattern = '\\.T0')))

str_count(starvation_induction_summary_tbl$gene_id, 'CopciAB_')

write_tsv(starvation_induction_summary_tbl, file = str_c(out_dir, "/starvation_induced_cleanSummary.tsv"))



# Light induction ---------------------------------------------------------------------------

# aerial_light_induced_genes

names(comp_lev_lt)
comp <- "Light_induction"
comp_lt <- comp_lev_lt[[comp]]
sigtype <- "upreg" # upreg, downreg
comp_lt <- comp_lt[str_detect(names(comp_lt), pattern=sigtype)]
comp_lt <- comp_lt[str_detect(names(comp_lt), pattern = "AE\\.")]
comp_lt$LI_24h_AE_vs_DM_156h_AE.upreg <- NULL 
names(comp_lt)

knitr::kable(sapply(comp_lt, length))

lAE_DEGs_union <- unique(unlist(comp_lt,use.names=FALSE))
lAE_DEGs_union_tbl <- tibble(gene_id = lAE_DEGs_union)
lAE_DEGs_union_tbl <- left_join(lAE_DEGs_union_tbl, ipro_simple, by=c("gene_id"="Protein_accession"))
dim(lAE_DEGs_union_tbl)

write_tsv(lAE_DEGs_union_tbl, file = str_c(out_dir, "/aerial_light_induced_genes.tsv"))

filter(lAE_DEGs_union_tbl, gene_id=="CopciAB_447627.T0") # Roc1
filter(lAE_DEGs_union_tbl, gene_id=="CopciAB_369621.T0") # WC-1
filter(lAE_DEGs_union_tbl, gene_id=="CopciAB_8610.T0") # WC-2
filter(lAE_DEGs_union_tbl, gene_id=="CopciAB_501445.T0") # nob1 early light

# WithoutCoreStarv
lAE_DEGs_unionNoCoreStarv_tbl <- filter(lAE_DEGs_union_tbl, !gene_id %in% Core_starvation_response_gene_list$gene_id)
dim(lAE_DEGs_unionNoCoreStarv_tbl) 

write_tsv(lAE_DEGs_unionNoCoreStarv_tbl, file = str_c(out_dir, "/aerial_light_induced_genes_WithoutCoreStarv.tsv"))


# attached_light_induced_genes

names(comp_lev_lt)
comp <- "Light_induction"
comp_lt <- comp_lev_lt[[comp]]
sigtype <- "upreg" 
comp_lt <- comp_lt[str_detect(names(comp_lt), pattern=sigtype)]
comp_lt <- comp_lt[str_detect(names(comp_lt), pattern = "AM\\.")]
comp_lt$LI_24h_AM_vs_DM_156h_AM.upreg <- NULL
names(comp_lt)

knitr::kable(sapply(comp_lt, length))

lAM_DEGs_union <- unique(unlist(comp_lt,use.names=FALSE))
lAM_DEGs_union_tbl <- tibble(gene_id = lAM_DEGs_union)
lAM_DEGs_union_tbl <- left_join(lAM_DEGs_union_tbl, ipro_simple, by=c("gene_id"="Protein_accession"))

write_tsv(lAM_DEGs_union_tbl, file = str_c(out_dir, "/attached_light_induced_genes.tsv"))
dim(lAM_DEGs_union_tbl)

filter(lAM_DEGs_union_tbl, gene_id=="CopciAB_447627.T0") # Roc1
filter(lAM_DEGs_union_tbl, gene_id=="CopciAB_369621.T0") # WC-1
filter(lAM_DEGs_union_tbl, gene_id=="CopciAB_8610.T0") # WC-2
filter(lAM_DEGs_union_tbl, gene_id=="CopciAB_501445.T0") # nob1 early light

# WithoutCoreStarv
lAM_DEGs_unionNoCoreStarv_tbl <- filter(lAM_DEGs_union_tbl, !gene_id %in% Core_starvation_response_gene_list$gene_id)
dim(lAM_DEGs_unionNoCoreStarv_tbl) 

write_tsv(lAM_DEGs_unionNoCoreStarv_tbl, file = str_c(out_dir, "/attached_light_induced_genes_WithoutCoreStarv.tsv"))


# light induced intersection
length(intersect(lAE_DEGs_union, lAM_DEGs_union)) 
length(intersect(lAE_DEGs_unionNoCoreStarv_tbl$gene_id, lAM_DEGs_unionNoCoreStarv_tbl$gene_id)) 

# light induced union.

light_union <- union(lAM_DEGs_union, lAE_DEGs_union)
light_union_tbl <- tibble(gene_id = light_union)
light_union_tbl <- left_join(light_union_tbl, ipro_simple, by=c("gene_id"="Protein_accession"))

dim(light_union_tbl) 

write_tsv(light_union_tbl, file = str_c(out_dir, "/light_induced_union.tsv"))

filter(light_union_tbl, gene_id=="CopciAB_447627.T0") # Roc1
filter(light_union_tbl, gene_id=="CopciAB_369621.T0") # WC-1
filter(light_union_tbl, gene_id=="CopciAB_8610.T0") # WC-2
filter(light_union_tbl, gene_id=="CopciAB_501445.T0") # nob1 early light
dim(light_union_tbl)

length(union(lAE_DEGs_unionNoCoreStarv_tbl$gene_id, lAM_DEGs_unionNoCoreStarv_tbl$gene_id))

genes_specific_to_light_induction <- filter(light_union_tbl, gene_id %in% setdiff(light_union_tbl$gene_id, Core_starvation_response_gene_list$gene_id))
write_tsv(genes_specific_to_light_induction, file = str_c(out_dir, "/genes_specific_to_light_induction.tsv"))

filter(genes_specific_to_light_induction, gene_id=="CopciAB_447627.T0") # Roc1
filter(genes_specific_to_light_induction, gene_id=="CopciAB_369621.T0") # WC-1
filter(genes_specific_to_light_induction, gene_id=="CopciAB_8610.T0") # WC-2
filter(genes_specific_to_light_induction, gene_id=="CopciAB_501445.T0") # nob1 early light
dim(genes_specific_to_light_induction)

### clean summary table -----------------------------------------------------


light_induction_summary_tbl <- tibble(comp=c('light_induced_AM', 'light_induced_AE', 'light_induced_AMAEunion', 'core_light_induced'),
                                      group_n =c(length(lAM_DEGs_union), length(lAE_DEGs_union), length(light_union), length(genes_specific_to_light_induction$gene_id)),
                                      gene_id = c(str_remove_all(str_c(lAM_DEGs_union, collapse = ','), pattern = '\\.T0'),
                                                  str_remove_all(str_c(lAE_DEGs_union, collapse = ','), pattern = '\\.T0'),
                                                  str_remove_all(str_c(light_union, collapse = ','), pattern = '\\.T0'),
                                                  str_remove_all(str_c(genes_specific_to_light_induction$gene_id, collapse = ','), pattern = '\\.T0')))

write_tsv(light_induction_summary_tbl, file = str_c(out_dir, "/light_induced_cleanSummary.tsv"))


# Mycelium formation/specification -----------------------------------------------------------

names(comp_lev_lt)
comp <- "Mycelium_specification"
comp_lt <- comp_lev_lt[[comp]]
sigtype <- "upreg" 
comp_lt <- comp_lt[str_detect(names(comp_lt), pattern=sigtype)]

names(comp_lt)
sapply(comp_lt, length) %>% sort() %>% knitr::kable()


mycspec_diffreg_unique_geneList <- unique(unlist(comp_lt, use.names = FALSE))
length(mycspec_diffreg_unique_geneList)

str_remove_all(str_c(mycspec_diffreg_unique_geneList, collapse=','), pattern = '\\.T0')

# upSetR plots

comp_new_lt <- comp_lt
names(comp_new_lt) <- str_remove(names(comp_lt), pattern = '\\..*')
upreg_genes_mat <- make_comb_mat(comp_new_lt)
UpSet(upreg_genes_mat)

# upSetR plots (regroup)
groups <- 'DM_156h_AE_vs_DM_156h_AM.upreg,ELI1hAE_vs_ELI1hAM.upreg,ELI2hAE_vs_ELI2hAM.upreg,LI2hAE_vs_LI2hAM.upreg,LI6hAE_vs_LI6hAM.upreg,LI12hAE_vs_LI12hAM.upreg,LI18hAE_vs_LI18hAM.upreg,LI_24h_AE_vs_LI_24h_AM.upreg'
groups <- str_remove_all(groups, '\\.upreg')

groups <- unlist(str_split(groups, pattern = ','))

ComplexHeatmap::UpSet(upreg_genes_mat, set_order = groups)

m1 <- upreg_genes_mat
comb_name(m1, readable = TRUE)
comb_name(m1)
comb_size(m1)
comb_degree(m1)
length(extract_comb(m1, "00100000"))
comb_index <- extract_comb(m1, "00100000")

length(extract_comb(m1, "11111111"))
comb_index <- extract_comb(m1, "11111111")

View(filter(ipro_simple, Protein_accession %in% comb_index))

m1 <- upreg_genes_mat


index <- sapply(comb_name(m1), function(x) {
  return(sum(as.integer(unlist(str_split(x, pattern = '')))) >=5)
})
index_names <- names(index[index])

multyUpregGlist <- lapply(index_names, function(x) {
  return(extract_comb(m1, x))
})

multyUpregGlist <- unique(unlist(multyUpregGlist))
length(multyUpregGlist)
filter(ipro_simple, Protein_accession %in% multyUpregGlist) %>% View()

str_remove_all(str_c(multyUpregGlist, collapse=','), pattern = '\\.T0')

# create data for stem -------------------------------------------------

# AM

light_AM_stam_index <- str_replace_all(unlist(AM_lt), pattern = "-", replacement = "_")
setdiff(light_AM_stam_index, colnames(cpm_wide2))
light_AM_stam <- cpm_wide2[light_AM_stam_index]
light_AM_stam <- rownames_to_column(light_AM_stam, var="gene_id")
write_tsv(light_AM_stam, str_c(out_dir,"/light_AM_stam.tsv"))

# AE

light_AE_stam_index <- str_replace_all(unlist(AE_lt), pattern = "-", replacement = "_")
setdiff(light_AE_stam_index, colnames(cpm_wide2))
light_AE_stam <- cpm_wide2[light_AE_stam_index]
light_AE_stam <- rownames_to_column(light_AE_stam, var="gene_id")
write_tsv(light_AE_stam, str_c(out_dir,"/light_AE_stam.tsv"))

###############################################################################################################

# Heatmaps -------------------------------------------------

## Hydrophobins ---------------------------------------------------------------
filter(ipro_simple, str_detect(uni_ipro_access, pattern = 'IPR001338'))
hydrophobins <- read_tsv('additionals/Hydrophobins.tsv')

geneset <- str_c(hydrophobins$Hydrophobins, '.T0')
setdiff(geneset, rownames(log2cpm))
geneset_mx <-  log2cpm[rownames(log2cpm) %in% geneset,]

geneset_mx %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id") %>% 
  gather(key = "sample", value = "log2CPM", -gene_id) -> geneset_short_tbl

geneset_short_tbl <- left_join(geneset_short_tbl, sample_info, by="sample")
geneset_short_tbl %>% 
  group_by(sample_group, gene_id) %>% 
  summarise(mean_log2cpm = mean(log2CPM)) %>% 
  ungroup() %>% 
  spread(key = sample_group, value = mean_log2cpm) -> geneset_wide_mean_tbl

geneset_wide_mean_mx <- as.matrix(column_to_rownames(geneset_wide_mean_tbl, var = "gene_id"))

col_order <- 'DM_156h_AM,ELI1hAM,ELI2hAM,LI2hAM,LI6hAM,LI12hAM,LI18hAM,LI_24h_AM,DM_156h_AE,ELI1hAE,ELI2hAE,LI2hAE,LI6hAE,LI12hAE,LI18hAE,LI_24h_AE,LI_24h_HK,L_D_6h,L_D_12h,L_D_18h,L_D_24h'
col_order <- unlist(str_split(col_order, pattern = ','))
col_order <- str_trim(col_order)

geneset_wide_mean_mx_ordered <- geneset_wide_mean_mx[,match(col_order,colnames(geneset_wide_mean_mx))]

geneset_wide_mean_mx_scale_ordered <- t(scale(t(geneset_wide_mean_mx_ordered)))
geneset_wide_mean_mx_scale_ordered <- as.matrix(geneset_wide_mean_mx_scale_ordered)


ht <- Heatmap(geneset_wide_mean_mx_scale_ordered,
              cluster_columns = FALSE, 
              show_row_names=TRUE, show_row_dend = TRUE, show_column_dend = FALSE,
              row_dend_width = unit(4, "cm"), row_dend_gp = gpar(lwd=1.5),
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D",
              heatmap_legend_param = list(title = "Z-score", legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                          labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize=18)), column_names_gp = gpar(fontsize = 20),
              column_names_rot = 60
)

ComplexHeatmap::draw(ht)

pdf(file = paste0(out_dir, "/Heatmap_Hydrophobins_",SD,".pdf"),
    width=12,
    height=8, 
    pointsize = 10, 
    family = "sans", 
    bg = "white")
ComplexHeatmap::draw(ht)
dev.off()

## Core Starvation ------------------------------------------------------------
geneset <- Core_starvation_response_gene_list$gene_id
setdiff(geneset, rownames(log2cpm))
geneset_mx <-  log2cpm[rownames(log2cpm) %in% geneset,]

geneset_short_mx <- geneset_mx[,str_detect(colnames(geneset_mx), pattern = "DM_60|DM_84|DM_132|DM_156h_AM|SWAT|^LI.*AM")]
geneset_short_mx <- geneset_short_mx[,str_detect(colnames(geneset_short_mx), pattern = "LI2", negate = TRUE)]

geneset_short_mx %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id") %>% 
  gather(key = "sample", value = "log2CPM", -gene_id) -> geneset_short_tbl

geneset_short_tbl <- left_join(geneset_short_tbl, sample_info, by="sample")
geneset_short_tbl %>% 
  group_by(sample_group, gene_id) %>% 
  summarise(mean_log2cpm = mean(log2CPM)) %>% 
  ungroup() %>% 
  spread(key = sample_group, value = mean_log2cpm) -> geneset_wide_mean_tbl

geneset_wide_mean_mx <- as.matrix(column_to_rownames(geneset_wide_mean_tbl, var = "gene_id"))
genset_col_order <- c("DM_60h","SWAT2h","SWAT4h","SWAT8h","SWAT16h","SWAT24h", 
                      "DM_84h", "DM_132h","DM_156h_AM", "LI6hAM", "LI12hAM", "LI18hAM", "LI_24h_AM")
match(colnames(geneset_wide_mean_mx), genset_col_order)
geneset_wide_mean_mx <- geneset_wide_mean_mx[,genset_col_order]


geneset_wide_mean_mx_1 <- geneset_wide_mean_mx[,1:6]
geneset_wide_mean_mx_2 <- geneset_wide_mean_mx[,c(1,8:13)]

library(ComplexHeatmap)

geneset_wide_mean_mx_scale_1 <- t(scale(t(geneset_wide_mean_mx_1)))
geneset_wide_mean_mx_scale_2 <- t(scale(t(geneset_wide_mean_mx_2)))

ht1 <- Heatmap(geneset_wide_mean_mx_scale_1,cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = FALSE, 
               heatmap_legend_param = list(title = NULL, legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                           labels_gp = gpar(fontsize = 15)), column_names_gp = gpar(fontsize = 20))
ht2 <-  Heatmap(geneset_wide_mean_mx_scale_2,cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = FALSE, 
                heatmap_legend_param = list(title = NULL, legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                            labels_gp = gpar(fontsize = 15)), column_names_gp = gpar(fontsize = 20))

pdf(file = paste0(out_dir, "/Heatmap_Core_starvation_response_gene_list_SWAT_",SD,".pdf"),width=16,height=12, pointsize = 18, family = "sans", bg = "white")
draw(ht1)
dev.off()

pdf(file = paste0(out_dir, "/Heatmap_Core_starvation_response_gene_list_Standard_",SD,".pdf"),width=16,height=12, pointsize = 18, family = "sans", bg = "white")
draw(ht2)
dev.off()


ht_list <- ht1 + ht2
draw(ht_list)


### cazymes, transporters -------------------------------------------------------------

#### transporters dataset -----------------------------------------
tranporters_list <- read_tsv("additionals/Transporters_IPR.tsv")

interpro_tbl <- read_tsv("additionals/CopciABprotT0_iproscan93_header_clean.tsv")
GO_transmembrane_transport <- filter(interpro_tbl, str_detect(uni_GO_annot, pattern="GO:0055085"))

interpro_tport_index <- sapply(interpro_tbl$uni_ipro_access, function(x) {
  return(any(unlist(str_split(x, pattern = "__")) %in% tranporters_list$IPR_ID))
})
interpro_tport_tbl <- interpro_tbl[interpro_tport_index,]
dim(interpro_tport_tbl) 
write_tsv(interpro_tport_tbl, "additionals/CopciAB_Transporters_IPR.tsv")

length(setdiff(GO_transmembrane_transport$Protein_accession,interpro_tport_tbl$Protein_accession))
length(setdiff(interpro_tport_tbl$Protein_accession, GO_transmembrane_transport$Protein_accession))

#### plot transporters --------------------------------------------------------

#### core + PCWDE + FCW + Transporters v1 --------------------------------------------
cazymes_path <- "additionals/CopciAB_new_dbcan.xlsx"
cazymes_sheets <- excel_sheets(cazymes_path)
cazymes <- read_excel(cazymes_path, cazymes_sheets[1])
cazymes <- dplyr::select(cazymes, 1:6)

cazymes_noNA <- filter(cazymes, !is.na(Category))
table(cazymes_noNA$Category, useNA="ifany")

pcwde <- cazymes_noNA$ProtID[str_detect(cazymes_noNA$Category, pattern = "PCWDE")]
length(pcwde)
fcw <- cazymes_noNA$ProtID[str_detect(cazymes_noNA$Category, pattern = "FCW")]
length(fcw)

cazymes_noNA %>% 
  mutate(swat_starv = as.numeric(ProtID %in% str_remove(SWAT_DEGs_union_tbl$gene_id, pattern="\\.T0$"))) %>% 
  mutate(standard_starv = as.numeric(ProtID %in% str_remove(DMstarv_tbl$gene_id, pattern="\\.T0$"))) %>% 
  mutate(swat_AND_standard_starv = as.numeric(ProtID %in% str_remove(all_starv, pattern="\\.T0$"))) %>% 
  mutate(core_starv = as.numeric(ProtID %in% str_remove(Core_starvation_response_gene_list$gene_id, pattern="\\.T0$"))) -> cazymes_noNA_plus_new

cazymes_noNA_plus_new %>% 
  filter(str_detect(Category, pattern = 'PCWDE')) %>% 
  summarise(all = n(),
            sumPCWDE_swatStarv = sum(swat_starv),
            sumPCWDE_standardStarv = sum(standard_starv),
            sumPCWDE_allStarv = sum(swat_AND_standard_starv),
            sumPCWDE_coreStarv = sum(core_starv))

cazymes_noNA_plus_new %>% 
  filter(str_detect(Category, pattern = 'FCW')) %>% 
  summarise(all = n(),
            sumFCW_swatStarv = sum(swat_starv),
            sumFCW_standardStarv = sum(standard_starv),
            sumFCW_allStarv = sum(swat_AND_standard_starv),
            sumFCW_coreStarv = sum(core_starv))


interpro_tport_tbl %>% 
  mutate(swat_starv = as.numeric(Protein_accession %in% SWAT_DEGs_union_tbl$gene_id),
         standard_starv = as.numeric(Protein_accession %in% DMstarv_tbl$gene_id),
         swat_AND_standard_starv = as.numeric(Protein_accession %in% union(SWAT_DEGs_union_tbl$gene_id, DMstarv_tbl$gene_id)),
         core_starv  = as.numeric(Protein_accession %in% Core_starvation_response_gene_list$gene_id)) -> interpro_tport_tbl_plus

interpro_tport_tbl_plus %>% 
  summarise(all = n(),
            sumTport_swatStarv = sum(swat_starv),
            sumTport_standardStarv = sum(standard_starv),
            sumTport_allStarv = sum(swat_AND_standard_starv),
            sumTport_coreStarv = sum(core_starv))

write_tsv(cazymes_noNA_plus_new, str_c(out_dir, "/CopciAB_Cazymes_pcwdeDistribution.tsv"))
write_tsv(interpro_tport_tbl_plus, str_c(out_dir, "/CopciAB_transportersDistribution.tsv"))

# core substrate 
cazymes_noNA_plus_new %>% 
  filter(str_detect(Category, 'PCWDE')) %>% 
  filter(core_starv == 1) %>% 
  .$Substrate %>%
  str_split(pattern = ',') %>% 
  unlist() %>% 
  str_trim() %>% 
  table()


# swat substrate
cazymes_noNA_plus_new %>% 
  filter(str_detect(Category, 'PCWDE')) %>% 
  filter(swat_starv == 1) %>% 
  .$Substrate %>%
  str_split(pattern = ',') %>% 
  unlist() %>% 
  str_trim() %>% 
  table()


# standard substrate
cazymes_noNA_plus_new %>% 
  filter(str_detect(Category, 'PCWDE')) %>% 
  filter(standard_starv == 1) %>% 
  .$Substrate %>%
  str_split(pattern = ',') %>% 
  unlist() %>% 
  str_trim() %>% 
  table()

## core star gene matrix
geneset_wide_mean_noref_mx <- geneset_wide_mean_mx[,!colnames(geneset_wide_mean_mx) %in% c("DM_84h")]
geneset_wide_mean_noref_mx_scale <- t(scale(t(geneset_wide_mean_noref_mx)))
geneset_wide_mean_noref_mx_scale <- as.matrix(geneset_wide_mean_noref_mx_scale)

## core pcwde gene matrix
pcwde_mx <- matrix(as.numeric(str_remove(rownames(geneset_wide_mean_noref_mx_scale), pattern="\\.T0$") %in% pcwde), ncol=1)
rownames(pcwde_mx) <- rownames(geneset_wide_mean_noref_mx_scale)
colnames(pcwde_mx) <- "PCWDE" 
colSums(pcwde_mx) # 96

# core FCW gene matrix
fcw_mx <- matrix(as.numeric(str_remove(rownames(geneset_wide_mean_noref_mx_scale), pattern="\\.T0$") %in% fcw), ncol=1)
rownames(fcw_mx) <- rownames(geneset_wide_mean_noref_mx_scale)
colnames(fcw_mx) <- "FCW" 
colSums(fcw_mx) # 40

# core transporters gene matrix
intersect(setdiff(GO_transmembrane_transport$Protein_accession,interpro_tport_tbl$Protein_accession),rownames(geneset_wide_mean_noref_mx_scale)) # 1 CopciAB_374148.T0
filter(interpro_tport_tbl, Protein_accession %in% rownames(geneset_wide_mean_noref_mx_scale))

tporters_mx <- matrix(as.numeric(rownames(geneset_wide_mean_noref_mx_scale) %in% interpro_tport_tbl$Protein_accession), ncol=1)
rownames(tporters_mx) <- rownames(geneset_wide_mean_noref_mx_scale)
colnames(tporters_mx) <- "Transporters"
colSums(tporters_mx)

temp_mx <- cbind(pcwde_mx, fcw_mx, tporters_mx)
temp_mx[,2][temp_mx[,2] == 1] <- 2
temp_mx[,3][temp_mx[,3] == 1] <- 3
temp_mx <- as.matrix(apply(temp_mx,1,max))
colnames(temp_mx) <- 'source'

sum(temp_mx[,1]==1)
sum(temp_mx[,1]==2)
sum(temp_mx[,1]==3)

temp_mx[which(pcwde_mx[,1]==1),]

intersect(fcw,pcwde)

# plot heatmap

ht1 <-  Heatmap(geneset_wide_mean_noref_mx_scale,cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE, row_dend_width = unit(2, "cm"), 
                heatmap_legend_param = list(title = NULL, legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                            labels_gp = gpar(fontsize = 15)), column_names_gp = gpar(fontsize = 20),
                column_names_rot = 60)

ht2 <-  Heatmap(pcwde_mx, col=c("white", "red"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                show_heatmap_legend=FALSE,
                column_names_gp = gpar(fontsize = 20, col = "red"), column_names_rot = 60, width = unit(0.8, "cm"))

ht3 <-  Heatmap(fcw_mx, col=c("white", "blue"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                show_heatmap_legend=FALSE,
                column_names_gp = gpar(fontsize = 20, col = "blue"), column_names_rot = 60, width = unit(0.8, "cm"))

ht4 <-  Heatmap(tporters_mx, col=c("white", "darkgreen"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                show_heatmap_legend=FALSE,
                column_names_gp = gpar(fontsize = 20, col = "darkgreen"), column_names_rot = 60, width = unit(0.8, "cm"))


pdf(file = paste0(out_dir, "/Heatmap_Core_starvation_response_gene_list_SWAT_Standard_Withref_pcwde_fcw_transporter_",SD,".pdf"),width=15,height=15, pointsize = 18, family = "sans", bg = "white")
draw(ht1 + ht2 + ht3 + ht4)
dev.off()

#### core + PCWDE + FCW + Transporters v2 -----------------------------------------------
geneset <- Core_starvation_response_gene_list$gene_id
setdiff(geneset, rownames(log2cpm))
geneset_mx <-  log2cpm[rownames(log2cpm) %in% geneset,]

geneset_short_mx <- geneset_mx[,str_detect(colnames(geneset_mx), pattern = "DM_60|DM_84|DM_132|DM_156h_AM|SWAT|^LI.*AM")]
geneset_short_mx <- geneset_short_mx[,str_detect(colnames(geneset_short_mx), pattern = "LI2", negate = TRUE)]

geneset_short_mx %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id") %>% 
  gather(key = "sample", value = "log2CPM", -gene_id) -> geneset_short_tbl

geneset_short_tbl <- left_join(geneset_short_tbl, sample_info, by="sample")
geneset_short_tbl %>% 
  group_by(sample_group, gene_id) %>% 
  summarise(mean_log2cpm = mean(log2CPM)) %>% 
  ungroup() %>% 
  spread(key = sample_group, value = mean_log2cpm) -> geneset_wide_mean_tbl

geneset_wide_mean_mx <- as.matrix(column_to_rownames(geneset_wide_mean_tbl, var = "gene_id"))
genset_col_order <- c("DM_60h","SWAT2h","SWAT4h","SWAT8h","SWAT16h","SWAT24h", 
                      "DM_132h","DM_156h_AM", "LI6hAM", "LI12hAM", "LI18hAM", "LI_24h_AM")
match(colnames(geneset_wide_mean_mx), genset_col_order)
geneset_wide_mean_mx <- geneset_wide_mean_mx[,genset_col_order]

geneset_wide_mean_scale_mx <- as.matrix(t(scale(t(geneset_wide_mean_mx))))



geneset_wide_mean_mx_scale_1 <- geneset_wide_mean_scale_mx[,1:6]
geneset_wide_mean_mx_scale_2 <- geneset_wide_mean_scale_mx[,c(1,7:12)]

set.seed(1234)
kclusR <- geneset_wide_mean_scale_mx %>% dist(method = "euclidean") %>% hclust(method = "ward.D")

# complex heatmap megoldás
library(ComplexHeatmap)

set.seed(1234)

ht1 <- Heatmap(geneset_wide_mean_mx_scale_1,
               cluster_columns = FALSE, cluster_rows = kclusR,
               show_row_names=FALSE, show_row_dend = TRUE,
               column_names_rot = 60,
               heatmap_legend_param = list(title = NULL, legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                           labels_gp = gpar(fontsize = 15)), column_names_gp = gpar(fontsize = 20))

ht2 <-  Heatmap(geneset_wide_mean_mx_scale_2,cluster_columns = FALSE, cluster_rows = kclusR,
                show_row_names=FALSE, show_row_dend = FALSE,
                column_names_rot = 60,
                heatmap_legend_param = list(title = NULL, legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                            labels_gp = gpar(fontsize = 15)), column_names_gp = gpar(fontsize = 20))

ht3 <-  Heatmap(pcwde_mx, col=c("white", "red"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                show_heatmap_legend=FALSE,
                column_names_gp = gpar(fontsize = 20, col = "red"), column_names_rot = 60, width = unit(0.8, "cm"))

ht4 <-  Heatmap(fcw_mx, col=c("white", "darkgreen"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                show_heatmap_legend=FALSE,
                column_names_gp = gpar(fontsize = 20, col = "darkgreen"), column_names_rot = 60, width = unit(0.8, "cm"))

ht5 <-  Heatmap(tporters_mx, col=c("white", "blue"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                show_heatmap_legend=FALSE,
                column_names_gp = gpar(fontsize = 20, col = "blue"), column_names_rot = 60, width = unit(0.8, "cm"))


pdf(file = paste0(out_dir, "/Heatmap_Core_starvation_response_gene_list_SWAT_Standard_Withref_pcwde_fcw_transporter_v2_",SD,".pdf"),width=20,height=8, pointsize = 10, family = "sans", bg = "white")
draw(ht1+ht2+ht3+ht4+ht5)
dev.off()

#### core + PCWDE + FCW + Transporters v3 -----------------------------------------------

rordmax_fct <- function(m,w,b) {
  pcwde_roll_mx <- m[kclusR$order,,drop=FALSE]
  
  temp_roll <- zoo::rollapply(pcwde_roll_mx[,1], width=w, by=b, FUN=sum, align='left', fill=0)
  
  temp <- 0
  temp_v <- vector(mode = 'integer')
  for(i in temp_roll) {
    if(i==0) {
      temp_v <- append(temp_v, temp)
    } else {
      temp <- i
      temp_v <- append(temp_v, temp)
    }
  }
  
  temp_roll <- temp_v
  names(temp_roll) <- clabel
  temp_roll <- temp_roll[match(olabel,names(temp_roll))]
  return(temp_roll)
}

olabel <- kclusR$labels
clabel <- kclusR$labels[kclusR$order]

ht1 <- Heatmap(geneset_wide_mean_mx_scale_1,
               cluster_columns = FALSE, cluster_rows = kclusR,
               show_row_names=FALSE, show_row_dend = TRUE,
               column_names_rot = 60,
               heatmap_legend_param = list(title = NULL, legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                           labels_gp = gpar(fontsize = 15)), column_names_gp = gpar(fontsize = 20))

ht2 <-  Heatmap(geneset_wide_mean_mx_scale_2,cluster_columns = FALSE, cluster_rows = kclusR,
                show_row_names=FALSE, show_row_dend = FALSE,
                column_names_rot = 60,
                heatmap_legend_param = list(title = NULL, legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                            labels_gp = gpar(fontsize = 15)), column_names_gp = gpar(fontsize = 20))

ht3 <-  Heatmap(pcwde_mx, col=c("white", "red"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                show_heatmap_legend=FALSE,
                column_names_gp = gpar(fontsize = 20, col = "red"), column_names_rot = 60, width = unit(0.8, "cm"))


pcwde_mx_reord <- rordmax_fct(m=pcwde_mx,w = 20, b=20)
max(pcwde_mx_reord)

ht3_hist= rowAnnotation(foo = anno_lines(pcwde_mx_reord, ylim = c(0, 13),
                                         width = unit(0.8, "cm"),
                                         axis_param = list(
                                           side = "bottom",
                                           at = c(0, 13),
                                           labels = c("0", "13"),
                                           labels_rot = 45
                                         ),
                                         gp = gpar(col = "black", lty = 2, lwd=1.5))
)

ht4 <-  Heatmap(fcw_mx, col=c("white", "darkgreen"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                show_heatmap_legend=FALSE,
                column_names_gp = gpar(fontsize = 20, col = "darkgreen"), column_names_rot = 60, width = unit(0.8, "cm"))

fcw_mx_reord <- rordmax_fct(m=fcw_mx, w=20, b=20)
max(fcw_mx_reord)

ht4_hist= rowAnnotation(foo = anno_lines(fcw_mx_reord, ylim = c(0, 3),
                                         width = unit(0.8, "cm"),
                                         axis_param = list(
                                           side = "bottom",
                                           at = c(0, 3),
                                           labels = c("0", "3"),
                                           labels_rot = 45
                                         ),
                                         gp = gpar(col = "black", lty = 2, lwd=1.5))
)


ht5 <-  Heatmap(tporters_mx, col=c("white", "blue"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                show_heatmap_legend=FALSE,
                column_names_gp = gpar(fontsize = 20, col = "blue"), column_names_rot = 60, width = unit(0.8, "cm"))

tporters_mx_reord <- rordmax_fct(m=tporters_mx, w=20, b=20)
max(tporters_mx_reord)

ht5_hist= rowAnnotation(foo = anno_lines(tporters_mx_reord, ylim = c(0, 5),
                                         width = unit(0.8, "cm"),
                                         axis_param = list(
                                           side = "bottom",
                                           at = c(0, 5),
                                           labels = c("0", "5"),
                                           labels_rot = 45
                                         ),
                                         gp = gpar(col = "black", lty = 2, lwd=1.5))
)



pdf(file = paste0(out_dir,  "/Heatmap_Core_starvation_response_gene_list_SWAT_Standard_Withref_pcwde_fcw_transporter_v3_",SD,".pdf"),width=10,height=8, pointsize = 10, family = "sans", bg = "white")
draw(ht1+ht2+ht3+ht3_hist+ht4+ht4_hist+ht5+ht5_hist)
dev.off()

### All Starvation induced genes AND marked Core, SWAT only, Standat only genes --------------

SWAT_only <- setdiff(SWAT_DEGs_union_tbl$gene_id, Core_starvation_response_gene_list$gene_id)
DM_only <- setdiff(DMstarv_tbl$gene_id, Core_starvation_response_gene_list$gene_id)
core <- Core_starvation_response_gene_list$gene_id
all_starv <- unique(c(SWAT_DEGs_union_tbl$gene_id, DMstarv_tbl$gene_id))

geneset <- all_starv
geneset_mx <-  log2cpm[rownames(log2cpm) %in% geneset,]

geneset_mx %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id") %>% 
  gather(key = "sample", value = "log2CPM", -gene_id) -> geneset_short_tbl

geneset_short_tbl <- left_join(geneset_short_tbl, sample_info, by="sample")
geneset_short_tbl %>% 
  group_by(sample_group, gene_id) %>% 
  summarise(mean_log2cpm = mean(log2CPM)) %>% 
  ungroup() %>% 
  spread(key = sample_group, value = mean_log2cpm) -> geneset_wide_mean_tbl

geneset_wide_mean_mx <- as.matrix(column_to_rownames(geneset_wide_mean_tbl, var = "gene_id"))

col_order <- 'DM_12h,DM_36h,DM_60h,DM_84h,DM_108h,DM_132h,DM_156h_AE,ELI1hAE,ELI2hAE,LI2hAE,LI6hAE,LI12hAE,LI18hAE,LI_24h_AE,DM_156h_AM,ELI1hAM,ELI2hAM,LI2hAM,LI6hAM,LI12hAM,LI18hAM,LI_24h_AM,LI_24h_HK,L_D_6h,L_D_12h,L_D_18h,L_D_24h'
col_order <- unlist(str_split(col_order, pattern = ','))

geneset_wide_mean_mx_ordered <- geneset_wide_mean_mx[,match(col_order,colnames(geneset_wide_mean_mx))]

geneset_wide_mean_mx_scale_ordered <- t(scale(t(geneset_wide_mean_mx_ordered))) # row scale
geneset_wide_mean_mx_scale_ordered <- as.matrix(geneset_wide_mean_mx_scale_ordered)


# Core starvation
core_mx <- matrix(as.numeric(rownames(geneset_wide_mean_mx_scale_ordered) %in% core), ncol=1)
rownames(core_mx) <- rownames(geneset_wide_mean_mx_scale_ordered)
colnames(core_mx) <- "Core" 
colSums(core_mx)
dim(core_mx)

# SWAT only 
SWAT_only_mx <- matrix(as.numeric(rownames(geneset_wide_mean_mx_scale_ordered) %in% SWAT_only), ncol=1)
rownames(SWAT_only_mx) <- rownames(geneset_wide_mean_mx_scale_ordered)
colnames(SWAT_only_mx) <- "SWAT only" 
colSums(SWAT_only_mx)
dim(SWAT_only_mx)

# Standard only 
DM_only_mx <- matrix(as.numeric(rownames(geneset_wide_mean_mx_scale_ordered) %in% DM_only), ncol=1)
rownames(DM_only_mx) <- rownames(geneset_wide_mean_mx_scale_ordered)
colnames(DM_only_mx) <- "Standard only" 
colSums(DM_only_mx)
dim(DM_only_mx)

# plot 

ht1 <- Heatmap(geneset_wide_mean_mx_scale_ordered,
               cluster_columns = FALSE, 
               show_row_names=FALSE, show_row_dend = TRUE, show_column_dend = FALSE,
               row_dend_width = unit(4, "cm"), row_dend_gp = gpar(lwd=1.5),
               clustering_distance_rows = "euclidean",
               clustering_method_rows = "ward.D",
               heatmap_legend_param = list(title = "Z-score", legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                           labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize=18)), column_names_gp = gpar(fontsize = 20),
               column_names_rot = 60
)

ht2 <-  Heatmap(core_mx, col=c("white", "red"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                heatmap_legend_param = list(title = NULL),
                column_names_gp = gpar(fontsize = 20, col = "red"), column_names_rot = 60)

ht3 <-  Heatmap(SWAT_only_mx, col=c("white", "blue"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                heatmap_legend_param = list(title = NULL),
                column_names_gp = gpar(fontsize = 20, col = "blue"), column_names_rot = 60)

ht4 <-  Heatmap(DM_only_mx, col=c("white", "darkgreen"), cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE,
                heatmap_legend_param = list(title = NULL),
                column_names_gp = gpar(fontsize = 20, col = "darkgreen"), column_names_rot = 60)

pdf(file = paste0(out_dir, "/Heatmap_AllStarvationUpregs_CoreSWATStandardMarks_",SD,".pdf"),width=20,height=10, pointsize = 18, family = "sans", bg = "white")
ComplexHeatmap::draw(ht1+ht2+ht3+ht4)
dev.off()

### All Starvation transporters --------------------------------------------------------------------

### Autophagy ------------------------------------------------------------------------------

autophagy_genes <- read_table('additionals/autophagy_genes.tsv')

autophagy_genes$gene_id[autophagy_genes$gene_id %in% SWAT_DEGs_union_tbl$gene_id] # SWAT ban

geneset <- autophagy_genes$gene_id
setdiff(geneset, rownames(log2cpm))
geneset_mx <-  log2cpm[rownames(log2cpm) %in% geneset,]

geneset_mx %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id") %>% 
  gather(key = "sample", value = "log2CPM", -gene_id) -> geneset_short_tbl

geneset_short_tbl <- left_join(geneset_short_tbl, sample_info, by="sample")
geneset_short_tbl %>% 
  group_by(sample_group, gene_id) %>% 
  summarise(mean_log2cpm = mean(log2CPM)) %>% 
  ungroup() %>% 
  spread(key = sample_group, value = mean_log2cpm) -> geneset_wide_mean_tbl

geneset_wide_mean_mx <- as.matrix(column_to_rownames(geneset_wide_mean_tbl, var = "gene_id"))

col_order <- 'DM_60h,SWAT2h,SWAT4h,SWAT8h,SWAT16h,SWAT24h'
col_order <- unlist(str_split(col_order, pattern = ','))
col_order <- str_trim(col_order)

geneset_wide_mean_mx_ordered <- geneset_wide_mean_mx[,match(col_order,colnames(geneset_wide_mean_mx))]

geneset_wide_mean_mx_scale_ordered <- t(scale(t(geneset_wide_mean_mx_ordered))) # row scale
geneset_wide_mean_mx_scale_ordered <- as.matrix(geneset_wide_mean_mx_scale_ordered)


ht <- Heatmap(geneset_wide_mean_mx_scale_ordered,
              cluster_columns = FALSE, 
              show_row_names=TRUE, show_row_dend = TRUE, show_column_dend = FALSE,
              row_dend_width = unit(4, "cm"), row_dend_gp = gpar(lwd=1.5),
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D",
              heatmap_legend_param = list(title = "Z-score", legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                          labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize=18)), column_names_gp = gpar(fontsize = 20),
              column_names_rot = 60
)

pdf(file = paste0(out_dir, "/Heatmap_autophagySWAT_",SD,".pdf"),width=8,height=5, pointsize = 18, family = "sans", bg = "white")
ComplexHeatmap::draw(ht)
dev.off()

### deltaCRE upregs (Manish) -----------------------------------------------------------------------

deltaCRE_union <- union(comp_lev_lt$deltaCRE$C43_vs_DM_84h.upreg, comp_lev_lt$deltaCRE$C91_vs_DM_84h.upreg)

#### Transporters -------------------------------------------------------------------
temp_t <- interpro_tport_tbl_plus
temp_t$delta_CRE <- as.numeric(temp_t$Protein_accession %in% deltaCRE_union)

t_deltaCRE_summary <- temp_t %>% count(delta_CRE) %>% mutate(p = (n/sum(n))*100)
nt_deltaCRE <- t_deltaCRE_summary[2,2,drop=TRUE]

t_deltaCRE_summary_lt <- list(swat_AND_standard_starv = filter(temp_t, swat_AND_standard_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/nt_deltaCRE)*100),
                              core_starv = filter(temp_t, core_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/nt_deltaCRE)*100),
                              swat_starv = filter(temp_t, swat_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/nt_deltaCRE)*100),
                              standard_starv = filter(temp_t, standard_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/nt_deltaCRE)*100))
t_deltaCRE_summary_tbl <- bind_rows(t_deltaCRE_summary_lt, .id = 'index')

t_deltaCRE_summary_tbl %>% 
  filter(delta_CRE == 1) %>% 
  mutate(all_deltaCRE_n=nt_deltaCRE) %>% 
  select(1,5,3,4) %>% 
  rename('ingroup_deltaCRE_n'='n', 'ingroup_deltaCRE_p' = 'p')

t_deltaCRE_summary_tbl %>% 
  filter(delta_CRE == 1) %>% 
  mutate(index = factor(index, levels = unique(index), ordered = TRUE)) %>% 
  ggplot(aes(x=index, y=n, fill=p)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw()


filter(temp_t, delta_CRE == 1, 
       str_detect(uni_ipro_access, pattern = 'IPR036259'))

#### Cazymes -------------------------------------------------------------------

temp_p <- cazymes_noNA_plus_new
temp_p$delta_CRE <- as.numeric(temp_p$ProtID %in% str_remove(deltaCRE_union, pattern='\\.T0$'))

p_deltaCRE_summary <- temp_p %>% count(delta_CRE) %>% mutate(p = (n/sum(n))*100) 
np_deltaCRE <- p_deltaCRE_summary[2,2,drop=TRUE] 

p_deltaCRE_summary_lt <- list(swat_AND_standard_starv = filter(temp_p, swat_AND_standard_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/np_deltaCRE)*100),
                              core_starv = filter(temp_p, core_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/np_deltaCRE)*100),
                              swat_starv = filter(temp_p, swat_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/np_deltaCRE)*100),
                              standard_starv = filter(temp_p, standard_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/np_deltaCRE)*100))
p_deltaCRE_summary_tbl <- bind_rows(p_deltaCRE_summary_lt, .id = 'index')

p_deltaCRE_summary_tbl %>% 
  filter(delta_CRE == 1) %>% 
  mutate(all_deltaCRE_n=np_deltaCRE) %>% 
  select(1,5,3,4) %>% 
  rename('ingroup_deltaCRE_n'='n', 'ingroup_deltaCRE_p' = 'p')

p_deltaCRE_summary_tbl %>% 
  filter(delta_CRE == 1) %>% 
  mutate(index = factor(index, levels = unique(index), ordered = TRUE)) %>% 
  ggplot(aes(x=index, y=n, fill=p)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw()


#### PCWDE -------------------------------------------------------------------

temp_p <- filter(temp_p, str_detect(Category, pattern = 'PCWDE'))

table(temp_p$delta_CRE)

p_deltaCRE_summary <- temp_p %>% count(delta_CRE) %>% mutate(p = (n/sum(n))*100)
np_deltaCRE <- p_deltaCRE_summary[2,2,drop=TRUE]

p_deltaCRE_summary_lt <- list(swat_AND_standard_starv = filter(temp_p, swat_AND_standard_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/np_deltaCRE)*100),
                              core_starv = filter(temp_p, core_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/np_deltaCRE)*100),
                              swat_starv = filter(temp_p, swat_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/np_deltaCRE)*100),
                              standard_starv = filter(temp_p, standard_starv == 1) %>% count(delta_CRE) %>% mutate(p = (n/np_deltaCRE)*100))
p_deltaCRE_summary_tbl <- bind_rows(p_deltaCRE_summary_lt, .id = 'index')

p_deltaCRE_summary_tbl %>% 
  filter(delta_CRE == 1) %>% 
  mutate(all_deltaCRE_n=np_deltaCRE) %>% 
  select(1,5,3,4) %>% 
  rename('ingroup_deltaCRE_n'='n', 'ingroup_deltaCRE_p' = 'p')

p_deltaCRE_summary_tbl %>% 
  filter(delta_CRE == 1) %>% 
  mutate(index = factor(index, levels = unique(index), ordered = TRUE)) %>% 
  ggplot(aes(x=index, y=n, fill=p)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_bw()


#### TFs -------------------------------------------------------------------------------------
deltaCRE_union 
TFs_coreStarv 
table(TFs_coreStarv_index) 

str_remove_all(str_c(TFs_coreStarv, collapse = ','), pattern = '\\.T0')

TFs_coreStarvDeltaCre_index <- TFs_coreStarv %in% deltaCRE_union
table(TFs_coreStarvDeltaCre_index)

TFs_coreStarvDeltaCre <- TFs_coreStarv[TFs_coreStarvDeltaCre_index]

str_remove_all(str_c(TFs_coreStarvDeltaCre, collapse = ','), pattern = '\\.T0')

## Light induction ------------------------------------------------------------------------

light_induced_genes <- read_csv('additionals/lightinduced_genes.csv')

geneset <- genes_specific_to_light_induction$gene_id
setdiff(geneset, rownames(log2cpm))
geneset_mx <-  log2cpm[rownames(log2cpm) %in% geneset,]

gene_groups <- c("Dark_mycelium156h_AM__DM_156h_AM","Early_light_1h_AM__ELI1hAM","Early_light_2h_AM__ELI2hAM","2h_pLI_AM__LI2hAM","6h_pLI_AM__LI6hAM","12h_pLI_AM__LI12hAM","18h_pLI_AM__LI18hAM")
gene_groups <- str_split(gene_groups, pattern = "__", simplify = TRUE) %>% 
  as.data.frame() %>% 
  setNames(c("name", "sampleAM")) %>% 
  mutate(sampleAE = str_replace(sampleAM, pattern='AM', replacement='AE'))

library(purrr)

gene_groups_lt <- list(AM=gene_groups$sampleAM, AE=gene_groups$sampleAE)
gene_groupsClustering_lt <- list(AM=c("Cluster_1","Cluster_2"), AE=c("Cluster_2", "Cluster_1"))

gene_groupsExtra_lt <- map2(gene_groups_lt, gene_groupsClustering_lt, function(x, y) {
  gene_groups_pattern <- str_c(x, collapse="|")
  geneset_short_mx <- geneset_mx[,str_detect(colnames(geneset_mx), pattern = gene_groups_pattern)]
  
  geneset_short_mx %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene_id") %>% 
    gather(key = "sample", value = "log2CPM", -gene_id) -> geneset_short_tbl
  
  geneset_short_tbl <- left_join(geneset_short_tbl, sample_info, by="sample")
  geneset_short_tbl %>% 
    group_by(sample_group, gene_id) %>% 
    summarise(mean_log2cpm = mean(log2CPM)) %>% 
    ungroup() %>% 
    spread(key = sample_group, value = mean_log2cpm) -> geneset_wide_mean_tbl
  
  geneset_wide_mean_mx <- as.matrix(column_to_rownames(geneset_wide_mean_tbl, var = "gene_id"))
  geneset_wide_mean_mx <- geneset_wide_mean_mx[,match(x, colnames(geneset_wide_mean_mx))]
  
  geneset_wide_mean_mx_scale <- t(scale(t(geneset_wide_mean_mx))) # row scale
  
  geneset_wide_mean_mx_scale <- as.matrix(geneset_wide_mean_mx_scale)
  
  light_induced_genes %>% mutate(namesLong = str_c(str_remove(geneID, pattern='\\.T0'), '  ', names)) -> light_induced_genes
  table(light_induced_genes$geneID %in% rownames(geneset_wide_mean_mx_scale))
  light_induced_genes_f <- light_induced_genes[light_induced_genes$geneID %in% rownames(geneset_wide_mean_mx_scale),]
  ha = rowAnnotation(mark = anno_mark(at=match(light_induced_genes_f$geneID, rownames(geneset_wide_mean_mx_scale)), 
                                      labels = light_induced_genes_f$namesLong, 
                                      which = 'row'))
  
  # clustering
  kclusR <- geneset_wide_mean_mx_scale %>% dist(method = "euclidean") %>% hclust(method = "ward.D") %>% cutree(k=2)
  splitR <- factor(paste0("Cluster_", kclusR), levels=y) 
  
  kclusR_tbl <- as.data.frame(kclusR) %>% setNames('Cluster') %>% rownames_to_column(var = 'gene_id')

  # plot
  ht <- Heatmap(geneset_wide_mean_mx_scale,
                cluster_columns = FALSE, show_row_names=FALSE, show_row_dend = TRUE, 
                row_dend_width = unit(1, "cm"), column_dend_height = unit(1, "cm"), 
                row_dend_gp = gpar(lwd=1.5), column_dend_gp = gpar(lwd=1.5),
                clustering_distance_rows = "euclidean",
                clustering_method_rows = "ward.D",
                border = TRUE, row_gap = unit(2, "mm"),
                row_split = splitR,
                cluster_row_slices = FALSE,
                cluster_column_slices = FALSE,
                heatmap_legend_param = list(title = "Z-score", legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                            labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize=18)), column_names_gp = gpar(fontsize = 20),
                column_names_rot = 60,
                right_annotation = ha
  )
  
  result_lt <- list(mx=geneset_wide_mean_mx_scale, 
                    ha=ha,
                    kclusR_tbl=kclusR_tbl,
                    kclusR=kclusR,
                    splitR=splitR,
                    ht=ht)
  
  return(result_lt)
  
})

table(gene_groupsExtra_lt$AM$kclusR)
table(gene_groupsExtra_lt$AE$kclusR)

AMAE_GeneClusters_tbl <- bind_rows(list(AM=gene_groupsExtra_lt$AM$kclusR_tbl, 
                                        AE=gene_groupsExtra_lt$AE$kclusR_tbl), .id = 'lightInduced') %>% 
  mutate(Cluster = str_c('Cluster_', Cluster))

write_tsv(AMAE_GeneClusters_tbl, 
          str_c(out_dir, "/Heatmap_AMAE_light_induced_genes_clusters_bordered_AllCoreLight_v2_",SD,".tsv"))

AMAE_GeneClusters_tbl %>% count(lightInduced, Cluster)

for(i in names(gene_groupsExtra_lt)) {
  pdf(file = str_c(out_dir, "/Heatmap_attached_light_induced_genes_clusters_bordered_",i,"_AllCoreLight_v2_",SD,".pdf"),
      width=12,
      height=8, 
      pointsize = 10, 
      family = "sans", 
      bg = "white")
  
  ComplexHeatmap::draw(gene_groupsExtra_lt[[i]]$ht)
  dev.off()
}

### Plots for interesting GO groups ---------------------------------------------------------------------- 

GenesByGO_tbl <- read_tsv('All_comparisons_20230928_temp_GSEA/GOlist_DNA_repair_20231010.tsv')
unique(GenesByGO_tbl$geneInGO)

genes_specific_to_light_induction_DNArepair <- filter(genes_specific_to_light_induction, gene_id %in% unique(GenesByGO_tbl$geneInGO))

length(lAE_DEGs_unionNoCoreStarv_tbl$gene_id)
length(lAM_DEGs_unionNoCoreStarv_tbl$gene_id)
table(unique(GenesByGO_tbl$geneInGO) %in% lAE_DEGs_unionNoCoreStarv_tbl$gene_id) 
table(unique(GenesByGO_tbl$geneInGO) %in% lAM_DEGs_unionNoCoreStarv_tbl$gene_id) 

GenesByGO_AEupreg <- unique(GenesByGO_tbl$geneInGO[GenesByGO_tbl$geneInGO %in% lAE_DEGs_unionNoCoreStarv_tbl$gene_id])
GenesByGO_AMupreg <- unique(GenesByGO_tbl$geneInGO[GenesByGO_tbl$geneInGO %in% lAM_DEGs_unionNoCoreStarv_tbl$gene_id])

intersect(GenesByGO_AEupreg, GenesByGO_AMupreg) 
union(GenesByGO_AEupreg, GenesByGO_AMupreg)


geneset <- genes_specific_to_light_induction_DNArepair$gene_id
setdiff(geneset, rownames(log2cpm))
geneset_mx <-  log2cpm[rownames(log2cpm) %in% geneset,]


gene_groups <- c("Dark_mycelium156h_AM__DM_156h_AM","Early_light_1h_AM__ELI1hAM","Early_light_2h_AM__ELI2hAM","2h_pLI_AM__LI2hAM","6h_pLI_AM__LI6hAM","12h_pLI_AM__LI12hAM","18h_pLI_AM__LI18hAM")
gene_groups <- str_split(gene_groups, pattern = "__", simplify = TRUE) %>% 
  as.data.frame() %>% 
  setNames(c("name", "sampleAM")) %>% 
  mutate(sampleAE = str_replace(sampleAM, pattern='AM', replacement='AE')) %>% 
  select(-name) %>% 
  mutate(order = seq_len(n())) %>% 
  pivot_longer(cols = c('sampleAM', 'sampleAE'), names_to = 'samples', values_to = 'group') %>% 
  arrange(samples, order)

gene_groups_pattern <- str_c(gene_groups$group, collapse="|")
geneset_short_mx <- geneset_mx[,str_detect(colnames(geneset_mx), pattern = gene_groups_pattern)]

geneset_short_mx %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "gene_id") %>% 
  gather(key = "sample", value = "log2CPM", -gene_id) -> geneset_short_tbl

geneset_short_tbl <- left_join(geneset_short_tbl, sample_info, by="sample")
geneset_short_tbl %>% 
  group_by(sample_group, gene_id) %>% 
  summarise(mean_log2cpm = mean(log2CPM)) %>% 
  ungroup() %>% 
  spread(key = sample_group, value = mean_log2cpm) -> geneset_wide_mean_tbl

geneset_wide_mean_mx <- as.matrix(column_to_rownames(geneset_wide_mean_tbl, var = "gene_id"))
geneset_wide_mean_mx <- geneset_wide_mean_mx[,match(gene_groups$group, colnames(geneset_wide_mean_mx))]

geneset_wide_mean_mx_scale <- t(scale(t(geneset_wide_mean_mx))) # row scale

geneset_wide_mean_mx_scale <- as.matrix(geneset_wide_mean_mx_scale)

ht <- Heatmap(geneset_wide_mean_mx_scale,
              cluster_columns = FALSE, show_row_names=TRUE, show_row_dend = TRUE, 
              row_dend_width = unit(1, "cm"), column_dend_height = unit(1, "cm"), 
              row_dend_gp = gpar(lwd=1.5), column_dend_gp = gpar(lwd=1.5),
              clustering_distance_rows = "euclidean",
              clustering_method_rows = "ward.D",
              border = TRUE, row_gap = unit(2, "mm"), 
              cluster_row_slices = FALSE,
              cluster_column_slices = FALSE,
              heatmap_legend_param = list(title = "Z-score", legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                          labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize=18)), column_names_gp = gpar(fontsize = 20),
              column_names_rot = 60)

pdf(file = str_c(out_dir, "/Heatmap_attached_light_induced_genes_GO_DNArepair_AllCoreLight_v2_",SD,".pdf"),
    width=18,
    height=10, 
    pointsize = 10, 
    family = "sans", 
    bg = "white")

ComplexHeatmap::draw(ht)
dev.off()


### Fatty acid desaturase__Fatty acid desaturase domain -----------------------------------------
# IPR012171__IPR005804

AMAE_GeneClusters_tbl %>% 
  filter(lightInduced == 'AM') %>% 
  left_join(ipro_simple, by=c('gene_id'='Protein_accession')) %>% 
  filter(str_detect(.$uni_ipro_desc, pattern = 'Fatty acid desaturase')) %>% 
  View()

filter(ipro_simple, str_detect(uni_ipro_desc, pattern = 'Fatty acid desaturase'))

## Mycelium formation/specification ----------------------------------------------------------------------------

# Supplementary table COREAEAM

names(comp_lev_lt)
comp <- "Mycelium_specification"
comp_lt <- comp_lev_lt[[comp]]
sigtype <- c("upreg", "downreg") # upreg, downreg

compUpregDownreg_lt <- lapply(sigtype, function(x) return(comp_lt[str_detect(names(comp_lt), pattern=x)])) %>% setNames(sigtype)

suptab_AEAM_lt <- lapply(compUpregDownreg_lt, function(x) {
  result_lt <- lapply(x, function(y) {
    result_tbl <- tibble(gene_id = y)
    return(result_tbl)
  })
  result <- bind_rows(result_lt, .id='comp')
  return(result)
})

suptab_AEAM_lt$upreg$comp <- str_c(str_remove(suptab_AEAM_lt$upreg$comp, pattern = '_vs_.*'), '_upreg')
suptab_AEAM_lt$downreg$comp <- str_c(str_extract(suptab_AEAM_lt$downreg$comp, pattern = '(?<=_vs_)(.*)(?=\\.downreg)'), '_upreg')

suptab_AEAM_tbl <- bind_rows(suptab_AEAM_lt)

count(suptab_AEAM_tbl,comp) %>% 
  mutate(AE_index = str_detect(comp, pattern='AE')) %>% 
  arrange(desc(AE_index), comp)

suptab_AEAM_tbl %>% 
  mutate(gene_id = str_remove(gene_id, pattern = '\\.T0$')) %>% 
  group_by(comp) %>% 
  summarise(group_n = n(), 
            gene_id = str_c(gene_id, collapse = ',')) %>% 
  ungroup() %>% 
  arrange(desc(str_detect(comp, pattern='AE')), match(str_remove(comp, pattern='h.*'), c('DM_156', 'ELI1', 'ELI2', 'LI2', 'LI6', 'LI12', 'LI18', 'LI_24'))) -> suptab_AEAM_compact_tbl

write_tsv(suptab_AEAM_compact_tbl, paste0(out_dir, "/Supplementary_table_AEAM_",SD,".tsv"))

### Upseter comp Core AE / AM genes ------------------------------------------------
compUpregDownregCombMat_lt <- lapply(compUpregDownreg_lt, make_comb_mat)

UpSet(compUpregDownregCombMat_lt$upreg)

getMultyUpregGeneList_fct <- function(upsetCombM, nUpreg) {

  index <- sapply(comb_name(upsetCombM), function(x) {
    return(sum(as.integer(unlist(str_split(x, pattern = '')))) >=nUpreg)
  })
  index_names <- names(index[index])
  
  multyUpregGlist <- lapply(index_names, function(x) {
    return(extract_comb(upsetCombM, x))
  })
  
  multyUpregGlist <- unique(unlist(multyUpregGlist))
  return(multyUpregGlist)
}

getMultyUpregGeneList <- lapply(compUpregDownregCombMat_lt, 5,FUN = getMultyUpregGeneList_fct)

sapply(getMultyUpregGeneList, length)

getMultyUpregGeneIpro_lt <- lapply(getMultyUpregGeneList, function(x) {
  tibble(gene_id = x) %>% 
    left_join(ipro_simple, by=c('gene_id'='Protein_accession')) -> result_tbl
  return(result_tbl)
})

sapply(getMultyUpregGeneIpro_lt, nrow)

for(i in names(getMultyUpregGeneIpro_lt)) {
  write_tsv(getMultyUpregGeneIpro_lt[[i]], paste0(out_dir, "/Attached_mycelium_specificCoreGenesWithIproAnnot_", i, '_',SD,".tsv"))
}


geneset_wide_mean_mx_scale_list <- lapply(getMultyUpregGeneList, function(x) {

  geneset_mx <-  log2cpm[rownames(log2cpm) %in% x,]
  
  gene_groups_pattern <- str_c('DM_156|ELI|LI[0-9]|LI_24h_A')
  
  geneset_short_mx <- geneset_mx[,str_detect(colnames(geneset_mx), pattern = gene_groups_pattern)]
  
  geneset_short_mx %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "gene_id") %>% 
    gather(key = "sample", value = "log2CPM", -gene_id) -> geneset_short_tbl
  
  geneset_short_tbl <- left_join(geneset_short_tbl, sample_info, by="sample")
  geneset_short_tbl %>% 
    group_by(sample_group, gene_id) %>% 
    summarise(mean_log2cpm = mean(log2CPM)) %>% 
    ungroup() %>% 
    spread(key = sample_group, value = mean_log2cpm) -> geneset_wide_mean_tbl
  
  geneset_wide_mean_mx <- as.matrix(column_to_rownames(geneset_wide_mean_tbl, var = "gene_id"))
  
  geneset_wide_mean_mx_scale <- t(scale(t(geneset_wide_mean_mx))) # row scale
  geneset_wide_mean_mx_scale <- as.matrix(geneset_wide_mean_mx_scale)
  return(geneset_wide_mean_mx_scale)
})

sapply(geneset_wide_mean_mx_scale_list, nrow)


heatmap_lt <- lapply(geneset_wide_mean_mx_scale_list, function(x) {
  set.seed(123)
  
  ht <- Heatmap(x,
                cluster_columns = TRUE, 
                show_row_names=FALSE, show_row_dend = TRUE, 
                row_dend_width = unit(4, "cm"), row_dend_gp = gpar(lwd=1.5), column_dend_gp = gpar(lwd=1.5),
                clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
                clustering_method_rows = "ward.D", clustering_method_columns = "ward.D",
                border = TRUE, row_gap = unit(2, "mm"), 
                heatmap_legend_param = list(title = "Z-score", legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                            labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize=18)), column_names_gp = gpar(fontsize = 20),
                column_names_rot = 60
  )
  return(ht)
})

ComplexHeatmap::draw(heatmap_lt$upreg)
ComplexHeatmap::draw(heatmap_lt$downreg)


# plot into PDF
for(i in names(heatmap_lt)) {
  pdf(file = paste0(out_dir, "/Heatmap_attached_mycelium_specificGenes_", i, '_',SD,".pdf"),
      width=12,
      height=8, 
      pointsize = 10, 
      family = "sans", 
      bg = "white")
  ComplexHeatmap::draw(heatmap_lt[[i]])
  dev.off()
}

heatmapColOrdered_lt <- lapply(geneset_wide_mean_mx_scale_list, function(x) {
  col_order <- 'DM_156h_AE,ELI1hAE,ELI2hAE,LI2hAE,LI6hAE,LI12hAE,LI18hAE,LI_24h_AE,DM_156h_AM,ELI1hAM,ELI2hAM,LI2hAM,LI6hAM,LI12hAM,LI18hAM,LI_24h_AM'
  col_order <- unlist(str_split(col_order, pattern = ','))
  
  x_ordered <- x[,match(col_order,colnames(x))]
  
  set.seed(123)
  
  ht <- Heatmap(x_ordered,
                cluster_columns = FALSE, 
                show_row_names=FALSE, show_row_dend = TRUE, show_column_dend = FALSE,
                row_dend_width = unit(4, "cm"), row_dend_gp = gpar(lwd=1.5),
                clustering_distance_rows = "euclidean",
                clustering_method_rows = "ward.D",
                heatmap_legend_param = list(title = "Z-score", legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                            labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize=18)), column_names_gp = gpar(fontsize = 20),
                column_names_rot = 60
  )
  return(ht)
})

ComplexHeatmap::draw(heatmapColOrdered_lt$upreg)
ComplexHeatmap::draw(heatmapColOrdered_lt$downreg)

# plot into PDF
for(i in names(heatmapColOrdered_lt)) {
  pdf(file = paste0(out_dir, "/Heatmap_attached_mycelium_specificGenes_ColOrdered_", i, '_',SD,".pdf"),
      width=12,
      height=8, 
      pointsize = 10, 
      family = "sans", 
      bg = "white")
  ComplexHeatmap::draw(heatmapColOrdered_lt[[i]])
  dev.off()
}

### core starvation gének eltávolítása --------------------------------------------------


geneset_wide_mean_mx_scale_NOstarv_list <- lapply(geneset_wide_mean_mx_scale_list, function(x) {
  return(x[!rownames(x) %in% Core_starvation_response_gene_list$gene_id,])
})

sapply(geneset_wide_mean_mx_scale_NOstarv_list, nrow)

getMultyUpregGeneIpro_NOstarv_lt <- lapply(getMultyUpregGeneIpro_lt, function(x) {
  return(x[!x$gene_id %in% Core_starvation_response_gene_list$gene_id,])
})

sapply(getMultyUpregGeneIpro_NOstarv_lt, nrow)

for(i in names(getMultyUpregGeneIpro_NOstarv_lt)) {
  write_tsv(getMultyUpregGeneIpro_NOstarv_lt[[i]], paste0(out_dir, "/Attached_mycelium_specificCoreGenesWithIproAnnot_WoutCoreStarvation_", i, '_',SD,".tsv"))
}


#### plot Wout starv -------------------------
heatmapColOrdered_NOstarv_lt <- lapply(geneset_wide_mean_mx_scale_NOstarv_list, function(x) {
  col_order <- 'DM_156h_AE,ELI1hAE,ELI2hAE,LI2hAE,LI6hAE,LI12hAE,LI18hAE,LI_24h_AE,DM_156h_AM,ELI1hAM,ELI2hAM,LI2hAM,LI6hAM,LI12hAM,LI18hAM,LI_24h_AM'
  col_order <- unlist(str_split(col_order, pattern = ','))
  
  x_ordered <- x[,match(col_order,colnames(x))]
  
  set.seed(123)
  
  ht <- Heatmap(x_ordered,
                cluster_columns = FALSE, 
                show_row_names=FALSE, show_row_dend = TRUE, show_column_dend = FALSE,
                row_dend_width = unit(4, "cm"), row_dend_gp = gpar(lwd=1.5),
                clustering_distance_rows = "euclidean",
                clustering_method_rows = "ward.D",
                heatmap_legend_param = list(title = "Z-score", legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                            labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize=18)), column_names_gp = gpar(fontsize = 20),
                column_names_rot = 60
  )
  return(ht)
})

ComplexHeatmap::draw(heatmapColOrdered_NOstarv_lt$upreg)
ComplexHeatmap::draw(heatmapColOrdered_NOstarv_lt$downreg)

# plot into PDF
for(i in names(heatmapColOrdered_NOstarv_lt)) {
  pdf(file = paste0(out_dir, "/Heatmap_attached_mycelium_specificGenes_ColOrdered_WoutCoreStarv", i, '_',SD,".pdf"),
      width=12,
      height=8, 
      pointsize = 10, 
      family = "sans", 
      bg = "white")
  ComplexHeatmap::draw(heatmapColOrdered_NOstarv_lt[[i]])
  dev.off()
}


rAnnots <- bind_rows(tibble(gene_id = TFs, annot=str_c(gene_id, '_TFs'), annot_index = 'TFs'),
                     tibble(gene_id = str_c(hydrophobins$Hydrophobins, '.T0'), annot=str_c(gene_id, '_hydrophobins'), annot_index='hydrophobins'),
                     tibble(gene_id = str_c(fcw, '.T0'), annot=str_c(gene_id, '_fcw'), annot_index='fcw'))

geneComp_lt <- lapply(geneset_wide_mean_mx_scale_NOstarv_list, function(x) {
  filter(rAnnots, gene_id %in% rownames(x)) %>% 
    count(annot_index) -> result_tbl
  return(result_tbl)
})



#### plot Wout starv + with marked genes ----------------------------------------------------

heatmapColOrdered_NOstarvMarked_lt <- lapply(geneset_wide_mean_mx_scale_NOstarv_list, function(x) {
  
  col_order <- 'DM_156h_AE,ELI1hAE,ELI2hAE,LI2hAE,LI6hAE,LI12hAE,LI18hAE,LI_24h_AE,DM_156h_AM,ELI1hAM,ELI2hAM,LI2hAM,LI6hAM,LI12hAM,LI18hAM,LI_24h_AM'
  col_order <- unlist(str_split(col_order, pattern = ','))
  
  x_ordered <- x[,match(col_order,colnames(x))]
  
  set.seed(123)
  
  ha = rowAnnotation(mark = anno_mark(at=match(rAnnots$gene_id, rownames(x)), 
                                      labels = rAnnots$annot, 
                                      which = 'row'))
  
  ht <- Heatmap(x_ordered,
                cluster_columns = FALSE, 
                show_row_names=FALSE, show_row_dend = TRUE, show_column_dend = FALSE,
                row_dend_width = unit(4, "cm"), row_dend_gp = gpar(lwd=1.5),
                clustering_distance_rows = "euclidean",
                clustering_method_rows = "ward.D",
                heatmap_legend_param = list(title = "Z-score", legend_height = unit(6, "cm"), grid_width = unit(1, "cm"),
                                            labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize=18)), column_names_gp = gpar(fontsize = 20),
                column_names_rot = 60,
                right_annotation = ha
  )
  
  return(ht)
})

ComplexHeatmap::draw(heatmapColOrdered_NOstarvMarked_lt$upreg)
ComplexHeatmap::draw(heatmapColOrdered_NOstarvMarked_lt$downreg)

# plot into PDF
for(i in names(heatmapColOrdered_NOstarvMarked_lt)) {
  pdf(file = paste0(out_dir, "/Heatmap_attached_mycelium_specificGenes_ColOrdered_WoutCoreStarv_Marked_", i, '_',SD,".pdf"),
      width=12,
      height=8, 
      pointsize = 10, 
      family = "sans", 
      bg = "white")
  ComplexHeatmap::draw(heatmapColOrdered_NOstarvMarked_lt[[i]])
  dev.off()
}






