# load packages
suppressPackageStartupMessages({
  library(stringr)
  library(readr)
  library(tibble)
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(topGO)
  library(GO.db)
})


out_dir <- "All_comparisons_20230928_temp_GSEA"
dir.create(out_dir)
SD <- format(Sys.time(), "%Y%m%d")

cazymes_path <- "additionals/CopciAB_new_dbcan.xlsx"
cazymes_sheets <- excel_sheets(cazymes_path)
cazymes <- read_excel(cazymes_path, cazymes_sheets[1])
cazymes <- dplyr::select(cazymes, 1:6)

cazymes_noNA <- filter(cazymes, !is.na(Category)) 
table(cazymes_noNA$Category, useNA="ifany")

# GO  ---------------------------------------------------------------------

# Check GO database (version, types)
GO_dbInfo()
keytypes(GO.db)

# load Copci interproscan GO annotation
GO_annot_path <- "additionals/CopciABprotT0_iproscan93_header_clean.tsv"
GO_annot_long <- read_tsv(GO_annot_path)

GO_annot <- dplyr::select(GO_annot_long, Protein_accession, uni_GO_annot)

# clean GO annot

GO_annot_lt <- str_extract_all(GO_annot$uni_GO_annot, pattern = "GO:[0-9]+")
length(GO_annot_lt)
names(GO_annot_lt) <- GO_annot$Protein_accession
GO_annot_clean_lt <- GO_annot_lt[sapply(GO_annot_lt, length) > 0]
GO_annot_clean_lt <- lapply(GO_annot_clean_lt, unique)
length(GO_annot_clean_lt) 
GO_annot_clean_lt <- lapply(GO_annot_clean_lt, function(x) return(tibble(GO_id = x)))
GO_annot_clean_tbl <- bind_rows(GO_annot_clean_lt, .id = "transcript_id")

GO_O <- select(GO.db, keys=unique(GO_annot_clean_tbl$GO_id), columns=c("TERM","ONTOLOGY"), keytype="GOID")
sum(is.na(GO_O$ONTOLOGY)) 

GO_annot_cleanO_tbl <- left_join(GO_annot_clean_tbl, GO_O, by=c("GO_id"="GOID"))
sum(is.na(GO_annot_cleanO_tbl$ONTOLOGY)) 

write_tsv(GO_annot_cleanO_tbl, str_c(out_dir, "/GO_annot_cleanO.tsv"))

GO_annot_clean_short_tbl <- GO_annot_cleanO_tbl %>% 
  group_by(transcript_id) %>% 
  summarise(GOID = str_c(GO_id, collapse = ","))

write_tsv(GO_annot_clean_short_tbl, str_c(out_dir, "/GO_annot_cleanO_short.tsv"))

# topGO futtatása ---------------------------------------------------------------------

## load upreg geneset -------------------------------------------------

all_signif_path <- list.files("All_comparisons_20230928/", pattern = ".*._max_meanCPM_th2cpm.csv$", full.names = TRUE)

### Starvation ---------------------------------------------------------------------------------------

# (SWAT)
all_SWAT_upreg_path <- str_subset(all_signif_path, pattern = "SWAT.*upreg_max_meanCPM_th2cpm.csv$")
all_SWAT_upreg_lt <- lapply(all_SWAT_upreg_path, function(x) {
  return(read_csv(x)$gene_id)
})
names(all_SWAT_upreg_lt) <- str_remove(sapply(str_split(all_SWAT_upreg_path, pattern = "/"), function(x) rev(x)[1]), pattern = "\\.upreg_max_meanCPM_th2cpm.csv$")

sapply(all_SWAT_upreg_lt, length)

# (DM)
all_DM_upreg_path <- str_subset(all_signif_path, pattern = "DM.*DM_60h.*upreg_max_meanCPM_th2cpm.csv$")
all_DM_upreg_lt <- lapply(all_DM_upreg_path, function(x) {
  return(read_csv(x)$gene_id)
})
names(all_DM_upreg_lt) <- str_remove(sapply(str_split(all_DM_upreg_path, pattern = "/"), function(x) rev(x)[1]), pattern = "\\.upreg_max_meanCPM_th2cpm.csv$")

sapply(all_DM_upreg_lt, length)

# (LI)
all_LIstarv_upreg_path <- str_subset(all_signif_path, pattern = "LI.*DM_60h.*upreg_max_meanCPM_th2cpm.csv$")
all_LIstarv_upreg_lt <- lapply(all_LIstarv_upreg_path, function(x) {
  return(read_csv(x)$gene_id)
})
names(all_LIstarv_upreg_lt) <- str_remove(sapply(str_split(all_LIstarv_upreg_path, pattern = "/"), function(x) rev(x)[1]), pattern = "\\.upreg_max_meanCPM_th2cpm.csv$")

sapply(all_LIstarv_upreg_lt, length)

### Light induction ---------------------------------------------------------------------------------
# Light induction AM

all_AM_upreg_path <- str_subset(all_signif_path, pattern = ".*AM_vs_DM_156h_AM.upreg_max_meanCPM_th2cpm.csv$")
all_AM_upreg_lt <- lapply(all_AM_upreg_path, function(x) {
  return(read_csv(x)$gene_id)
})
names(all_AM_upreg_lt) <- str_remove(sapply(str_split(all_AM_upreg_path, pattern = "/"), function(x) rev(x)[1]), pattern = "\\.upreg_max_meanCPM_th2cpm.csv$")

sapply(all_AM_upreg_lt, length)
all_AM_upreg_lt$LI_24h_AM_vs_DM_156h_AM <- NULL 

# Light induction AE
all_AE_upreg_path <- str_subset(all_signif_path, pattern = ".*AE_vs_DM_156h_AE.upreg_max_meanCPM_th2cpm.csv$")
all_AE_upreg_lt <- lapply(all_AE_upreg_path, function(x) {
  return(read_csv(x)$gene_id)
})
names(all_AE_upreg_lt) <- str_remove(sapply(str_split(all_AE_upreg_path, pattern = "/"), function(x) rev(x)[1]), pattern = "\\.upreg_max_meanCPM_th2cpm.csv$")

sapply(all_AE_upreg_lt, length)
all_AE_upreg_lt$LI_24h_AE_vs_DM_156h_AE <- NULL 

### Mycelium formation/specification -----------------------------------------------------------------------------

# AE upreg (comp upreg)

all_micSpec_upreg_path <- str_subset(all_signif_path, pattern = ".*AE_vs_.*AM\\.upreg.*")
all_micSpec_upreg_lt <- lapply(all_micSpec_upreg_path, function(x) {
  return(read_csv(x)$gene_id)
})
names(all_micSpec_upreg_lt) <- str_remove(sapply(str_split(all_micSpec_upreg_path, pattern = "/"), function(x) rev(x)[1]), pattern = "\\.upreg_max_meanCPM_th2cpm.csv$")

sapply(all_micSpec_upreg_lt, length)

# AM upreg downreg (comp downreg) 

all_micSpec_downreg_path <- str_subset(all_signif_path, pattern = ".*AE_vs_.*AM\\.downreg.*")
all_micSpec_downreg_lt <- lapply(all_micSpec_downreg_path, function(x) {
  return(read_csv(x)$gene_id)
})
names(all_micSpec_downreg_lt) <- str_remove(sapply(str_split(all_micSpec_downreg_path, pattern = "/"), function(x) rev(x)[1]), pattern = "\\.downreg_max_meanCPM_th2cpm.csv$")

sapply(all_micSpec_downreg_lt, length)


## load upreg geneset  ----------------------
AllCompSignifDEGWithGeneID <- read_tsv('All_comparisons_20230928/AllCompSignifDEGWithGeneID.tsv')


## create topGO geneUniverse --------------------------------------------------------------------------------
geneID2GO <- topGO::readMappings(str_c(out_dir,"/GO_annot_cleanO_short.tsv"))
geneUniverse <- names(geneID2GO)

options(scipen=0)

# create GO function ----------------------------------------------------------------------------------

GO_enrichment_fct <- function(x) {
  genesOfInterest <- x
  
  geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
  names(geneList) <- geneUniverse
  allRes_lt <- lapply(c("BP", "MF", "CC"), function(z) {
    myGOdata <- new("topGOdata", description="My project", ontology=z, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize =10) # nodeSize =10
    resultW01Fisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
    topNodes <- sum(resultW01Fisher@score <= 0.05)
    if(topNodes == 1) {
      allRes <- GenTable(myGOdata, W01Fisher = resultW01Fisher, orderBy = "resultW01Fisher", ranksOf = "W01Fisher", topNodes = 2) # berakhatok olyan termet amely nem éri el a p<=0.05 ös küszöbértéket ezért késöbb ismételten szürnöm kell
    } else {
      allRes <- GenTable(myGOdata, W01Fisher = resultW01Fisher, orderBy = "resultW01Fisher", ranksOf = "W01Fisher", topNodes = topNodes)
    }
    
    GO_cazymes_p <- sapply(allRes$GO.ID, function(y) {
      geneInGO <- names(myGOdata@graph@nodeData@data[[y]]$genes)
      sg <- sigGenes(myGOdata)
      geneInGO_sg <- intersect(geneInGO, sg)
      cazymes_short <- cazymes_noNA$ProtID[str_detect(cazymes_noNA$Category, pattern = "PCWDE")]
      result <- mean(str_remove(geneInGO_sg, pattern = "\\.T0$") %in% cazymes_short)
      return(result)
    })
    allRes$Cazymes_PCWDE_p <- GO_cazymes_p
    allRes$Ontology <- z
    return(allRes)
  })
  allRes_tbl <- bind_rows(allRes_lt)
  allRes_tbl$W01Fisher <- as.numeric(str_replace(allRes_tbl$W01Fisher, pattern = "<.*", replacement = "0"))
  allRes_tbl <- filter(allRes_tbl, W01Fisher <= 0.05)
  return(allRes_tbl)
}

# GO Starvation ------------------------------------------------------------------------------------------------
## all starving samples separately -----------------------------------------------------------

all_starvation_lt <- c(all_SWAT_upreg_lt, all_DM_upreg_lt, all_LIstarv_upreg_lt)
GO_enrichment_lt <- lapply(all_starvation_lt, function(x) GO_enrichment_fct(x))

# plot
GO_enrichment_tbl <- bind_rows(GO_enrichment_lt, .id = "comp_name")

GO_enrichment_tbl %>% 
  mutate(comp_name = str_remove(comp_name, pattern = "\\.signif\\.upreg$")) %>%
  mutate(comp_name = str_remove(comp_name, pattern = "_vs_.*")) %>%
  mutate(comp_name = factor(comp_name, levels = c("SWAT2h", "SWAT4h", "SWAT8h", "SWAT16h", "SWAT24h",
                                                  "DM_132h", "DM_156h_AM",
                                                  "LI6hAM", "LI12hAM", "LI18hAM", "LI_24h_AM"), ordered = T)) %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(W01Fisher <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>% 
  dplyr::rename("pvalue"="W01Fisher", "GeneRatio"="signif_p", "PCWDE_Ratio"="Cazymes_PCWDE_p") %>% 
  ggplot(aes(y=Term, x= GeneRatio, size=PCWDE_Ratio, color=pvalue)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle=45, hjust = 1))

ggsave(str_c(out_dir,"/GOenrichment_starvation_", SD, ".pdf"), device = "pdf", width = 20, height = 15)

##  Core Starvation genes ------------------------------------------------------------------------------
coreStarvGenes <- read_tsv("All_comparisons_20230928/Core_starvation_response_gene_list.tsv")
coreStarvGenes <- coreStarvGenes$gene_id
coreStarvGenes_lt <- list(coreStarvGenes=coreStarvGenes)

GO_enrichment_lt <- lapply(coreStarvGenes_lt, function(x) GO_enrichment_fct(x))

# plot
GO_enrichment_tbl <- bind_rows(GO_enrichment_lt, .id = "comp_name")

GO_enrichment_tbl %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(W01Fisher <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>% 
  ggplot(aes(y=Term, x= signif_p, size=Cazymes_PCWDE_p, color=W01Fisher)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free")

ggsave(str_c(out_dir,"/GOenrichment_core_starvation_", SD, ".pdf"), device = "pdf", width = 8, height = 5)

coreStarvGO_tbl <-  GO_enrichment_tbl

## Swat and standard starvation ------------------------------------------------------------------------
swat_starvation <- read_tsv("All_comparisons_20230928/SWAT_DEGs_union.tsv")
swat_starvation <- swat_starvation$gene_id
length(swat_starvation)

stand_starvation <- read_tsv("All_comparisons_20230928/DM_LI_starvation.tsv")
stand_starvation <- stand_starvation$gene_id
length(stand_starvation)

swat_standard_StarvGenes_lt <- list(swat_starvation=swat_starvation, standard_starvation=stand_starvation)

GO_enrichment_lt <- lapply(swat_standard_StarvGenes_lt, function(x) GO_enrichment_fct(x))

# plot
GO_enrichment_tbl <- bind_rows(GO_enrichment_lt, .id = "comp_name")


GO_enrichment_tbl %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(W01Fisher <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>% 
  dplyr::rename("pvalue"="W01Fisher", "GeneRatio"="signif_p", "PCWDE_Ratio"="Cazymes_PCWDE_p") %>% 
  ggplot(aes(y=Term, x= GeneRatio, size=PCWDE_Ratio, color=pvalue)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free")

ggsave(str_c(out_dir,"/GOenrichment_swat_standard_starvation_", SD, ".pdf"), device = "pdf", width = 10, height = 10)

standardANDswatStarvGO_tbl <- GO_enrichment_tbl

# plot standard, swat, core together

bind_rows(standardANDswatStarvGO_tbl, coreStarvGO_tbl) %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(W01Fisher <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>% 
  dplyr::rename("pvalue"="W01Fisher", "GeneRatio"="signif_p", "PCWDE_Ratio"="Cazymes_PCWDE_p") %>% 
  mutate(comp_name = ifelse(comp_name=='coreStarvGenes','core_starvation',comp_name)) %>% 
  arrange(comp_name, GeneRatio) %>% 
  mutate(Term = factor(Term, levels = unique(Term), ordered = TRUE)) %>% 
  ggplot(aes(y=Term, x= GeneRatio, size=PCWDE_Ratio, color=pvalue)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free")

ggsave(str_c(out_dir,"/GOenrichment_core_standard_swat_starvation_Ordered_", SD, ".pdf"), device = "pdf", width = 10, height = 10)


# GO light induction --------------------------------------------------------------------------------

## All light induced samples -----------------------------------------------------------------------

all_light_lt <- c(all_AM_upreg_lt, all_AE_upreg_lt)

GO_enrichment_lt <- lapply(all_light_lt, function(x) GO_enrichment_fct(x))

GO_enrichment_tbl <- bind_rows(GO_enrichment_lt, .id = "comp_name")

GO_enrichment_tbl %>% 
  mutate(comp_name = str_remove(comp_name, pattern = "\\.signif\\.upreg$")) %>%
  mutate(comp_name = str_remove(comp_name, pattern = "_vs_.*")) %>%
  mutate(comp_name = factor(comp_name, levels = c("DM_156h_AE","ELI1hAM", "ELI2hAM", "LI2hAM", "LI6hAM", "LI12hAM", "LI18hAM", "LI_24h_AM",
                                                  "ELI1hAE", "ELI2hAE", "LI2hAE", "LI6hAE", "LI12hAE", "LI18hAE", "LI_24h_AE"), ordered = T)) %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(W01Fisher <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>% 
  ggplot(aes(y=Term, x= signif_p, size=Cazymes_PCWDE_p, color=W01Fisher)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

ggsave(str_c(out_dir,"/GOenrichment_LightInductionAll_PCWDWMark_", SD, ".pdf"), device = "pdf", width = 20, height = 15)

# plot 2 NOT mark PCWDE

GO_enrichment_tbl %>% 
  mutate(comp_name = str_remove(comp_name, pattern = "\\.signif\\.upreg$")) %>%
  mutate(comp_name = str_remove(comp_name, pattern = "_vs_.*")) %>%
  mutate(comp_name = factor(comp_name, levels = c("DM_156h_AE","ELI1hAM", "ELI2hAM", "LI2hAM", "LI6hAM", "LI12hAM", "LI18hAM", "LI_24h_AM",
                                                  "ELI1hAE", "ELI2hAE", "LI2hAE", "LI6hAE", "LI12hAE", "LI18hAE", "LI_24h_AE"), ordered = T)) %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(W01Fisher <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>% 
  ggplot(aes(y=Term, x=signif_p, size=Significant, color=W01Fisher)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  

ggsave(str_c(out_dir,"/GOenrichment_LightInductionAll_", SD, ".pdf"), device = "pdf", width = 20, height = 15)

## GO Core/Pure light induction (-Core starv) ----------------------------------------------------------------

pureLightInducedGenes <- read_tsv("All_comparisons_20230928/genes_specific_to_light_induction.tsv")
pureLightInducedGenes <- pureLightInducedGenes$gene_id
pureLightInducedGenes_lt <- list(pureLightInducedGenes=pureLightInducedGenes)

length(pureLightInducedGenes)

GO_enrichment_lt <- lapply(pureLightInducedGenes_lt, function(x) GO_enrichment_fct(x))

# plot
GO_enrichment_tbl <- bind_rows(GO_enrichment_lt, .id = "comp_name")

GO_enrichment_tbl %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(W01Fisher <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>% 
  ggplot(aes(y=Term, x= signif_p, size=Significant, color=W01Fisher)) +
  geom_point() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free") +
  theme_bw()

ggsave(str_c(out_dir,"/GOenrichment_core_light_induced_genes_cluster_", SD, ".pdf"), device = "pdf", width = 8, height = 5)

core_GO_enrichment <- GO_enrichment_tbl

## GO AM light induction ---------------------------------------------------------------
light_clust_tbl <- read_tsv("All_comparisons_20230928/Heatmap_AMAE_light_induced_genes_clusters_bordered_AllCoreLight_v2_20230928.tsv")

light_clust_tbl %>% 
  filter(lightInduced == 'AM') %>% 
  dplyr::select(gene_id, Cluster) -> am_light_clust


am_light_clust_lt <- split(am_light_clust$gene_id, am_light_clust$Cluster)

GO_enrichment_lt <- lapply(am_light_clust_lt, function(x) GO_enrichment_fct(x))

# plot
GO_enrichment_tbl <- bind_rows(GO_enrichment_lt, .id = "comp_name")

GO_enrichment_tbl %>% 
  dplyr::rename("pvalue"="W01Fisher") %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(pvalue <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>% 
  dplyr::rename("GeneRatio"="signif_p") %>% 
  arrange(desc(comp_name), GeneRatio) %>% 
  mutate(Term = factor(Term, levels = unique(Term), ordered = TRUE)) %>% 
  mutate(comp_name = factor(comp_name, levels=c('Cluster_2', 'Cluster_1'), ordered=TRUE)) %>% 
  ggplot(aes(y=Term, x= GeneRatio, size=Significant, color=pvalue)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free") +
  labs(size="Significant")

ggsave(str_c(out_dir,"/GOenrichment_attached_light_induced_genes_cluster_", SD, ".pdf"), device = "pdf", width = 9, height = 8)

LI_AM_GO_enrichment <- GO_enrichment_tbl

# plot 2
# AM + Core

core_GO_enrichment$comp_name <- 'Core'

bind_rows(LI_AM_GO_enrichment,core_GO_enrichment) %>% 
  dplyr::rename("pvalue"="W01Fisher") %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(pvalue <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>% 
  dplyr::rename("GeneRatio"="signif_p") %>% 
  arrange(comp_name, GeneRatio) %>% 
  mutate(Term = factor(Term, levels = unique(Term), ordered = TRUE)) %>% 
  ggplot(aes(y=Term, x= GeneRatio, size=Significant, color=pvalue)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free") +
  labs(size="Significant")

ggsave(str_c(out_dir,"/GOenrichment_attached_light_induced_genes_cluster_ANDcore_", SD, ".pdf"), device = "pdf", width = 9, height = 8)

## GO AE light induction ---------------------------------------------------------------

light_clust_tbl %>% 
  filter(lightInduced == 'AE') %>% 
  dplyr::select(gene_id, Cluster) -> ae_light_clust


ae_light_clust_lt <- split(ae_light_clust$gene_id, ae_light_clust$Cluster)

GO_enrichment_lt <- lapply(ae_light_clust_lt, function(x) GO_enrichment_fct(x))

# plot
GO_enrichment_tbl <- bind_rows(GO_enrichment_lt, .id = "comp_name")

GO_enrichment_tbl %>% 
  dplyr::rename("pvalue"="W01Fisher") %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(pvalue <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>% 
  dplyr::rename("GeneRatio"="signif_p") %>% 
  arrange(comp_name, GeneRatio) %>% 
  mutate(Term = factor(Term, levels = unique(Term), ordered = TRUE)) %>% 
  ggplot(aes(y=Term, x= GeneRatio, size=Significant, color=pvalue)) +
  geom_point() +
  theme_bw() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free") +
  labs(size="Significant")

ggsave(str_c(out_dir,"/GOenrichment_aerial_light_induced_genes_cluster_", SD, ".pdf"), device = "pdf", width = 9, height = 8)


## Check GO genes ----------------------------------------------------------------------------

getGenesByGO_fct <- function(geneSet, GOinterest) {
  geneList <- factor(as.integer(geneUniverse %in% geneSet))
  names(geneList) <- geneUniverse
  myGOdata <- new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize =10)
  result_tbl <- tibble(geneInGO = names(myGOdata@graph@nodeData@data[[GOinterest]]$genes)) 
  sg <- sigGenes(myGOdata)
  result_tbl %>% mutate(signifGenes = geneInGO%in% sg, GOinterest=GOinterest) -> result_tbl 
  
  return(result_tbl)
}

# GO:0006281 : DNA-repair

GenesByGO_lt <- lapply(am_light_clust_lt, FUN = getGenesByGO_fct, GOinterest='GO:0006281')
GenesByGO_tbl <- bind_rows(GenesByGO_lt, .id = 'groups')

write_tsv(GenesByGO_tbl, str_c(out_dir,"/GOlist_DNA_repair_", SD, ".tsv"))

# GO:0006629 : Lipid metabolic process

GenesByGO_lt <- lapply(am_light_clust_lt, FUN = getGenesByGO_fct, GOinterest='GO:0006629') # 177 lipid
GenesByGO_tbl <- bind_rows(GenesByGO_lt, .id = 'groups')

write_tsv(GenesByGO_tbl, str_c(out_dir,"/GOlist_LipidMetabolicProcess_", SD, ".tsv"))

count(GenesByGO_tbl, groups, signifGenes)
filter(LI_AM_GO_enrichment, GO.ID=='GO:0006629')

GenesByGOannot_tbl <- filter(GO_annot_long, Protein_accession %in% filter(GenesByGO_tbl, groups=='Cluster_2', signifGenes)$geneInGO)

# Mycelium specification --------------------------------------------------------------------------------

## upreg (AE upreg) -------------------------------------------------------------

GO_enrichment_lt <- lapply(all_micSpec_upreg_lt, function(x) GO_enrichment_fct(x))

# plot
GO_enrichment_tbl <- bind_rows(GO_enrichment_lt, .id = "comp_name")

GO_enrichment_tbl %>% 
  mutate(comp_name = str_remove(comp_name, pattern = "\\.signif\\.upreg$")) %>%
  mutate(comp_name = str_remove(comp_name, pattern = "_vs_.*")) %>%
  mutate(comp_name = factor(comp_name, levels = c("DM_156h_AE", "ELI1hAE", "ELI2hAE", "LI2hAE", "LI6hAE", "LI12hAE", "LI18hAE", "LI_24h_AE"), ordered = T)) %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(W01Fisher <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>%
  dplyr::rename("GeneRatio"="signif_p") %>% 
  arrange(comp_name, GeneRatio) %>% 
  mutate(Term = factor(Term, levels = unique(Term), ordered = TRUE)) %>% 
  ggplot(aes(y=Term, x= GeneRatio, size=Significant, color=W01Fisher)) +
  geom_point() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(str_c(out_dir,"/GOenrichment_Mycelium_specification_genes_AEupreg_Ordered_", SD, ".pdf"), device = "pdf", width = 13, height = 8)

## downreg (AM upreg) ----------------------------------------------------------------------------

GO_enrichment_lt <- lapply(all_micSpec_downreg_lt, function(x) GO_enrichment_fct(x))

# plot
GO_enrichment_tbl <- bind_rows(GO_enrichment_lt, .id = "comp_name")

GO_enrichment_tbl %>% 
  mutate(comp_name = str_remove(comp_name, pattern = "\\.signif\\.upreg$")) %>%
  mutate(comp_name = str_remove(comp_name, pattern = "_vs_.*")) %>%
  mutate(comp_name = str_replace(comp_name, pattern='AE$', replacement = 'AM')) %>% 
  mutate(comp_name = factor(comp_name, levels = c("DM_156h_AM", "ELI1hAM", "ELI2hAM", "LI2hAM", "LI6hAM", "LI12hAM", "LI18hAM", "LI_24h_AM"), ordered = T)) %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(W01Fisher <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>%
  
  dplyr::rename("GeneRatio"="signif_p") %>% 
  arrange(comp_name, GeneRatio) %>% 
  mutate(Term = factor(Term, levels = unique(Term), ordered = TRUE)) %>% 
  
  ggplot(aes(y=Term, x= GeneRatio, size=Significant, color=W01Fisher)) +
  geom_point() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(str_c(out_dir,"/GOenrichment_Mycelium_specification_genes_AMupreg_Ordered_", SD, ".pdf"), device = "pdf", width = 12, height = 18)


## Core AE/AM genes ----------------------------------------------------------------------------

coreAE_tbl <- read_tsv('All_comparisons_20230928/Attached_mycelium_specificCoreGenesWithIproAnnot_upreg_20230928.tsv') 
coreAM_tbl <- read_tsv('All_comparisons_20230928/Attached_mycelium_specificCoreGenesWithIproAnnot_downreg_20230928.tsv') 

coreAEAM_list <- list(coreAE = coreAE_tbl$gene_id,
                      coreAM = coreAM_tbl$gene_id)

GO_enrichment_lt <- lapply(coreAEAM_list, function(x) GO_enrichment_fct(x))

# plot
GO_enrichment_tbl <- bind_rows(GO_enrichment_lt, .id = "comp_name")

GO_enrichment_tbl %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(W01Fisher <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>%
  ggplot(aes(y=Term, x= signif_p, size=Significant, color=W01Fisher)) +
  geom_point() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(str_c(out_dir,"/GOenrichment_Mycelium_specification_genes_AEAMcore_", SD, ".pdf"), device = "pdf", width = 10, height = 6)


## Core AE/AM NO starv Core genes ----------------------------------------------------------------------------

coreAEAMnoCoreStarv_list <- lapply(coreAEAM_list, function(x) {
  return(x[! x %in% coreStarvGenes])
})

names(coreAEAMnoCoreStarv_list) <- str_c(names(coreAEAMnoCoreStarv_list), '_noStarv')

sapply(coreAEAMnoCoreStarv_list, length)

GO_enrichment_lt <- lapply(coreAEAMnoCoreStarv_list, function(x) GO_enrichment_fct(x))

# plot
GO_enrichment_tbl <- bind_rows(GO_enrichment_lt, .id = "comp_name")

GO_enrichment_tbl %>% 
  mutate(comp_name = str_remove(comp_name, pattern='_noStarv$')) %>% 
  mutate(signif_p = Significant/Annotated) %>% 
  filter(W01Fisher <= 0.05) %>% 
  mutate(Term = str_c(Term, ' ', GO.ID)) %>%
  dplyr::rename("GeneRatio"="signif_p") %>% 
  arrange(comp_name, GeneRatio) %>% 
  mutate(Term = factor(Term, levels = unique(Term), ordered = TRUE)) %>% 
  ggplot(aes(y=Term, x= GeneRatio, size=Significant, color=W01Fisher)) +
  geom_point() +
  scale_color_gradient(low="red", high="blue") +
  facet_grid(Ontology ~ comp_name, scales = "free", space = "free") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(str_c(out_dir,"/GOenrichment_Mycelium_specification_genes_AEAMcore_NOstarv_Ordered_", SD, ".pdf"), device = "pdf", width = 8, height = 7)


filter(coreAE_tbl, gene_id %in% coreAEAMnoCoreStarv_list$coreAE_noStarv) %>% View()
filter(coreAM_tbl, gene_id %in% coreAEAMnoCoreStarv_list$coreAM_noStarv) %>% View()

coreAEAMnoCoreStarv_Genelist <- lapply(coreAEAMnoCoreStarv_list, function(x) {
  result_tbl <- tibble(group_n=length(x),
                       gene_id=str_remove_all(str_c(x, collapse = ','), pattern = '\\.T0'))
  return(result_tbl)
})
coreAEAMnoCoreStarv_Genelist <- bind_rows(coreAEAMnoCoreStarv_Genelist, .id = 'comp')

write_tsv(coreAEAMnoCoreStarv_Genelist, str_c(out_dir,"/Mycelium_specification_genesList_AEAMcore_NOstarv_", SD, ".tsv"))

