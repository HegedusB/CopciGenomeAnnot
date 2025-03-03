# load packages -------------------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicFeatures)
  library(rtracklayer)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(knitr)
  library(tidyr)
})



# load data --------------------------------------------------------------------------------------------------
# load annots
longread_annot_path <- './CopciAB_new_annot_CDS_20220425m2_newIds_isoforms_mancur_filtered.gtf'
xie_annot_path <- "../ref_genome/fullgenome/Coprinopsis_cinerea_A43mutB43mut_pab1-1_326_mHB_2.gtf"
jgi_annot_path <- "../ref_genome/Copci_AmutBmut1_GeneModels_FrozenGeneCatalog_20160912_mHB.gtf"
okayama_annot_path <- '../ref_genome/Copci_Okayama/Copci1.filtered_proteins.BroadModels1.gff'

longread_annot_gr <- import(longread_annot_path)
xie_annot_gr <- import(xie_annot_path)
jgi_annot_gr <- import(jgi_annot_path)
okayama_annot_gr <- import(okayama_annot_path)

jgi_annot_gr$transcript_id <- str_c(jgi_annot_gr$transcript_id, '.T0')
jgi_annot_gr$gene_id <- str_c(jgi_annot_gr$gene_id, '.T0')
length(unique(str_remove(longread_annot_gr$gene_id, pattern = '\\.T[0-9]+$')))
length(unique(longread_annot_gr$transcript_id))

longread_annot_T0_gr <- longread_annot_gr[str_detect(longread_annot_gr$transcript_id, pattern = "T0$")]
xie_annot_T1_gr <- xie_annot_gr[str_detect(xie_annot_gr$transcript_id, pattern = "T1$")]

longread_annot_T0_df <- as.data.frame(longread_annot_T0_gr)
xie_annot_T1_df <- as.data.frame(xie_annot_T1_gr)
jgi_annot_df <- as.data.frame(jgi_annot_gr)
okayama_annot_df <- as.data.frame(okayama_annot_gr)

xie_annot_T1_df %>% filter(type=="CDS") %>% .$gene_id %>% unique() %>% length()

# annot summary ---------------------------------------------------------------------------------------

annot_lt <- list(longread = longread_annot_T0_gr,
                 xie = xie_annot_T1_gr,
                 jgi =  jgi_annot_gr)

annot_okayama_lt <- list(okayama = okayama_annot_gr)

make_annotation_summary_fct <- function(annot_gr) {
  annot_TxDb <- makeTxDbFromGRanges(annot_gr)
  
  # get data
  annot_genes_gr <- genes(annot_TxDb)
  annot_exon_gr <- exonsBy(annot_TxDb, by = "tx", use.names=TRUE)
  annot_cds_gr <- cdsBy(annot_TxDb, by = "tx", use.names=TRUE)
  annot_introns_gr <- intronsByTranscript(annot_TxDb, use.names=TRUE)
  annot_5UTR_gr <- fiveUTRsByTranscript(annot_TxDb, use.names=TRUE)
  annot_3UTR_gr <- threeUTRsByTranscript(annot_TxDb, use.names=TRUE)
  
  # summarise data
  ## gene
  gene_tbl <- tibble(number = length(annot_genes_gr),
                     average_number = NA,
                     length = sum(width(annot_genes_gr)),
                     average_length = mean(width(annot_genes_gr)))
  
  # exon
  exon_tbl <- tibble(number = sum(sapply(width(annot_exon_gr), length)),
                     average_number = number/gene_tbl$number,
                     length = sum(sum(width(annot_exon_gr))),
                     average_length = length/gene_tbl$number)
  
  # cds
  cds_tbl <- tibble(number = sum(sapply(width(annot_cds_gr), length)),
                    average_number = number/gene_tbl$number,
                    length = sum(sum(width(annot_cds_gr))),
                    average_length = length/gene_tbl$number)
  
  introns_tbl <- tibble(number = sum(sapply(width(annot_introns_gr), length)),
                        average_number = number/gene_tbl$number,
                        length = sum(sum(width(annot_introns_gr))),
                        average_length = length/gene_tbl$number)
  
  
  utr5_tbl <- tibble(number = length(annot_5UTR_gr), 
                     average_number = number/gene_tbl$number, 
                     length = sum(sum(width(annot_5UTR_gr))),
                     average_length = length/number) 
                    

  
  utr3_tbl <- tibble(number = length(annot_3UTR_gr), 
                     average_number = number/gene_tbl$number,
                     length = sum(sum(width(annot_3UTR_gr))),
                     average_length = length/number)

  result_lt <- list(gene = gene_tbl,
                    exon = exon_tbl,
                    cds = cds_tbl,
                    intron = introns_tbl,
                    utr3 = utr3_tbl,
                    utr5 = utr5_tbl)
  
  result_tbl <- bind_rows(result_lt, .id="source")
  return(result_tbl)
}


make_annotationOkayama_summary_fct <- function(annot_gr) {
  annot_TxDb <- makeTxDbFromGRanges(annot_gr)
  
  # get data
  annot_genes_gr <- genes(annot_TxDb)
  annot_exon_gr <- exonsBy(annot_TxDb, by = "tx", use.names=TRUE)
  annot_cds_gr <- cdsBy(annot_TxDb, by = "tx", use.names=TRUE) # 14750
  annot_introns_gr <- intronsByTranscript(annot_TxDb, use.names=TRUE)
  annot_5UTR_gr <- fiveUTRsByTranscript(annot_TxDb, use.names=TRUE)
  annot_3UTR_gr <- threeUTRsByTranscript(annot_TxDb, use.names=TRUE)
  
  # summarise data
  ## gene
  gene_tbl <- tibble(number = length(annot_genes_gr),
                     average_number = NA,
                     length = sum(width(annot_genes_gr)),
                     average_length = mean(width(annot_genes_gr)))
  
  # exon
  exon_tbl <- tibble(number = sum(sapply(width(annot_exon_gr), length)),
                     average_number = number/gene_tbl$number,
                     length = sum(sum(width(annot_exon_gr))),
                     average_length = length/gene_tbl$number)
  
  # cds
  cds_tbl <- tibble(number = sum(sapply(width(annot_cds_gr), length)),
                    average_number = number/gene_tbl$number,
                    length = sum(sum(width(annot_cds_gr))+3), 
                    average_length = length/gene_tbl$number)
  
  introns_tbl <- tibble(number = sum(sapply(width(annot_introns_gr), length)),
                        average_number = number/gene_tbl$number,
                        length = sum(sum(width(annot_introns_gr))),
                        average_length = length/gene_tbl$number)
  

  utr5_tbl <- tibble(number = length(annot_5UTR_gr), 
                     average_number = number/gene_tbl$number, 
                     length = sum(sum(width(annot_5UTR_gr))), 
                     average_length = length/number) 
  

  annot_3UTR_gr <- annot_3UTR_gr[sum(width(annot_3UTR_gr)) != 3] 
  utr3_tbl <- tibble(number = length(annot_3UTR_gr),
                     average_number = number/gene_tbl$number,
                     length = sum(sum(width(annot_3UTR_gr))-3), 
                     average_length = length/number)

  result_lt <- list(gene = gene_tbl,
                    exon = exon_tbl,
                    cds = cds_tbl,
                    intron = introns_tbl,
                    utr3 = utr3_tbl,
                    utr5 = utr5_tbl)
  
  result_tbl <- bind_rows(result_lt, .id="source")
  return(result_tbl)
}

annot_summary_lt <- mclapply(annot_lt, make_annotation_summary_fct, mc.cores = 1)
annot_okayama_summary_lt <- mclapply(annot_okayama_lt, make_annotationOkayama_summary_fct, mc.cores = 1)

annot_summary_lt <- c(annot_summary_lt, annot_okayama_summary_lt)

annot_summary_tbl <- bind_rows(annot_summary_lt, .id="annot_variant")

annot_summary_tbl %>% arrange(source) -> annot_summary_tbl


annot_summary_copy_tbl <- annot_summary_tbl
annot_summary_copy_tbl$annot_variant <- c('JGI_V2', 'JGI_V1', 'Xie', 'Okayama')[match(annot_summary_copy_tbl$annot_variant, c('longread', 'jgi', 'xie', 'okayama'))]

write_tsv(annot_summary_copy_tbl, 'compare_refannots_properties_20230619.tsv')

# microexon summary ------------------------------------


## short exon (<=51 nt) -------------------------------------------

make_short_exon_summary3_fct <- function(annot_gr) {
  annot_TxDb <- txdbmaker::makeTxDbFromGRanges(annot_gr)
  
  # get data
  annot_exon_gr <- exonsBy(annot_TxDb, by = "tx", use.names=TRUE)

  # summarise data
  exon_width_lt <- width(annot_exon_gr)
  genes_with_short_exons_lt <- lapply(exon_width_lt, function(y) {
    exon_length_index <- y<=51 
    exon_summary = list(short_terminal_exon_n = integer(0),
                        short_internal_exon_n = integer(0))
    if(length(unlist(y)) <= 2) { 
      if(any(exon_length_index)) { 
        exon_summary$short_terminal_exon_n <- unlist(y[exon_length_index], use.names = FALSE)
      }
    } else { 
      
      exon_summary = list(short_terminal_exon_n = unlist(y)[c(1,length(unlist(y)))][unlist(exon_length_index)[c(1,length(unlist(y)))]],
                          short_internal_exon_n = unlist(y)[-c(1,length(unlist(y)))][unlist(exon_length_index)[-c(1,length(unlist(y)))]])
    }
    return(exon_summary)
  })
  
  length(unlist(lapply(genes_with_short_exons_lt, "[[", "short_internal_exon_n")))
  length(unlist(lapply(genes_with_short_exons_lt, "[[", "short_terminal_exon_n"))) 
  
  result_lt <- list(short_internal_exon_n = unlist(lapply(genes_with_short_exons_lt, "[[", "short_internal_exon_n")),
                    short_terminal_exon_n = unlist(lapply(genes_with_short_exons_lt, "[[", "short_terminal_exon_n")))
  
  return(result_lt)
}

short_exon_summary3_lt <- mclapply(annot_lt, make_short_exon_summary3_fct, mc.cores = 1)

lapply(short_exon_summary3_lt, "[[", "short_internal_exon_n") %>% 
  lapply(as.data.frame) %>% 
  lapply(setNames, nm = "exon_length") %>% 
  bind_rows(.id = "annotation") %>% 
  remove_rownames() -> short_exon_summary3_internal_tbl

dplyr::count(short_exon_summary3_internal_tbl, annotation)

short_exon_summary3_internal_tbl %>% 
  ggplot(aes(exon_length, fill=annotation)) +
  geom_bar() +
  facet_wrap(.~annotation) +
  theme_bw()

short_exon_summary3_internal_tbl %>% 
  mutate(annotation = as.character(factor(annotation, labels = c('JGI_V1', 'JGI_V2', 'Xie')))) %>% 
  ggplot(aes(exon_length, color=annotation)) +
  geom_density() +
  xlab('internal_exon_length') +
  theme_bw()

ggsave('microExonUpTo51ntDensityPlot.pdf', device = 'pdf', width = 8, height = 5)

short_exon_summary3_internal_tbl %>% 
  mutate(annotation = as.character(factor(annotation, labels = c('JGI_V1', 'JGI_V2', 'Xie')))) %>% 
  dplyr::count(annotation, exon_length, name = 'exon_count') %>%
  ggplot(aes(x = exon_length, y = exon_count, color=annotation, group=annotation, fill=annotation)) +
  geom_point() +
  geom_line(linetype=1) +
  geom_smooth(se = FALSE,  method = 'loess', formula = 'y ~ x') +
  xlab('internal_exon_length') +
  ylab('internal_exon_count') +
  theme_bw()


ggsave('microExonUpTo51ntCountPointPlot.pdf', device = 'pdf', width = 9, height = 5)


short_exon_summary3_internal_tbl %>% 
  mutate(exon_length_triplet = exon_length %% 3) %>% 
  dplyr::count(annotation, exon_length_triplet)



short_exon_summary3_internal_tbl %>% 
  dplyr::count(annotation, exon_length) %>% 
  spread(key = exon_length, value = n)


lapply(short_exon_summary3_lt, "[[", "short_internal_exon_n") %>% 
  lapply(as.data.frame) %>% 
  lapply(setNames, nm = "exon_length") %>% 
  bind_rows(.id = "annotation") -> short_exon_summary3_internal_withRownames_tbl


short_exon_summary3_internal_withRownames_tbl %>% 
  rownames_to_column('geneID') %>% 
  mutate(geneID = str_remove(geneID, pattern = '\\.T[0-9]+$')) %>% 
  mutate(geneID = str_remove(geneID, pattern = '-T[0-9]+$')) -> short_exon_summary3_internal_withRownames_tbl

short_exon_summary3_internal_withRownames_tbl %>% 
  dplyr::count(annotation, geneID) %>% # 
  dplyr::count(annotation, n) %>% #
  spread(key = n, value = nn)

## microexons (<=15 nt)  -------------------------------------------------------------------

short_exon_summary3_internal_withRownames_tbl %>% 
  filter(exon_length <= 15) %>% 
  group_by(annotation) %>% 
  summarise(geneID_n = length(unique(geneID)))

short_exon_summary3_internal_withRownames_tbl %>% 
  group_by(annotation, exon_length) %>% 
  summarise(geneID_n = length(unique(geneID))) %>% 
  spread(key = exon_length, value = geneID_n)

write_tsv(short_exon_summary3_internal_withRownames_tbl, 'internalMicroexons_colection.tsv')


library(GenomicFeatures)
library(BSgenome)
library(ORFik)
library(parallel)

genome_path <- '../ref_genome/fullgenome_jgi/merged_genome/CopciAB_new_jgi_20220113.fasta' 
genome_seq <- readDNAStringSet(genome_path)

genomeO_path <- "../ref_genome/Copci_Okayama/Coprinopsis_cinerea.masked.fasta"
genomeO_seq <- readDNAStringSet(genomeO_path)

short_exon_summary3_internal_withRownames_15_tbl <- filter(short_exon_summary3_internal_withRownames_tbl, exon_length <= 15)
short_exon_summary3_internal_withRownames_15Longread_tbl <- filter(short_exon_summary3_internal_withRownames_15_tbl,annotation == 'longread')
short_exon_summary3_internal_withRownames_15Longread_lt <- split(short_exon_summary3_internal_withRownames_15Longread_tbl, short_exon_summary3_internal_withRownames_15Longread_tbl$geneID)

table(sapply(short_exon_summary3_internal_withRownames_15Longread_lt, function(x) any(duplicated(x$exon_length)))) # van olyan gén amely hasonló hosszúságú kis exonokat tartalmaz
short_exon_summary3_internal_withRownames_15Longread_lt[sapply(short_exon_summary3_internal_withRownames_15Longread_lt, function(x) any(duplicated(x$exon_length)))]


smallESum_lt <- lapply(short_exon_summary3_internal_withRownames_15Longread_lt, function(x) {

  tempGeneID <- unique(x$geneID)
  tempAnnot_gr <- annot_lt$longread[str_remove(annot_lt$longread$gene_id, pattern = '\\.T0$') %in% unique(tempGeneID)]
  tempAnnotShort_gr <- tempAnnot_gr[tempAnnot_gr$type %in% c('exon', 'CDS')]
  
  internalExon <- tempAnnotShort_gr[tempAnnotShort_gr$type == 'exon',]$exon_id
  internalExon_index <- as.integer(sapply(str_split(internalExon, pattern = "\\."), "[[", 3))
  internalExon <- internalExon[!internalExon_index %in% c(1, max(internalExon_index))]
  exonIndex_gr <- tempAnnotShort_gr[tempAnnotShort_gr$exon_id %in% internalExon] 
  exonIndex_gr <- exonIndex_gr[width(exonIndex_gr) <= 15] 
  exonIndex_gr <- exonIndex_gr[exonIndex_gr$type == "exon"] 
  
  exonIndex <- exonIndex_gr$exon_id
  
  result_inner_lt <- lapply(exonIndex, function(x) {
    exonIndex_inner_gr <- tempAnnotShort_gr[tempAnnotShort_gr$exon_id == x]
    if(any(exonIndex_inner_gr$type == "CDS")) {
      
      tempAnnot_fexon_gr <- tempAnnot_gr[!tempAnnot_gr$exon_id %in% x] 
      tempAnnot_txdb <- txdbmaker::makeTxDbFromGRanges(tempAnnot_fexon_gr)
      cbt <- cdsBy(tempAnnot_txdb, by='tx', use.name=TRUE) 
      cds_seq <- extractTranscriptSeqs(genome_seq, cbt)
      
      man_ORFs <- findORFs(cds_seq, 
                           startCodon = "ATG",
                           stopCodon = stopDefinition(1),
                           longestORF = TRUE, minimumLength = 30) 
      
      if(length(man_ORFs) > 0) {

        result_inner = tibble(exonID = x,
                              exon_length = width(tempAnnot_gr[tempAnnot_gr$exon_id %in% x & tempAnnot_gr$type == 'exon']),
                              cds_seq_width = width(cds_seq),
                              max_ORF_width = max(width(man_ORFs[[1]]))[[1]])
        return(result_inner)
      } else {
        result_inner = tibble(exonID = x,
                              exon_length = width(tempAnnot_gr[tempAnnot_gr$exon_id %in% x & tempAnnot_gr$type == 'exon']),
                              cds_seq_width = width(cds_seq),
                              max_ORF_width = 0)
        return(result_inner)
      }
    } else {
      result_inner = tibble(exonID = x,
                            exon_length = width(tempAnnot_gr[tempAnnot_gr$exon_id %in% x & tempAnnot_gr$type == 'exon']),
                            cds_seq_width = 0, 
                            max_ORF_width = 0)
      return(result_inner) 
    }
  })
  
  result_inner_tbl <- bind_rows(result_inner_lt)
  return(result_inner_tbl)
})


smallESum_tbl <- bind_rows(smallESum_lt)
smallESum_tbl$geneID <- str_remove(smallESum_tbl$exonID, pattern = "\\.T.*")


save(smallESum_lt, file = 'exon_cds_orf_test_data.RData')
write_tsv(smallESum_tbl, 'exon_cds_orf_test_data.tsv')


filter(smallESum_tbl, cds_seq_width == 0, max_ORF_width == 0) 
filter(smallESum_tbl, cds_seq_width != 0, max_ORF_width == 0) 

filter(smallESum_tbl, cds_seq_width != 0, cds_seq_width == max_ORF_width)
filter(smallESum_tbl, cds_seq_width != 0, cds_seq_width > max_ORF_width) 
filter(smallESum_tbl, cds_seq_width != 0, cds_seq_width < max_ORF_width) 

smallESum_tbl %>% 
  mutate(description = NA) %>% 
  mutate(description = ifelse(cds_seq_width == 0, 'noCDSonSmallExon', description)) %>% 
  mutate(description = ifelse(cds_seq_width != 0 & cds_seq_width == max_ORF_width, 'noalternativeStop', description)) %>% 
  mutate(description = ifelse(cds_seq_width != 0 & cds_seq_width > max_ORF_width, 'alternativeStop', description)) -> smallESumDesc_tbl

smallESumDesc_tbl %>% 
  mutate(divisible_3 = (exon_length %% 3 == 0)) -> smallESumDesc_tbl

dplyr::count(smallESumDesc_tbl, description, divisible_3)

write_tsv(smallESumDesc_tbl, 'exon_cds_orf_test_WithDescription_data.tsv')

smallESumDesc_tbl %>% 
  ggplot(aes(exon_length, fill=description)) +
  geom_bar(stat='count') +
  facet_grid(~description+divisible_3)


length(unique(smallESumDesc_tbl$geneID)) 
length(unique(filter(smallESumDesc_tbl, exon_length<=3)$geneID)) 

length(unique(filter(smallESumDesc_tbl, description == 'alternativeStop')$geneID))
length(unique(filter(smallESumDesc_tbl, description == 'alternativeStop')$geneID)) 
length(unique(filter(smallESumDesc_tbl, description == 'noalternativeStop')$geneID))

length(unique(filter(smallESumDesc_tbl, description == 'alternativeStop', exon_length <= 3)$geneID)) 


filter(smallESumDesc_tbl, exon_length==1)$geneID %>% length()

filter(smallESumDesc_tbl, exon_length==1)$geneID %>% unique() %>% length()


filter(smallESumDesc_tbl, exon_length==2)$geneID %>% length()

filter(smallESumDesc_tbl, exon_length==2)$geneID %>% unique() %>% length()


filter(smallESumDesc_tbl, exon_length==3)$geneID %>% length()

filter(smallESumDesc_tbl, exon_length==3)$geneID %>% unique() %>% length()



## microexon position ---------------------------------------------------------------------------------------------


smallEOSum_lt <- lapply(short_exon_summary3_internal_withRownames_15Longread_lt, function(x) {

  tempGeneID <- unique(x$geneID)
  tempAnnot_gr <- annot_lt$longread[str_remove(annot_lt$longread$gene_id, pattern = '\\.T0$') %in% tempGeneID]
  tempAnnotShort_gr <- tempAnnot_gr[tempAnnot_gr$type %in% c('exon', 'CDS')]
  

  internalIndex <- c(1, max(as.integer(tempAnnotShort_gr$exon_number))) # ha nem a transzkript elso vagy utolsó exonja 
  internal_gr <- tempAnnotShort_gr[!tempAnnotShort_gr$exon_number %in% as.character(internalIndex),]
  exonIndex <- internal_gr[internal_gr$type == 'exon',]$exon_id[width(internal_gr[internal_gr$type == 'exon']) %in% x$exon_length] # ne vegyem figyelembe a CDS hosszat is csak az exon hosszat
  

  GOIannot_txdb <- txdbmaker::makeTxDbFromGRanges(tempAnnotShort_gr)
  ebt <- exonsBy(GOIannot_txdb, by='tx', use.name=TRUE)
  ebt_df <- as.data.frame(ebt)
  ebt_df <- arrange(ebt_df, exon_rank)
  
  tscript_length <- sum(ebt_df$width)
  
  
  ebt_df$exonStart <-  c(1, (cumsum(ebt_df$width)[-length(ebt_df$width)])+1)
  ebt_df$exonEnd <- cumsum(ebt_df$width)
  

  if(!all((ebt_df$exonEnd - ebt_df$exonStart)+1 == ebt_df$width)) stop("Exon length problems!")
  
  microExon_tbl <- filter(ebt_df, exon_name %in% exonIndex)
  microExon_lt <- split(microExon_tbl, microExon_tbl$exon_name)
  
  

  tacriptOlap_lt <- lapply(microExon_lt, function(x) {
  
    
    result <- cut(seq(x$exonStart,x$exonEnd,1), round(quantile(seq(1,tscript_length,1), seq(0,1,0.1))), 
                  right = TRUE, include.lowest = TRUE, labels = str_c(seq(10,100,10), "%"))
    result <- str_c(unique(as.character(result)), collapse = "_")
    return(result)
  })
  
  tacriptOlap_tbl <- tibble(exon_name = names(tacriptOlap_lt),
                            tacriptOlap = unlist(tacriptOlap_lt))
  
  
  microExonResult_tbl <- left_join(microExon_tbl, tacriptOlap_tbl, by="exon_name")
  return(microExonResult_tbl)
})


smallEOSum_tbl <- bind_rows(smallEOSum_lt)
dplyr::count(smallEOSum_tbl, tacriptOlap)

write_tsv(smallEOSum_tbl, 'internalMicroexon_transcript_exoneverlap_percent.tsv')

tacriptOlapSplit_lt <- str_split(smallEOSum_tbl$tacriptOlap, pattern = "_")

tacriptOlapSplitRate_lt <- lapply(tacriptOlapSplit_lt, function(x) {
  result_tbl <- tibble(olapPcent = x)
  result_tbl <- mutate(result_tbl, olapPcentRate = 1/n())
  return(result_tbl)
})


tacriptOlapSplitRate_tbl <- bind_rows(tacriptOlapSplitRate_lt)
tacriptOlapSplitRate_tbl %>% 
  group_by(olapPcent) %>% 
  summarise(olapPcentRateSum = sum(olapPcentRate)) %>% 
  ungroup() %>% 
  arrange(nchar(olapPcent), olapPcent) %>% 
  mutate(olapPcent = factor(olapPcent, levels = olapPcent, ordered = TRUE)) %>% 
  ggplot(aes(x=olapPcent, y=olapPcentRateSum)) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("10% transcript segments") +
  ylab("Sum of overlapping microexon proportions") +
  theme_bw()

ggsave("overlappingMicroexPropsInTSegments.pdf", device = "pdf", height = 5, width = 5)


# microexon interpro enrichment -----------------------------------

ipro_annot_path <- "C:\\Users\\hbotond\\OneDrive - brc.hu\\Cikkek_Saját\\pj8_CopciAB\\0_Writing_Manuscript\\forCellGenomics\\CellGenomicsJuly1\\extraData\\CopciAB_new_annot_CDS_20220425m1_iproscan88_headered.txt"

ipro_annot <- read_tsv(ipro_annot_path)
ipro_annot_T0 <- filter(ipro_annot, str_detect(Protein_accession, pattern = 'T0$'))
ipro_annot_T0_ipr <- filter(ipro_annot_T0, !str_detect(InterPro_annotations_accession, pattern = '-'))
ipro_annot_T0_ipr <- mutate(ipro_annot_T0_ipr, Protein_accession = str_remove(Protein_accession, pattern='\\.T0$'))

ipro_annot_T0_ipr %>% 
  group_by(Protein_accession) %>% 
  filter(!duplicated(InterPro_annotations_accession)) %>% 
  ungroup() %>% 
  select(Protein_accession, InterPro_annotations_accession, InterPro_annotations_description) -> ipro_annot_T0_ipr_clean

short_internalExonOnlyLongR  <- short_exon_summary3_internal_withRownames_15Longread_tbl
short_internalExonOnlyLongR <- filter(short_internalExonOnlyLongR, exon_length<=15)

ipro_annot_T0_ipr_cleanP <- mutate(ipro_annot_T0_ipr_clean, shortExon_index = Protein_accession %in% short_internalExonOnlyLongR$geneID)

iprodict <- filter(ipro_annot_T0_ipr_cleanP, !duplicated(InterPro_annotations_accession)) %>% select(InterPro_annotations_accession, InterPro_annotations_description)

temp_lt <- mclapply(unique(ipro_annot_T0_ipr_cleanP$InterPro_annotations_accession), function(x) {

  NKIP <- filter(ipro_annot_T0_ipr_cleanP, shortExon_index==FALSE, InterPro_annotations_accession ==x)

  NKNIP <- filter(ipro_annot_T0_ipr_cleanP, shortExon_index==FALSE, InterPro_annotations_accession !=x)

  KIP <- filter(ipro_annot_T0_ipr_cleanP, shortExon_index==TRUE, InterPro_annotations_accession ==x)

  KNIP <- filter(ipro_annot_T0_ipr_cleanP, shortExon_index==TRUE, InterPro_annotations_accession !=x)
  
  cmatrix = matrix(c(nrow(KIP), nrow(KNIP), nrow(NKIP), nrow(NKNIP)), nrow = 2, byrow = TRUE)
  pvalue <- fisher.test(cmatrix, alternative = 'greater')$p.value
  
  result_tbl <- tibble(iproID = x, 
                       nonSmallExon = nrow(NKIP), 
                       SmallExon = nrow(KIP),
                       pvalue = pvalue)
  return(result_tbl)
}, mc.cores = 1)

bind_rows(temp_lt) %>% 
  arrange(pvalue) %>% 
  mutate(padjust = p.adjust(pvalue, method='fdr')) %>% 
  left_join(iprodict, by=c('iproID'='InterPro_annotations_accession')) -> temp_summary

filter(temp_summary, padjust<=0.05) %>% 
  arrange(desc(padjust)) %>% 
  mutate(iproID_plus = str_c(iproID, InterPro_annotations_description, sep = '__'),
         padjustm10 = -log10(padjust),
         ratio = SmallExon/(SmallExon+nonSmallExon)) -> temp_summary_signif

View(temp_summary_signif)

p1 <- temp_summary_signif %>% 
  mutate(iproID_plus = factor(iproID_plus, levels=iproID_plus, ordered = T)) %>% 
  ggplot(aes(y=iproID_plus, x=padjustm10)) +
  geom_point(aes(color=ratio, size=SmallExon)) +
  theme_bw()

ggsave(plot = p1, filename = 'smallexon15_ipri_Fishertest.pdf', device = 'pdf', width = 12, height = 6)


# Cytochrome P450 interpro

ipro_annot_T0_ipr_cleanP %>% 
  filter(InterPro_annotations_accession == "IPR017972") -> IPR017972_tbl

ipro_annot_T0_ipr_cleanP %>% 
  filter(InterPro_annotations_accession == "IPR002401") -> IPR002401_tbl

ipro_annot_T0_ipr_cleanP %>% 
  filter(InterPro_annotations_accession == "IPR001128") -> IPR001128_tbl

ipro_annot_T0_ipr_cleanP %>% 
  filter(InterPro_annotations_accession == "IPR036396") -> IPR036396_tbl


P450_tbl <- bind_rows(list(IPR017972_tbl, IPR002401_tbl, IPR001128_tbl, IPR036396_tbl))
write_csv(P450_tbl, "signifMicroexonContainingCytochromeP450.tsv")


P450_tbl %>% dplyr::count(InterPro_annotations_accession, shortExon_index)


# CopciAB V2 vs. CopciO ---------------------------------------------------------


okayama_annotgff3_path <- '../ref_genome/Copci_Okayama/Copci1.filtered_proteins.BroadModels1.gff3'
okayama_annotgff3_gr <- import(okayama_annotgff3_path)
okayama_annotgff3_df <- as.data.frame(okayama_annotgff3_gr)

okayama_annotgff3_df$ParentNew <- as.character(str_replace(okayama_annotgff3_df$Parent, pattern="gene", replacement="mRNA"))
okayama_annotgff3_df <- okayama_annotgff3_df[okayama_annotgff3_df$ParentNew != "character(0)",]
okayama_annotgff3_df %>% 
  group_by(ParentNew) %>% 
  mutate(proteinId = proteinId[!is.na(proteinId)]) %>% 
  filter(type == "exon") -> okayama_annotgff3Exon_df


okayama_aa_seq <- readAAStringSet("Coprinopsis_cinerea.proteins.fasta")
okayamaLookup_tbl <- tibble(name = sapply(str_split(names(okayama_aa_seq), pattern = "\\|"), "[[", 4),
                            protID = sapply(str_split(names(okayama_aa_seq), pattern = "\\|"), "[[", 3))

okayamaLookup_tbl %>% 
  filter(str_detect(name, pattern = "T0$")) %>% 
  mutate(name=str_remove(name, pattern = "T0$")) -> okayamaLookup_tbl

CopciAB_dict <- read_tsv("Copci_AB_id_dictionary_t0.tsv")
unique(smallESumDesc_tbl$geneID)
CopciABSmallExon_dict <- filter(CopciAB_dict, gene_id %in% unique(smallESumDesc_tbl$geneID))

CopciABSmallExonNoNA_dict <- filter(CopciABSmallExon_dict, !is.na(CopciABSmallExon_dict$target)) %>% 
  select("gene_id", "transcript_id", "target")

CopciABSmallExonNoNA_dict <- left_join(CopciABSmallExonNoNA_dict, okayamaLookup_tbl, by=c("target"="name"))

okayama_annotgff3ExonMicroInAB_df <- okayama_annotgff3Exon_df[okayama_annotgff3Exon_df$proteinId %in% CopciABSmallExonNoNA_dict$protID,]
okayama_annotgff3ExonMicroInAB_lt <- split(okayama_annotgff3ExonMicroInAB_df, okayama_annotgff3ExonMicroInAB_df$proteinId)


# quick look on the microexon contetnt of the Okayama annotation (terminal+internal)
table(sapply(okayama_annotgff3ExonMicroInAB_lt, function(x) any(x$width <= 15)))
okayama_annotgff3ExonMicroInAB_lt[!sapply(okayama_annotgff3ExonMicroInAB_lt, function(x) any(x$width <= 15))]


###########################
library("BSgenome")
library("pwalign")


smallESumDesc_tbl 
nrow(smallESumDesc_tbl)
length(unique(smallESumDesc_tbl$geneID))

compMicroEx_fct <- function(x) {

  Okayama_proteinId <- unique(x$proteinId)
  
  ABV2_transcriptId <- CopciABSmallExonNoNA_dict %>% filter(protID == Okayama_proteinId) %>% .$transcript_id
  
  ab <- longread_annot_T0_df %>% filter(type == "exon", transcript_id==ABV2_transcriptId) %>% 
    filter(exon_id %in% smallESumDesc_tbl$exonID) 
  o <- okayama_annotgff3Exon_df %>% 
    filter(proteinId == Okayama_proteinId) %>% 
    arrange(ID) %>% 
    filter(!seq(n()) %in% c(1,n())) %>% 
    filter(width<=15)
  
  ABV2_data <- ab %>% 
    select(transcript_id, exon_id, width) %>% 
    setNames(c("ABV2_transcriptId",	"ABV2_microexonId",	"ABV2_microexonwidth"))
  
  O_data <- o %>% 
    ungroup() %>% 
    select(proteinId, ID, width) %>% 
    setNames(c("Okayama_proteinId",	"Okayama_microexonId",	"Okayama_microexonwidth"))
  
  if(nrow(O_data) == 0) O_data[1,] <- NA
  
  if(nrow(o) > 0) {
    ab_exon_seq <- BSgenome::getSeq(genome_seq, makeGRangesFromDataFrame(ab, keep.extra.columns = TRUE))
    ABV2_data$seq_nt_ABV2 <- as.character(ab_exon_seq)
    
    o_exon_seq <- BSgenome::getSeq(genomeO_seq, makeGRangesFromDataFrame(o, keep.extra.columns = TRUE))
    O_data$seq_nt_Okayama <- as.character(o_exon_seq)
    
  } else {
    ab_exon_seq <- BSgenome::getSeq(genome_seq, makeGRangesFromDataFrame(ab, keep.extra.columns = TRUE))
    ABV2_data$seq_nt_ABV2 <- as.character(ab_exon_seq)
    
    O_data$seq_nt_Okayama <- NA
  }

  result <- list(ABV2 = ABV2_data,
                 Okayama = O_data)
  
  return(result)
  
}

compMicroEx_lt <- lapply(okayama_annotgff3ExonMicroInAB_lt, FUN = compMicroEx_fct)

getMicroExonIdentity_fct <- function(x) {

  if(!all(is.na(x$Okayama))) {
    pid_vct <- vector("numeric")
    for(i in x$ABV2$seq_nt_ABV2) {
      s1 <- DNAString(i)
      for(j in x$Okayama$seq_nt_Okayama) {
        s2 <- DNAString(j)
        palign_AB_O <- pwalign::pairwiseAlignment(s1, s2, type="global")
        pid_vct <- append(pid_vct, pid(palign_AB_O, type = "PID4"))
      }
    }
    
    pname_vct <- vector("numeric")
    for(i in x$ABV2$ABV2_microexonId) {
      for(j in x$Okayama$Okayama_microexonId) {
        palign_AB_O_comp <- str_c(i, j, sep = "__")
        pname_vct <- append(pname_vct, palign_AB_O_comp)
      }
    }
    
    names(pid_vct) <- pname_vct
    
    pid_vct %>% as.data.frame() %>% 
      setNames("pident") %>% 
      rownames_to_column(var = "comps") %>% 
      separate(comps, sep = "__", into = c("AB", "O")) %>% 
      group_by(AB) %>% 
      filter(pident == max(pident)) %>% 
      filter(seq(n())==1) %>%
      ungroup() %>% 
      rename("max_pident" = "pident") -> max_pid_tbl

    joinedTable <- left_join(x$ABV2, max_pid_tbl, by=c("ABV2_microexonId"="AB")) %>% 
      left_join(x$Okayama, by=c("O"="Okayama_microexonId"))
    
    return(joinedTable)
  } else {
    noMatch_tbl <- cbind(x[[1]], x[[2]])
    noMatch_tbl %>% 
      mutate(O=NA, max_pident = NA) %>% 
      select(ABV2_transcriptId, ABV2_microexonId, ABV2_microexonwidth, seq_nt_ABV2, 
             O, max_pident, Okayama_proteinId,Okayama_microexonwidth, seq_nt_Okayama) -> noMatch_tbl
    return(noMatch_tbl)
  }
}

compMicroExwithidentity_lt <- lapply(compMicroEx_lt, FUN = getMicroExonIdentity_fct)

table(sapply(compMicroExwithidentity_lt, ncol))

compMicroExwithidentity_lt[sapply(compMicroExwithidentity_lt, function(x) {any(duplicated(x$O))})]

cleanDuplicatedO_fct <- function(x) {

  if(any(duplicated(x$O))) {
    x %>% 
    group_by(O) %>% 
    mutate(O = ifelse(max_pident != max(max_pident), NA, O),
           Okayama_proteinId = ifelse(max_pident != max(max_pident), NA, Okayama_proteinId),
           Okayama_microexonwidth = ifelse(max_pident != max(max_pident), NA, Okayama_microexonwidth),
           seq_nt_Okayama = ifelse(max_pident != max(max_pident), NA, seq_nt_Okayama),
           max_pident = ifelse(max_pident != max(max_pident), NA, max_pident)) %>% 
      ungroup() -> cleanX
    return(cleanX)
  } else {
    return(x)
  }
}

compMicroExwithidentityNoMultyComp_lt <- lapply(compMicroExwithidentity_lt, FUN=cleanDuplicatedO_fct)

compMicroExwithidentityNoMultyComp_lt[sapply(compMicroExwithidentityNoMultyComp_lt, function(x) {any(duplicated(na.exclude(x$O)))})][[1]] %>% View()

compMicroExwithidentityNoMultyComp_tbl <- bind_rows(compMicroExwithidentityNoMultyComp_lt, .id = "Okayama_proteinIdId")
compMicroExwithidentityNoMultyComp_tbl %>% 
  mutate(PresenceApsenceIndex = (max_pident == 100) & !is.na(max_pident)) -> compMicroExwithidentityNoMultyComp_tbl


table(compMicroExwithidentityNoMultyComp_tbl$max_pident == 100, useNA="ifany")

compMicroExwithidentityNoMultyComp_tbl %>% 
  dplyr::count(PresenceApsenceIndex) 

compMicroExwithidentityNoMultyComp_tbl %>% 
  group_by(PresenceApsenceIndex) %>% 
  summarise(nGenes = length(unique(ABV2_transcriptId)))

write_tsv(compMicroExwithidentityNoMultyComp_tbl, "microexonCopciABV2_vs_OkayamaCopm.tsv")

length(unique(compMicroExwithidentityNoMultyComp_tbl$ABV2_transcriptId)) 
dim(compMicroExwithidentityNoMultyComp_tbl) 


length(unique(smallESumDesc_tbl$geneID))
dim(smallESumDesc_tbl)

CopciABSmallExon_dict %>% 
  filter(is.na(target)) -> CopciABSmallExonNA_dict

length(CopciABSmallExonNA_dict$gene_id)
View(CopciABSmallExonNA_dict)


smallESumDesc_tbl %>% 
  filter(geneID %in% CopciABSmallExonNA_dict$gene_id) -> CopciABSmallExonNAtranscript_dict

nrow(CopciABSmallExonNAtranscript_dict)


temp <- filter(compMicroExwithidentityNoMultyComp_tbl, PresenceApsenceIndex)



toPlot_df <- tibble(l = c("l1", "l2", "l2","l3", "l3", "l3"),
                    gene = c("a","b", "c", "d", "e", "f"),
                    n = c(1580, 1403, 177, 862, 541, 177))

toPlot_df %>% 
  ggplot(aes(x=l, y=n, fill=gene)) +
  geom_bar(color = "white",
           stat = "identity",
           position = "stack",
           width = 1,
           lwd = 2,
           show.legend = FALSE) +
  scale_fill_manual(values = c("gray50", "#ABD9E9", "#FEE090", "#4575B4", "#F68A69", "rosybrown4")) +
  ylab("Number of internal microexons") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  

ggsave("CopciABCopciO_internalMicroexonComp.pdf", device="pdf", height = 5, width = 5)
