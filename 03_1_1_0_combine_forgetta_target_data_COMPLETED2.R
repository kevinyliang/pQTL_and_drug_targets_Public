#' for each of the genes that was listed as positive controls (drugs or Mendelian disease genes)
#' 1. obtain the ensembl ID using the proteome gene list.
#'  -> Any genes without ensembl ID (hence were not in the gene list) are removed
library(data.table)
library(tidyverse)
library(glue)


drug_target_info  <- read.table(
  "/drug_targets.txt",sep="\t",header=T,quote=""
) %>% mutate(src='drug')
mendel_gene_info  <- read.table(
  "/mendelian_disease_genes.txt",sep="\t",header=T,quote=""
)%>%mutate(src = 'mendel')

gene_annotation  <- readRDS("/relevant_ensembl_canonical_protein_coding.Rds")

gene_id  <- setNames(
  gene_annotation$gene_id_novrs,nm=gene_annotation$gene_name
)


all_targets  <- rbind(
  drug_target_info,
  mendel_gene_info
)
all_targets  <- all_targets %>% 
  aggregate(
    by = src ~ Trait + Gene.Symbol,
    FUN = paste,collapse=','
  )
all_targets$ensembl_id  <- gene_id[all_targets$Gene.Symbol]
all_targets <- all_targets %>% filter(!is.na(ensembl_id))

saveRDS(
  all_targets,
  "/forgetta_target.Rds"
)