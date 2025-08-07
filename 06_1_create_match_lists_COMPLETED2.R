#' UKB-PPP goes by UniprotID and gene name
#' Match it to the proteome list and resolve synonyms manually
#' 
#' For proteins that are measured in more than 1 panels -> select one
#'    order it by OLINK ID and pick the first one as correlation between levels were reported to be good.
#' 
#' 

library(tidyverse)
library(data.table)
library(glue)
library(yaml)


args  <- read_yaml("/config/config.yaml")[['proj_resource']]

outdir  <- file.path(args$rt_dir,'results','GWAS_interactors')

gene_annotation  <- readRDS("/relevant_ensembl_canonical_protein_coding_wo_strand.Rds")



ukb_independent_qtls  <- read.table(
  "/UKB_PPP_independent_pQTLs.tsv",sep="\t",header=T,quote=""
) %>% dplyr::select(UKBPPP.ProteinID) %>% unique()


ukb_independent_qtls$protein_name  <- unlist(
  lapply(
    strsplit(ukb_independent_qtls$UKBPPP.ProteinID,":"),
    FUN = function(x){x[1]}
  )
)
ukb_independent_qtls$uniprot_id  <- unlist(
  lapply(
    strsplit(ukb_independent_qtls$UKBPPP.ProteinID,":"),
    FUN = function(x){x[2]}
  )
)
ukb_independent_qtls$olink_ID  <- unlist(
  lapply(
    strsplit(ukb_independent_qtls$UKBPPP.ProteinID,":"),
    FUN = function(x){x[3]}
  )
)

# keep only one per group
# these are proteins that were measured in more than 1 panels so
# same protein but 2 different uniprot IDs
# since will be querying by protein and uniprot IDs, so keep one entry
ukb_independent_qtls  <- ukb_independent_qtls %>% 
  group_by(protein_name,uniprot_id) %>% 
  arrange(olink_ID) %>% slice_max(olink_ID)

in_ref  <- ukb_independent_qtls %>% filter(protein_name %in% gene_annotation$gene_name)
not_in_ref  <- ukb_independent_qtls %>% filter(!protein_name %in% gene_annotation$gene_name)


selected_ukb_pqtl_file_w_anno  <- inner_join(
  in_ref,
  gene_annotation,
  by = c("protein_name"="gene_name")
)

assertthat::assert_that(
  nrow(selected_ukb_pqtl_file_w_anno)==nrow(in_ref)

)
fwrite(
  selected_ukb_pqtl_file_w_anno,
  "/UKB_PPP_matched_data.tsv",
  sep="\t",quote=F
)




fwrite(
  not_in_ref,
  "/UKB_PPP_matched_data_not_in_ref.tsv",
  sep="\t",quote=F
)



