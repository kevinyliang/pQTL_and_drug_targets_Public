#' Read in the opentarget molecule list for all the drug - drug target pairs.
#'  1. read in all the openTarget info
#'  2. unnest the linked Targets part where each drugs can have more than 1 entry, 1 for each target
#'  3. Filter to those in the gene annotation file
#' 
#' 
#'  tables (https://community.opentargets.org/t/get-marketed-drugs-for-a-set-of-targets-with-bigquery/101): 
#'    molecules = drugs
#'    linkedtargets = what the drugs act on
#' linkedTargets: The molecule dataset comprises in linkedTargets.rows all the targets (in the form of Ensembl gene IDs) on which a drug acts according to its mechanism of action
#' 

library(tidyverse)
library(glue)
library(data.table)
library(yaml)
library(rtracklayer)
library(sparklyr)
library(purrr)
library(sparklyr.nested)

args  <- read_yaml("/config/config.yaml")[['proj_resource']]


outdir  <- "/GWAS_interactors"


gene_annotation_file  <- readRDS("/relevant_ensembl_canonical_protein_coding.Rds")

# obtain all the information

sc <- spark_connect(master = "local")
drug_info <- spark_read_parquet(
  sc = sc, 
  name = 'opentargets_drug',
  path = "/opentarget_info/drugs_vrs24.03/molecule/", 
  memory = FALSE
) %>% filter(!is.na(linkedTargets)) 
rel_drug_info  <- drug_info %>% 
  sdf_select(
    id,name,tradeNames,description,linkedDiseases,linkedTargets
  ) %>% sdf_unnest(
    linkedDiseases
  ) %>% dplyr::rename(
    disease = rows,
    disease_count = count
  ) %>% sdf_unnest(
    linkedTargets
  ) %>% dplyr::rename(
    target = rows,
    target_count = count
  ) %>% collect()
# target is a list of lists -> unnest twice to make each row a single drug target pair
rel_drug_info  <- rel_drug_info %>% filter(target_count > 0)
rel_drug_info  <- rel_drug_info %>% unnest(target)
rel_drug_info  <- rel_drug_info %>% unnest(target)

rel_drug_info  <- rel_drug_info %>% filter(target %in% gene_annotation_file$gene_id_novrs)


saveRDS(
  rel_drug_info,
  file.path(outdir,'opentarget_drug_target_info.Rds')
)



spark_disconnect(sc)