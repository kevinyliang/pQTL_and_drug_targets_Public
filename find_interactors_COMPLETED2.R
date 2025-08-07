#' Find interactors for a given set of gene list
#' 
#' 
library(tidyverse)
library(glue)
library(data.table)
library(yaml)
library(sparklyr)
library(sparklyr.nested)



gene_annotation  <- readRDS("/relevant_ensembl_canonical_protein_coding.Rds")

find_interactors  <- function(gene_list){
  sc <- spark_connect(master = "local")
  intact_interactions  <- spark_read_parquet(
    sc = sc,
    name = 'intact_interactions',
    path = "/opentarget_info/interactions_vrs24.03/interaction",memory=F
  ) 


  intact_interaction_evidence  <- spark_read_parquet(
    sc = sc,
    name = 'intact_interactions_evidence',
    path = "/opentarget_info/interactions_vrs24.03/interactionEvidence",memory=F
  )
  rel_intact_interaction_evidence  <- intact_interaction_evidence %>% 
    sdf_unnest(speciesB) %>%  sdf_select(-taxonId,-mnemonic) %>% dplyr::rename(
      speciesB = scientificName
    ) %>%
    sdf_unnest(speciesA) %>%  sdf_select(-taxonId,-mnemonic) %>% dplyr::rename(
      speciesA = scientificName
    ) %>% 
    sdf_unnest(interactionResources)
  
  assertthat::assert_that(
    class(gene_list) == 'data.frame' & 
    "ensembl_id" %in% colnames(gene_list)
  )
  interactors_df  <- data.frame()
  for (r in 1:nrow(gene_list)){
    print(glue("{r} / {nrow(gene_list)}"))
    row  <- gene_list[r,]
    assertthat::assert_that(
      class(row$ensembl_id)=='character'
    )
    n_interactors  <- intact_interactions %>% 
      filter(
        targetA == row$ensembl_id &
        sourceDatabase == "intact" &
        targetA != targetB
      )  %>%
      sdf_select(
        targetB
      ) %>% collect() %>% filter(targetB %in% gene_annotation$gene_id_novrs) %>% unique()

    n_filtered_interactors  <- intact_interactions %>% 
      filter(
        targetA == row$ensembl_id &
        sourceDatabase == "intact" &
        targetA != targetB &
        scoring > args$intact_mi_threshold
      )  %>%
      sdf_select(
        targetB
      ) %>% collect() %>% filter(targetB %in% gene_annotation$gene_id_novrs) %>% unique()

    direct_interactors  <- rel_intact_interaction_evidence %>% 
      filter(
        targetA == row$ensembl_id &
        sourceDatabase == 'intact' &
        interactionTypeShortName == 'direct interaction' &
        targetA != targetB
      ) %>% 
      sdf_select(
        targetB
      ) %>% collect()  %>% filter(targetB %in% gene_annotation$gene_id_novrs) %>% unique()
    assertthat::assert_that(
      all(direct_interactors$targetB %in% n_interactors$targetB)
    )
    row$n_direct  <- nrow(direct_interactors)
    row$n_interactors  <- nrow(n_interactors)
    row$n_filtered  <- nrow(n_filtered_interactors)
    interactors_df  <- rbind(
      interactors_df,
      row
    )
  }
  spark_disconnect(sc)
  return(interactors_df)
}
