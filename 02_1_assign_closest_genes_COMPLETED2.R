#' For each set of predefined traits -> assigned to it the closest gene using the list of gene that was curated in step 1
#' 
#' 
library(tidyverse)
library(glue)
library(data.table)
library(yaml)

args  <- read_yaml("/config/config.yaml")[['proj_resource']]

source("./x1_find_nearest_gwas_genes.R")

outdir  <- file.path('/nearest_genes')



trait_to_process  <- args$ei_12_traits

# this is with strand information
protein_annotation_file <- readRDS("/relevant_ensembl_canonical_protein_coding_tmp.Rds")

gwas_clumped_dir <- "/clumped_GWAS"

for (each_trait in trait_to_process){
  outfile  <- file.path(
      outdir,glue("{each_trait}_GWAS_hit_nearest_gene.Rds")
    )
  print(each_trait)
  trait_clumped_res  <- readRDS(file.path(gwas_clumped_dir,glue("{each_trait}_processed.clumped_gwas.Rds")))
  trait_closest_gene_info  <- find_nearest_gwas_genes(
    trait_clumped_res = trait_clumped_res,
    protein_annotation_file = protein_annotation_file
  )
  saveRDS(  
    trait_closest_gene_info,
    outfile
  )
}