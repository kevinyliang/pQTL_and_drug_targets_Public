#' find_nearest_gwas_genes
#' Find the nearest genes to some presepcified SNPs
#' 
#' 1. For each SNP, obtain the position and chromosome
#' 2. from the Proteome gene list -> retain only those in the same chromosome
#' 3. Compute the distance to SNP (abs(start of Transcript - snp position))
#' 4. FIlter for those with distance less than or equal to 1,000,000
#' 5. Obtain the smallest distance (if tie, keep both)
find_nearest_gwas_genes <- function(trait_clumped_gwas_res,protein_annotation_file){
  closest_gene_info  <- data.frame()
  for (each_snp in trait_clumped_gwas_res$SNP){
    snp_info  <- unlist(strsplit(each_snp,":"))
    snp_chr  <- snp_info[1]
    snp_pos  <- as.numeric(snp_info[2])

    closest_gene  <- protein_annotation_file %>% filter(seqid == glue("chr{snp_chr}"))
    closest_gene  <- closest_gene %>% mutate(
      TSS = ifelse(strand == "+", start, end)
    )
    closest_gene  <- closest_gene %>% 
      mutate(
        closest_tss = abs(TSS - snp_pos)
      ) %>% filter(closest_tss <= 1000000) %>% slice_min(closest_tss)
    
    closest_gene$SNP  <- each_snp
    closest_gene$snp_chr  <- snp_chr
    closest_gene$snp_pos  <- snp_pos
    closest_gene_info  <- rbind(
      closest_gene_info,
      closest_gene
    )
  }
  return(closest_gene_info)
}