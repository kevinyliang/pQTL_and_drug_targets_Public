#' Create a list of ~20K genes with ensembl and entrez ID so can map to various files
#' 
#' Filter criteria:
#'  1. transcripts only
#'  2. protein coding (transcript and gene)
#'  3. fully remapped to GRCh37
#'  4. omit chromosome Y and M
#'  5. Canonical transcripts (Ensembl_canonical only)
#'  6. Break ties with longest if exists
#' 
library(tidyverse)
library(rtracklayer)
library(data.table)
library(glue)
library(yaml)
library(biomaRt)


args  <- read_yaml("/config/config.yaml")[['proj_resource']]

outdir  <- "/results/GWAS_interactors"

#%% read in the gene code annotation file and find longest transcript from HAVANA 
#' Obtained from: https://www.gencodegenes.org/human/release_39lift37.html
gene_transcript_annotations  <- rtracklayer::readGFF(
  args$gene_annotation_data
)

#' Clean up the IDs
gene_transcript_annotations$gene_id_novrs  <- gsub("\\..*$","",gene_transcript_annotations$gene_id)
gene_transcript_annotations$gene_id_vrs  <- gsub("^.*?\\.","",gene_transcript_annotations$gene_id)
gene_transcript_annotations$transcript_id_novrs  <- gsub("\\..*$","",gene_transcript_annotations$transcript_id)


#%% obtain
# obtain indicies for entries we need
ensembl_canonical_protein_coding <- gene_transcript_annotations[
  which(
    gene_transcript_annotations$type == 'transcript' &
    gene_transcript_annotations$transcript_type == 'protein_coding' &
    gene_transcript_annotations$gene_type == 'protein_coding' &
    gene_transcript_annotations$remap_status == 'full_contig' &
    !grepl(as.character(gene_transcript_annotations$seqid),pattern='chrY|chrM') &
    unlist(
      lapply(
        gene_transcript_annotations$tag, # tag is a list of characters
        FUN=function(x){any(grepl("Ensembl_canonical",x))} 
      )
    )
  ),
]

# retain only information we will use
relevant_ensembl_canonical_protein_coding = ensembl_canonical_protein_coding[,c("seqid","type","start","end","ID","gene_id","gene_id_novrs","gene_type","gene_name","transcript_type","transcript_name",'strand')] %>% as.data.frame()

relevant_ensembl_canonical_protein_coding$seqid <- as.character(relevant_ensembl_canonical_protein_coding$seqid)

# recode chrX as chr23 (the 23) to match with other sources of data
relevant_ensembl_canonical_protein_coding$seqid  <- ifelse(
  relevant_ensembl_canonical_protein_coding$seqid == 'chrX','chr23',relevant_ensembl_canonical_protein_coding$seqid
)



#%% Obtain transcript length
relevant_ensembl_canonical_protein_coding$ts_length  <- relevant_ensembl_canonical_protein_coding$end - relevant_ensembl_canonical_protein_coding$start 
assertthat::assert_that(all(relevant_ensembl_canonical_protein_coding$ts_length > 0))

# if there are duplicated entries based on criteria we set and information we consider
# we will obtain the longest one
dup_genes  <- relevant_ensembl_canonical_protein_coding$gene_name[duplicated(relevant_ensembl_canonical_protein_coding$gene_name)]
length(dup_genes)

relevant_ensembl_canonical_protein_coding <- relevant_ensembl_canonical_protein_coding %>% group_by(gene_name) %>% slice_max(ts_length) %>% ungroup()

# make sure no duplicates
assertthat::assert_that(
  sum(duplicated(relevant_ensembl_canonical_protein_coding$gene_id_novrs)) == 0 &
  sum(duplicated(relevant_ensembl_canonical_protein_coding$gene_id)) == 0 &
  sum(duplicated(relevant_ensembl_canonical_protein_coding$gene_name)) == 0 &
  sum(is.na(relevant_ensembl_canonical_protein_coding$gene_id)) == 0 &
  sum(is.na(relevant_ensembl_canonical_protein_coding$gene_name)) == 0
)


saveRDS(
  relevant_ensembl_canonical_protein_coding,
  file.path(
    outdir,'relevant_ensembl_canonical_protein_coding.Rds'
  )
)
fwrite(
  relevant_ensembl_canonical_protein_coding,
  file.path(outdir,'relevant_ensembl_canonical_protein_coding.tsv'),
  sep="\t",quote=F
)

#%%

