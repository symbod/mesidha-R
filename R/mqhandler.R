#!/usr/bin/env Rscript

# Load required libraries ----
## Load R libraries ----

.onLoad <- function(libname, pkgname) {
  library("reticulate", character.only = TRUE, quietly = TRUE)
  reticulate::configure_environment(pkgname)
  ## Load python library ----
  py_install("mqhandler", pip = TRUE, ignore_installed=TRUE, python_version="0.0.22")
}

# Main functions ----

## Filter protein IDs ----

#' @title Filter Protein IDs
#' @description 
#' Remove decoy IDs, contaminated IDs and / or filter IDs by organism.
#' 
#' @param data a Dataframe with a column containing the protein IDs
#' @param protein_column name of column with protein IDs
#' @param organism (optional) specify organism the ids should match to
#' @param rev_con (optional) bool to indicate if protein ids from decoy fasta (REV__, CON__) should be kept
#' @param res_column name of column of filtered protein IDs. If NULL, the protein_column will be overridden
#' @param keep_empty bool to indicate if empty reduced gene names cells should be kept or deleted
#' @param reviewed (optional) bool to indicate if newly retrieved protein IDs should be reduced to reviewed ones
#'
#' @return filtered Dataframe
#' @export
filter_protein_ids <- function(data, protein_column, organism = NULL, rev_con = FALSE, 
                               keep_empty = FALSE, res_column = NULL, reviewed = TRUE) {
  mq <- import("mqhandler")
  return(mq$filter_ids$filter_protein_ids(data = data, protein_column = protein_column, 
                                          organism = organism, rev_con = rev_con, 
                                          keep_empty = keep_empty, res_column = res_column,
                                          reviewed = reviewed))
}

## Remap Gene Names ----

#' @title Remap Gene Names
#' @description 
#' Use protein IDs to get associated Gene Names and fill exmpty entries 
#' and optionally replace already existing Names. Following modes are possible.
#' 
#' Modes:
#' all - Use primarly fasta infos and additionally uniprot infos.
#' fasta - Use information extracted from fasta headers.
#' uniprot - Use mapping information from uniprot and use all gene names.
#' uniprot_primary - Use mapping information from uniprot and only all primary gene names.
#' uniprot_one - Use mapping information from uniprot and only use most frequent single gene name.
#'
#' @param data a Dataframe with a column containing the protein IDs
#' @param mode mode of refilling (See in description for more infos)
#' @param protein_column name of column with protein IDs
#' @param gene_column (optional) name of column with gene names
#' @param res_column name of column of reduced gene names results. If NULL, the gene_column will be overridden
#' @param keep_empty bool to indicate if empty reduced gene names cells should be kept or deleted
#' @param organism (optional) specify organism the ids should match to
#' @param fasta (optional) fasta file when mode all or fasta
#'
#' @return remapped Dataframe
#' @export
remap_genenames <- function(data, mode, protein_column, gene_column, res_column = NULL,
                            skip_filled = FALSE, organism = NULL, fasta = NULL, keep_empty = FALSE){
  mq <- import("mqhandler")
  return(mq$remap_genenames$remap_genenames(data = data, mode = mode, protein_column = protein_column,
                                            gene_column = gene_column, res_column = res_column,
                                            skip_filled = skip_filled, organism = organism,
                                            fasta = fasta, keep_empty = keep_empty))
}

## Reduce Gene Names ----

#' @title Reduce Gene Names
#' @description 
#' Use existing Gene Names to reduce them and replace them by a single Name 
#' to redundant Names. Following modes are possible.
#' 
#' Modes:
#' ensembl - Use gProfiler to reduce gene names to those having a Ensembl ID
#' HGNC - Use HGNC database to reduce gene names to those having an entry in HGNC (only for human)
#' mygeneinfo - Use mygeneinfo database to reduce gene names to those having an entry in mygeneinfo
#' enrichment - Use gProfiler to reduce gene names to those having a functional annotation
#'
#' @param data a Dataframe with a column containing the protein IDs
#' @param mode mode of refilling (See in description for more infos)
#' @param gene_column name of column with gene names
#' @param res_column name of column of reduced gene names results. If NULL, the gene_column will be overridden
#' @param keep_empty bool to indicate if empty reduced gene names cells should be kept or deleted
#' @param organism (optional) specify organism the ids should match to
#' @param HGNC_mode setting on how to selected the gene names in HGNC (mostfrequent, all)
#'
#' @return reduced Dataframe
#' @export
reduce_genenames <- function(data, mode, gene_column, organism, 
                            res_column = NULL, keep_empty = FALSE, HGNC_mode = "mostfrequent"){
  mq <- import("mqhandler")
  return(mq$reduce_genenames$reduce_genenames(data = data, mode = mode, gene_column = gene_column,
                                              res_column = res_column, keep_empty = keep_empty, 
                                              organism = organism, HGNC_mode = HGNC_mode))
}

## Map Orthologs ----

#' @title Map Orthologs
#' @description 
#' Get ortholog gene names from origin organism to target organism.
#'
#' @param data a Dataframe with a column containing the gene names
#' @param gene_column name of column with gene names
#' @param organism specify organism the ids currently match to
#' @param tar_organism specify organism the ids should match to
#' @param res_column name of column of reduced gene names results. If NULL, the gene_column will be overridden
#' @param keep_empty bool to indicate if empty reduced gene names cells should be kept or deleted
#'
#' @return dataframe with orthologs
#' @export
map_orthologs <- function(data, gene_column, organism, tar_organism,
                          res_column, keep_empty) {
  mq <- import("mqhandler")
  return(mq$map_orthologs$map_orthologs(data = data, gene_column = gene_column, 
                                        organism = organism, tar_organism = tar_organism,
                                        res_column = res_column, keep_empty = keep_empty))
}

# Smaller functions ----

#' @title Grep Header Information
#' @description 
#' Grep the information,that is stored in the headers inside a 
#' fasta file, and return them inside a Dataframe. 
#'
#' @param fasta the path to the fasta file
#'
#' @return a Dataframe with the collected information
#' @export
grep_header_info <- function(fasta) {
  mq <- import("mqhandler")
  return(mq$fasta_grepper$grep_header_info(fasta))
}