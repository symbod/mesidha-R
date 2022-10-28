#!/usr/bin/env Rscript

# Load required libraries ----
## Load R libraries ----
suppressPackageStartupMessages({
  required_packages <- c("reticulate","data.table", "roxygen2", "devtools")
  for(package in required_packages){
    if(!require(package,character.only = TRUE, quietly = TRUE)) install.packages(package, dependencies = TRUE, quietly = TRUE)
    library(package, character.only = TRUE, quietly = TRUE)
  }
})
## Load python library ----
py_install("mqhandler", pip = TRUE)
if (!py_module_available("mqhandler")) {
  stop("Python package could not be installed")
}
mq <- import("mqhandler")


# roxygenise();      # Builds the help files

## Main functions ----

### Filter protein IDs ----

#' @title Filter Protein IDs
#' 
#' Remove decoy IDs, contaminated IDs and / or filter IDs by organism.
#'
#' @param data a Dataframe with a column containing the protein IDs
#' @param protein_column name of column with protein IDs
#' @param organism (optional) specify organism the ids should match to
#' @param decoy (optional) bool to indicate if protein ids from decoy fasta (REV__, CON__) should be kept
#' @param action (optional) what to do, if IDs cell is empty after filtering. Keep empty cell, delete it or fill it based on gene name.
#' @param gene_column (optional) name of column with gene names
#' @param reviewed (optional) bool to indicate if newly retrieved protein IDs should be reduced to reviewed ones
#'
#' @return filtered Dataframe
#' @export
#'
#' @examples
filter_protein_ids <- function(data, protein_column, organism = NULL, decoy = FALSE, 
                               action = "delete", gene_column = NULL, reviewed = TRUE) {
  
  return(mq$filter_ids$filter_protein_ids(data = data, id_column = protein_column, 
                                          organism = organism, decoy = decoy, 
                                          action = action, gene_column = gene_column,
                                          reviewed = reviewed))
}

## Remap Gene Names

#' @title Remap Gene Names
#' 
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
#' @param skip_filled (optional) bool to indicate if already filled gene names should be skipped
#' @param organism (optional) specify organism the ids should match to
#' @param fasta (optional) fasta file when mode all or fasta
#'
#' @return remapped Dataframe
#' @export
#'
#' @examples
remap_genenames <- function(data, mode, protein_column, gene_column = NULL,
                            skip_filled = FALSE, organism = NULL, fasta = NULL){
  return(mq$remap_genenames$remap_genenames(data = data, mode = mode, protein_column = protein_column,
                                            gene_column = gene_column, skip_filled = skip_filled, 
                                            organism = organism, fasta = fasta))
}

### Get Orthologs

#' @title Get Orthologs
#' 
#' Get ortholog gene names from origin organism to target organism.
#'
#' @param data a Dataframe with a column containing the gene names
#' @param gene_column name of column with gene names
#' @param organism specify organism the ids currently match to
#' @param tar_organism specify organism the ids should match to
#'
#' @return dataframe with orthologs
#' @export
#'
#' @examples
get_orthologs <- function(data, gene_column, organism, tar_organism) {
  return(mq$map_orthologs$get_orthologs(data = data, gene_column = gene_column, 
                                        organism = source_organism, tar_organism = target_organism))
}

## Smaller functions ----

#' @title Grep Header Information
#'
#' Grep the information,that is stored in the headers inside a 
#' fasta file, and return them inside a Dataframe. 
#'
#' @param fasta the path to the fasta file
#'
#' @return a Dataframe with the collected information
#' @export
#'
#' @examples
grep_header_info <- function(fasta) {
  return(mq$fasta_grepper$grep_header_info(fasta))
}
