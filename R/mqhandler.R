#!/usr/bin/env Rscript

#' Main method used for running mqhandler
#'
#' Runs KeyPathwayMiner localy via the standalone
#' or remotely via the RESTful API of the KeyPathwayMiner
#' website. Input parameters are indicator matrices and a graph_file.
#' You can view or change the run parameters through the
#' kpm_options() function.
#'
#' @param indicator_matrices List of data frames.
#' @param graph Path of the graph file or an igraph object.
#' NULL if you want to use a graph from the web service (only for remote runs).
#' Use get_networks()
#' to see all networks.
#' @export

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

#' Filter Protein IDs
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

## Smaller functions ----

#' Grep Header Information
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
