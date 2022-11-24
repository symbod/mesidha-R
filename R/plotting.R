# Plotting Functions

# Load required libraries ----
suppressPackageStartupMessages({
  required_packages <- c("data.table","ggplot2", "dplyr")
  for(package in required_packages){
    if(!require(package,character.only = TRUE, quietly = TRUE)) install.packages(package, dependencies = TRUE, quietly = TRUE)
    library(package, character.only = TRUE, quietly = TRUE)
  }
})


#' Overview plot for logging data of mqhandler's functions
#'
#' @param logging overview logging Dataframe that has been returned by one of the four methods of mqhandler
#' @param out_dir (optional) output directory to save the plots
#' @param file_type (optional) file type for the plots (png, pdf, jpg, ...)
#'
#' @return list of two ggplot objects
#' @export
#'
#' @examples
create_overview_plot <- function(logging, out_dir=NULL, file_type = "png"){
  dt <- logging %>% dplyr::select(starts_with("Nr"))
  melted_dt <- melt(dt, variable.name = "type", value.name = "nr", id.vars = NULL)
  boxplot <- ggplot(melted_dt, aes(x = type, y = nr, fill=type)) + geom_boxplot() + xlab("") + 
    ylab("Number of entries per line") + ggtitle("Distribution of the number of entries per line") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) + 
    scale_fill_manual(guide = "none", values = c("#2c6a94", "#dc7530", "#36863a", "#b63737")) +
    stat_boxplot(geom ='errorbar', width = 0.4)
  sum_dt <- colSums(dt)
  sum_dt <- data.frame("type" = names(sum_dt), "nr" = as.vector(sum_dt))
  sum_dt$type <- factor(sum_dt$type, levels=sum_dt$type)
  barplot <- ggplot(sum_dt, aes(x=type, y=nr, fill=type)) + geom_bar(stat="identity") +
    geom_text(aes(label = nr), vjust=-0.5) + xlab("") + ylab("Number of entries") +
    theme(axis.text.x = element_text(angle = 45, vjust=0.5)) + ggtitle("Number of entries added and removed") + 
    scale_fill_manual(guide = "none", values = c("#2c6a94", "#dc7530", "#36863a","#b63737"))
  if(!is.null(out_dir)){
    path <- paste0(out_dir,"overview_log_box.", file_type)
    ggsave(path, boxplot, dpi=80, width=6, height=6)
    path <- paste0(out_dir, "overview_log_bar.", file_type)
    ggsave(path, barplot, dpi=80, width=6, height=6)
  }
  return(list(boxplot, barplot))
}


#' Detailed plot for logging data of mqhandler's filter_protein_ids method
#'
#' @param logging detailed logging Dataframe that has been returned by the filtered_protein_ids method of mqhandler
#' @param organism organism that has been specified in the filter_protein_ids method of mqhandler
#' @param reviewed specified reviewed parameter of the filter_protein_ids method of mqhandler
#' @param decoy specified rev_con parameter of the filter_protein_ids method of mqhandler
#' @param out_dir (optional) output directory to save the plot
#' @param file_type (optional) file type for the plot
#'
#' @return ggplot object
#' @export
#'
#' @examples
create_filter_detailed_plot <- function(logging, organism, reviewed, decoy, out_dir=NULL, file_type="png"){
  organisms <- c("Homo Sapiens (Human)", "Rattus norvegicus (Rat)", "Mus musculus (Mouse)", "Oryctolagus cuniculus (Rabbit)")
  names(organisms) <- c("human", "rat", "mouse", "rabbit")
  nr <- c()
  names <- c()
  if(decoy == FALSE){
    nr_decoys <- nrow(logging[logging$Organism == "Decoy", ])
    nr <- c(nr, nr_decoys)
    names <- c(names, "Decoys")
    if(reviewed == TRUE){
      nr_unreviewed <- nrow(logging[(logging$Reviewed == "unreviewed") & (logging$Organism == organisms[organism]),])
      nr <- c(nr, nr_unreviewed)
      names <- c(names, "Unreviewed")
    }
    nr_not_found <- nrow(logging[logging$Organism == "Not found",])
    nr_wrong_organism <- nrow(logging[!logging$Organism %in% c("Not found", "Decoy", organisms[organism]),])
    nr <- c(nr, nr_not_found, nr_wrong_organism)
    names <- c(names, "Not found IDs", "Wrong organism")
    dt <- data.frame(type = names, nr = nr)
    p <- ggplot(dt, aes(x=type, y = nr, fill=type)) + geom_bar(stat="identity") +
      geom_text(aes(label = nr), vjust=-0.5) + xlab("") + ylab("Number of entries") +
      ggtitle("Number of entries removed by cause") + 
      scale_fill_manual(guide = "none", values = c("#2c6a94", "#dc7530", "#36863a", "#b63737"))
    if(!is.null(out_dir)){
      path <- paste0(out_dir,"detailed_log.", file_type)
      ggsave(path, p, dpi=80, width=6, height=6)
    }
    return(p)
  }
}

#' Detailed plot for logging data of mqhandler's reduce_genenames method
#'
#' @param logging detailed logging Dataframe that has been returned by the reduce_genenames method of mqhandler
#' @param out_dir (optional) output directory to save the plot
#' @param file_type (optional) file type for the plot
#'
#' @return ggplot object
#' @export
#'
#' @examples
create_reduced_detailed_plot <- function(logging, out_dir = NULL, file_type = "png"){
  nr_not_found <- nrow(logging[logging$`Reduced Gene Name` == "Not found", ])
  nr_different_name <- nrow(logging[logging$`Reduced Gene Name` != "Not found", ])
  nr <- c(nr_not_found, nr_different_name)
  names <- c("Not found IDs", "Different Name")
  dt <- data.frame(type=names, nr=nr)  
  dt$type <- factor(dt$type, levels=dt$type)
  p <- ggplot(dt, aes(x=type, y = nr, fill=type)) + geom_bar(stat="identity") +
    geom_text(aes(label = nr), vjust=-0.5) + xlab("") + ylab("Number of entries") +
    ggtitle("Number of entries removed by cause") + 
    scale_fill_manual(guide = "none", values = c("#2c6a94", "#dc7530", "#36863a", "#b63737"))
  if(!is.null(out_dir)){
    path <- paste0(out_dir,"detailed_log.", file_type)
    ggsave(path, p, dpi=80, width=6, height=6)
  }
  return(p)
}
  
#' Detailed plot for logging data of mqhandler's map_orthologs method
#'
#' @param logging detailed logging Dataframe that has been returned by the map_orthologs method of mqhandler
#' @param organism organism that has been specified in the map_orthologs method of mqhandler
#' @param out_dir (optional) output directory to save the plot
#' @param file_type (optional) file type for the plot
#'
#' @return
#' @export
#'
#' @examples
create_ortholog_detailed_plot <- function(logging, organism, out_dir = NULL, file_type = "png"){
  logging[is.na(logging)] <- ""
  nr_target_names <- nrow(logging[(logging$ortholog_ensg != "") & (logging$target_symbol == ""),])
  nr_target_ids <- nrow(logging[(logging$ortholog_ensg == "") & (logging$target_symbol == ""),])
  nr_source_ids <- nrow(logging[logging$ortholog_ensg == "Not found",])
  nr_wrong_organism <- nrow(logging[!logging$source_organism %in% c("Not found", organism),])
  nr <- c(nr_target_names, nr_target_ids, nr_source_ids, nr_wrong_organism)
  names <- c("Not found target names", "Not found target IDs", "Not found source IDs", "Wrong organism")
  dt <- data.frame(type = names, nr = nr)
  dt$type <- factor(dt$type, levels = dt$type)
  p <- ggplot(dt, aes(x=type, y = nr, fill=type)) + geom_bar(stat="identity") +
    geom_text(aes(label = nr), vjust=-0.5) + xlab("") + ylab("Number of entries") +
    ggtitle("Number of entries removed by cause") + 
    scale_fill_manual(guide = "none", values = c("#2c6a94", "#dc7530", "#36863a", "#b63737")) +
    theme(axis.text.x = element_text(angle = 45, vjust=0.5))
  if(!is.null(out_dir)){
    path <- paste0(out_dir,"detailed_log.", file_type)
    ggsave(path, p, dpi=80, width=6, height=6)
  }
  return(p)
}
