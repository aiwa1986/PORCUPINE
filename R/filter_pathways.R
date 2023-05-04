#' Filter a pathway list
#'
#' This function filters a list of pathways to include only genes in pathways
#' present in networks
#' 
#' @import plyr
#' @import dplyr
#' 
#' @param pathways_list list of pathways
#' @param edges Table, containing information on "reg" and "tar"
#' 
#' @return A list of filtered pathways
#' 
#' @export

filter_pathways <- function(pathways_list,
                    edges) {
  pathways_filt <- plyr::ldply(pathways_list, data.frame) %>%
                   dplyr::rename(pathway = ".id", gene = "X..i..") %>%
                   dplyr::filter(gene %in% edges$tar) %>%
                   with(., split(gene, pathway))
  return(pathways_filt)
}

#' Filter a pathway list based on a number of genes in a pathway
#'
#' This function filters a list of pathways based on specified minimum and 
#' maximum size for number of genes in a pathway
#' 
#' @import purrr
#' 
#' @param pathways_list list of pathways
#' @param minSize Minimum size for number of genes in a pathway (default: 5)
#' @param maxSize Maximum size for a number of genes in a pathway (default: 150)
#' 
#' @return A list of filtered pathways
#' 
#' @export

filter_pathways_size <- function(
                        pathways_list,
                        minSize = 5,
                        maxSize = 150) {
    pathways_filt <- purrr::keep(pathways_list, function(x) length(x) >= minSize)
    pathways_filt <- purrr::keep(pathways_filt, function(x) length(x) <= maxSize)
    return(pathways_filt)
}
