#' Visualize and select an optimal number of clusters of individuals based
#' on the network information for a pathway of interest
#'
#' This function fertemines and visualizes the optimal number of clusters using
#' kmeans clustering
#' 
#' @import factoextra
#' 
#' @param pathway_of_interest List with genes in a pathway of interest
#' @param reg_net Numeric matrix with samples in columns, features in rows
#' @param edges Data frame containing information on "reg" and "tar"
#' @param scale_data Whether to scale the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' @param center_data Whether to scale the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' @param kmax Integer specifying the maximum nuber of clusters to consider
#' 
#' @return ggplot object
#' @export


select_number_clusters <- function(pathway_of_interest,
                    reg_net,
                    edges,
                    scale_data = TRUE,
                    center_data = TRUE
                    kmax = 8) {
    results <- get_pathway_ind_scores(pathway_of_interest,
                    reg_net,
                    edges,
                    scale_data = scale_data,
                    center_data = center_data)
    g <- factoextra::fviz_nbclust(results,
                    FUNcluster = kmeans,
                    k.max = kmax)
    plot(g)
}

#' Visualize clustering results
#'
#' Visuzalition of clustering of patients into specified number of clusters
#' 
#' @import factoextra
#' 
#' @param pathway_of_interest List with genes in a pathway of interest
#' @param reg_net Numeric matrix with samples in columns, features in rows
#' @param edges Data frame containing information on "reg" and "tar"
#' @param scale_data Whether to scale the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' @param center_data Whether to scale the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' @param number_of_clusters Integer, specifying the number of clusters
#' 
#' @return ggplot object
#' @export 
#' 
visualize_clusters <- function(pathway_of_interest,
                    reg_net,
                    edges,
                    scale_data = TRUE,
                    center_data = TRUE,
                    number_of_clusters) {
    if (!is.integer(number_of_clusters) || number_of_clusters <= 0 ||
        number_of_clusters >= ncol(reg_net)) {
        stop("number_of_clusters must be a positive integer
        less than the number of patients.")
    }
    results <- get_pathway_ind_scores(pathway_of_interest,
                    reg_net,
                    edges,
                    scale_data = scale_data,
                    center_data = center_data)
    final_clusters <- kmeans(results, number_of_clusters, nstart = 25)
    g <- factoextra::fviz_cluster(final_clusters,
                    data = results,
                    geom = "point",
                    ggtheme = theme_bw(),
                    ellipse = TRUE)
    plot(g)
 }