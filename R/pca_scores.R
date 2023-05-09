#' Extract patient heterogeneity scores on the fist two principal components
#' 
#' This function extracts patient heterogeneity scores on the fist
#' two principal components.
#' 
#' @param data Numeric matrix with samples in columns, and features in rows
#' @param scale_data Whether to scale the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' @param center_data Logical, whether to center the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' 
#' @return patient heterogeneity scores for PC1 and PC2
#' @export

get_ind_scores <- function(data, scale_data = TRUE, center_data = TRUE) {
    data_t <- Matrix::t(data)
    # Perform scaling and centering of the data
    res_pca <- irlba::prcomp_irlba(data_t,
                scale. = scale_data,
                center = center_data)
    individual_scores <- factoextra::get_pca_ind(res_pca)$coord
    individual_scores <- individual_scores[, 1:2]
    rownames(individual_scores) <- colnames(data)
    return(individual_scores)
}

#' Extract patient heterogeneity scores on the fist two principal components
#' for a list of pathways
#' 
#' This function extracts patient heterogeneity scores on the fist
#' two principal components for a list of pathways
#' 
#' @param pathways_list list of pathways
#' @param reg_net Numeric matrix with samples in columns, features in rows
#' @param edges Table, containing information on "reg" and "tar"
#' @param scale_data Whether to scale the data (TRUE) or not (FALSE).
#' Default is TRUE.
#' @param center_data Logical, whether to center the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' 
#' @return patient heterogeneity scores for each pathway
#' @export
#' 
get_pathway_ind_scores <- function(pathways_list,
                        reg_net,
                        edges,
                        scale_data = TRUE,
                        center_data = TRUE) {
    res <- lapply(pathways_list, function(pathway) {
        idx <- which(edges$tar %in% pathway)
        subnet <- reg_net[idx, ]
        individual_scores <- get_ind_scores(subnet,
                            scale_data = scale_data,
                            center_data = center_data)
    })
    names(res) <- names(pathways_list)
    return(res)
}


#' Extract feature scores on the first two principal components
#' 
#' This function extracts feature scores on the first two principal components
#' 
#' @import irlba
#' 
#' @param data Numeric matrix with samples in columns, and features in rows
#' @param scale_data Whether to scale the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' @param center_data Logical, whether to center the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' 
#' @return features scores on the PC1 and PC2
#' @export
#' 
get_features_scores <- function(data,
                    scale_data = TRUE,
                    center_data = TRUE) {
    data_t <- Matrix::t(data)
    res_pca <- irlba::prcomp_irlba(data_t,
                scale. = scale_data,
                center = center_data)
    edges_scores <- res_pca$rotation
    edges_scores <- as.data.frame(edges_scores[, 1:2])
    return(edges_scores)
}

#' Extract feature scores on the first two principal components for a list of
#' pathways
#' 
#' This functions extracts feature scores on the first two principal
#' components for a list of pathways
#' 
#' @param pathways_list list of pathways
#' @param reg_net Numeric matrix with samples in columns, features in rows
#' @param edges Table, containing information on "reg" and "tar"
#' @param scale_data Whether to scale the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' @param center_data Logical, whether to center the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' 
#' @return edges scores on PC1 and PC2 for each pathway
#' @export
#' 

get_pathway_features_scores <- function(pathways_list,
                            reg_net,
                            edges,
                            scale_data = TRUE,
                            center_data = TRUE) {
    res <- lapply(pathways_list, function(pathway) {
    idx <- which(edges$tar %in% pathway)
    reg  <- edges$reg[edges$tar %in% pathway]
    tar  <- edges$tar[edges$tar %in% pathway]
    subnet <- reg_net[idx, ]
    edges_scores <- get_features_scores(subnet,
                    scale_data = scale_data,
                    center_data = center_data)
    edges_scores$reg <- reg
    edges_scores$tar <- tar
    return(edges_scores)
    })
    names(res) <- names(pathways_list)
    return(res)
}