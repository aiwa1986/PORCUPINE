#' Extract patient heterogeneity scores on the fist two principal components
#' 
#' This function extracts patient heterogeneity scores on the fist
#' two principal components.
#'
#' @import irlba, factoextra
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


#' PCA plot of individuals colored according to a clinical feature
#' 
#' This function makes a PCA plot of individuals colored according to a
#' clinical feature
#' 
#' @import ggplot2
#' 
#' @param individual_scores individual heterogeneity scores on PC1 and PC2 
#' @param clin_feature clinical feature for individuals (e.g. age)
#' 
#' @return PCA plot 
#' @export
#' 

plot_clinical_association <- function(individual_scores, clin_col) {
    g <- ggplot(individual_scores, aes(scores[, 1], scores[, 2], 
        col = clin_feature)) + geom_point(size = 2.5) +
        xlab("PC1") + ylab("PC2") + theme_bw()
    plot(g)
}

#' Significance of association of clinical feature with patient heterogeneity
#' scores (PC1 or PC2 scores)
#' 
#' This function performs a correlation test or Kruskal-Wallis test (for numeric
#' or categorical variables, respectively) between patient heterogeneity scores 
#' on PC1 or PC2

#' @import stats
#' 
#' @param individual_scores PC1 and PC2 heterogeneity scores for individuals
#' @param component principal component to look at (1,2)
#' @param clin_feature clinical feature for individuals (e.g. age)
#' 
#' @return pvalue
#' @export


clinical_correlation <- function(individual_scores, component, clin_feature) {
    if (is.numeric(clin_feature)) {
    # Perform correlation test if the column is numeric
    pval <- 
        stats::cor.test(individual_scores[, component], clin_feature)$"p.value"
        } else {
    # Perform Kruskal-Wallis test if the column is not numeric
    pval <- 
        stats::kruskal.test(individual_scores[, component] ~ clin_feature)$"p.value"
        }
    return(pval)
}
