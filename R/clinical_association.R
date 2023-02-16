#' Extract individual scores from the PCA results
#'
#' @import irlba
#' 
#' @param data Numeric matrix with samples in columns, and features in rows
#' 
#' @return individual scores for PC1 to PC10
#' 
#' @export

get_ind_scores <- function(data, scale_data = TRUE) {
    data_t <- Matrix::t(data)
    # Perform scaling and centering of the data
    if (scale_data == TRUE) {
        res_pca <- irlba::prcomp_irlba(data_t, scale. = TRUE, center = TRUE)
        } else {
        res_pca <- irlba::prcomp_irlba(data_t)
        }
    individual_scores <- factoextra::get_pca_ind(res_pca)$coord
    individual_scores <- individual_scores[, 1:2]
    rownames(individual_scores) <- colnames(data)
    return(individual_scores)
}

#' Extract individual scores from PCA results for a list of pathways
#' 
#' @param data Numeric matrix with samples in columns, and features in rows
#' @param pathways_list list of pathways
#' @param reg_net Numeric matrix with samples in columns, features in rows
#' @param edges Table, containing information on "reg" and "tar"
#' @param ncores A number of cores to use (default: 1)
#' 
#' @return individual scores for each pathway
#' 
#' @export
#' 
get_pathway_ind_scores <- function(pathways_list,
    reg_net,
    edges,
    scale_data = TRUE) {
    res <- lapply(pathways_list, function(pathway) {
    idx <- which(edges$tar %in% pathway)
    subnet <- reg_net[idx, ]
    individual_scores <- get_ind_scores(subnet, scale_data = scale_data)
    })
    names(res) <- names(pathways_list)
    return(res)
    }

#' Plot PCA colored individuals based on the clinical feature
#' 
#' @param individual_scores PC1 and PC2 scores for individuals
#' @param clin_col clinical feature for individuals (e.g. age)
#' 
#' @return PCA plot with individuals colored according to clinical feature
#' 
#' @export
#' 

plot_clinical_association <- function(scores, clin_col) {
    ggplot(scores, aes(scores[,1], scores[,2], col = clin_col)) +
    geom_point(size = 2.5) + xlab("PC1") + ylab("PC2") +  theme_bw()
}

#' Correlation test or Kruskal Wallis test (numeric or categorical feature)
#' 
#' @param individual_scores PC1 and PC2 scores for individuals
#' @param component principal component to look at (1,2)
#' @param clin_col clinical feature for individuals (e.g. age)
#' 
#' @return pvalue
#' @export


clinical_correlation <- function(individual_scores, component, clin_col) {
    if (is.numeric(clin_col)) {
    # Perform correlation test if the column is numeric
    pval <- cor.test(individual_scores[,component], clin_col)$"p.value"
        } else {
    # Perform Kruskal-Wallis test if the column is not numeric
    pval <- kruskal.test(individual_scores[,component] ~ clin_col)$"p.value"
    }
    return(pval)
    }
