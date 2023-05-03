#' Run PCA analysis on the data
#'
#' This function perform PCA analysis on the data
#'
#' @import irlba
#' 
#' @param data Numeric matrix with samples in columns, and features in rows
#' @param scale_data Logical, whether to scale the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' @param center_data Logical, whether to center the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' @return Dataframe with pca results

run_pca <- function(data, scale_data = TRUE, center_data = TRUE) {
    data_t <- t(data)
    # Perform scaling and centering of the data
    res_pca <- irlba::prcomp_irlba(data_t, scale. = scale_data,
                center = center_data)
    # Extract variance explained by the first PC
    pc1 <- summary(res_pca)$importance[2, 1] * 100
    pca_result <- data.frame("pc1" = pc1, "n_edges" = ncol(data_t))
    return(pca_result)
}
