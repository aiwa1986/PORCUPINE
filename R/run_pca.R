#' Run PCA analysis on the data
#'
#' This function perform PCA analysis on the data
#' 
#' @import irlba
#' 
#' @param data Numeric matrix with samples in columns, and features in rows
#' 
#' @return Dataframe with pca results
#' 
#' @export

run_pca <- function(data, scale_data = TRUE) {
    data_t <- Matrix::t(data)
    # Perform scaling and centering of the data
    if (scale_data == TRUE) {
        res_pca <- irlba::prcomp_irlba(data_t, scale. = TRUE, center = TRUE)
        } else {
        res_pca <- irlba::prcomp_irlba(data_t)
        }
    # Extract variance explained by the first PC
    pc1 <- summary(res_pca)$importance[2, 1] * 100
    pca_result <- data.frame("pc1" = pc1, "n_edges" = ncol(data_t))
    return(pca_result)
    }
