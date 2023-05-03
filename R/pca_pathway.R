#' Run PCA analysis for a list of pathways performed on network edges
#'
#' This function runs PCA analysis for a list of pathways
#' 
#' @import parallel
#' 
#' @param pathways_list list of pathways 
#' @param reg_net Numeric matrix with samples in columns, features in rows
#' @param edges Table, containing information on "reg" and "tar"
#' @param ncores A number of cores to use (default: 1)
#' @param scale_data Logical, whether to scale the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' @param center_data Logical, whether to center the data (TRUE) or not (FALSE),
#' Default is TRUE.
#' 
#' @return Dataframe with pca results for pathways in a pathway file
#' @export

pca_pathway <- function(pathways_list,
                reg_net,
                edges,
                ncores = 1,
                scale_data = TRUE,
                center_data = TRUE) {
    res <- parallel::mclapply(pathways_list, function(pathway) {
                idx <- which(edges$tar %in% pathway)
                subnet <- reg_net[idx, ]
                pca_result <- run_pca(subnet,
                              scale_data = scale_data,
                              center_data = center_data)
    }, mc.cores = ncores)
    res <- as.data.frame(do.call("rbind", res))
    res$pathway <- rownames(res)
    res$pathway_size <- lengths(pathways_list)
    res <- res[, c("pathway", "pc1", "n_edges", "pathway_size")]
    return(res)
}
