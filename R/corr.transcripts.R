#' @title Calculate bootrstapped correlation coefficents for genes from SEURAT object
#'
#' @details This function calculates bootstrapped correlation coefficients for a target gene against all other genes.
#' @param object Seurat object
#' @param target Gene for which correlation analyses against all other genes should be calculated
#' @param ident.use Only include particular clusters, as defined in Seurat object
#' @param thresh Only consider genes that are expressed in a percentage of cells greater than thresh
#' @param method Method to calculate correlation (default is "pearson")
#' @param nbootstrap Number of bootstrap iterations to perform (default is 30)
#' @param nThreads Number of threads to perform computation on (default is 10)
#' @param slot Slot to pull data from (default is "data")
#' @importFrom Seurat subset LayerData
#' @return Data frame containing bootstrapped correlation coefficients
#' @export

corr.transcripts<-function(object,target,ident.use=NULL,thresh=0,method="pearson",nbootstrap=30,nThreads=10, slot = "data"){
  stopifnot(target %in% rownames(object))
  stopifnot(is.null(ident.use) || all(ident.use %in% levels(object)))
  message(paste("Clusters found: ", ident.use[ident.use %in% levels(object)]))

  message(paste("Dimensions of imported matrix:", dim(object)[1],"x",dim(object)[2]))
  if (!is.null(ident.use)){
    object<-subset(object,ident.use=ident.use)
    message(paste("Dimensions of imported matrix after filtering for identitiy classes:", dim(object)[1],"x",dim(object)[2]))
  }
  mtx<-as.matrix(LayerData(object, slot = slot))
  valid_methods <- c("pearson", "spearman", "kendall")  # Add other valid methods
  method <- match.arg(method, valid_methods)

  out<-NULL
  message("Bootstrapping correlation coefficients")
  out<-cor.bootstrap(mtx,target=target,nbootstrap=nbootstrap,method=method,nThreads=nThreads)
  return(out)
}
