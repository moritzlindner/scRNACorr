#' @title Calculate bootrstapped correlation coefficents for genes from matrix
#'
#' @details This function creates an ROC plot. Returns or prints a ggplot object.
#' @param mtx matrix (transcripts x cells) containing scaled expression values
#' @param target Gene for wich correlation analyses against all other genes should be calculated
#' @param method Method to calculate correlation
#' @param nbootstrap Number of bootrstrap iterations to perform
#' @param nThreads number of Threads to perform computation on
#' @importFrom WGCNA cor
#' @return data.frame containing Tukey Box and Wishker Plot statistics (columns) for all transcripts (rows) in the matrix
#' @export

cor.bootstrap<-function(mtx,target,nbootstrap=nbootstrap,method="pearson",nThreads=10, silent = F) {
  mtx<-t(mtx)
  corr<- matrix(NA, nrow = ncol(mtx), ncol = nbootstrap)
  rownames(corr)<-colnames(mtx)
  normdist<-corr

  if (!silent){
    pb <- txtProgressBar(min = 1, max = nbootstrap, style = 3)
  }
  for (i in 1:nbootstrap){
    curr<-mtx[sample(1:nrow(mtx),replace = T),]
    corr[,i]<-cor(curr,curr[,colnames(curr)==target], method=method,nThreads = nThreads)
    normdist[,i]<-pnorm(scale(corr[,i]))
    if (!silent){
      setTxtProgressBar(pb, i)
    }
    gc(full=T)
  }
  close(pb)
  rm(curr)
    if (!silent){
    rm(pb)
  }  
  corr<-t(apply(corr, 1, function(x){boxplot.stats(x,do.conf=FALSE,do.out = FALSE)$stats}))
  colnames(corr)<-c("min","p25","median","p75","max")

  normdist<-t(apply(normdist, 1, function(x){boxplot.stats(x,do.conf=FALSE,do.out = FALSE)$stats}))
  colnames(normdist)<-c("pct.smaller.min","pct.smaller.p25","pct.smaller.median","pct.smaller.p75","pct.smaller.max") #,"pct.smaller.p"
  gc(full=TRUE)
  return(cbind(corr,normdist))
}
