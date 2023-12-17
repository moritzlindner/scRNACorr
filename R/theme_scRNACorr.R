#' Theme optimization for scRNACorr
#' @importFrom ggpubr theme_pubr
#' @importFrom ggplot2 theme element_text unit
#' @import ggplot2
#' @return none

theme_scRNACorr<-function(font.size=8){
  out<-theme_pubr()
  out<-out+theme(axis.text.x = element_text(angle = 0, vjust = 0.5,
                                                              size = font.size, hjust = 1),
                          axis.text.y = element_text(angle = 0, vjust = 0.5,
                                                              size = font.size, hjust = 1),
                          plot.margin = unit(c(0,0,0,0), "cm"),
                          axis.title =  element_text(size = font.size,face= "bold"),
                          legend.text = element_text(angle = 0, vjust = 0.5,
                                                              size = font.size, hjust = 1),
                          legend.title = element_text(angle = 0, vjust = 0.5,
                                                               size = font.size, hjust = 1),
                          legend.position="bottom",
                          legend.justification = "center",
                          legend.title.align = 0.5
  )
  return(out)
}
