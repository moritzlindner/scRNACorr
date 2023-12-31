#' Calculate and visualize bootstrapped correlation coefficients for genes from a SEURAT object
#'
#' This function generates a histogram of bootstrapped correlation coefficients for genes from a correlation table, along with an optional box and whisker plot for a Gene of Interest (GOI).
#'
#' @param corrtbl A correlation table as generated by the corr.transcripts function.
#' @param GOI Gene of Interest (GOI) to be highlighted in the plot. Default is NULL.
#' @param bins Number of bins to split data in for the histogram.
#' @param alpha Transparency of the histogram bars (0 for fully transparent, 1 for opaque).
#' @param probs Percentile borders to plot on the histogram (default is c(0.75, 0.95, 0.99)).
#' @param font.size Font size for axis labels and text annotations.
#' @return A ggplot2 object containing the correlation coefficient histogram and, if specified, a box and whisker plot for the GOI.
#' @import ggplot2
#' @importFrom ggplot2 ggplot aes geom_vline geom_hline geom_histogram geom_boxplot geom_text theme element_text element_blank xlim ylab
#' @importFrom cowplot plot_grid
#' @seealso
#' \code{\link{corr.transcripts}}
#' @export

cor.plot <-
  function(corrtbl,
           GOI = NULL,
           bins = 400,
           alpha = 1,
           probs = c(.75, .95, .99),
           font.size = 8) {
    stopifnot(is.matrix(corrtbl))

    if (!is.null(GOI)) {
      stopifnot(GOI %in% rownames(corrtbl))
    }

    corrtbl <- as.data.frame(corrtbl)
    corrtbl$gn <- rownames(corrtbl)
    low <- min(corrtbl$median, na.rm = T)
    up <- max(corrtbl$median[corrtbl$median < 1], na.rm = T)

    if (!is.null(probs)) {
      warning(
        "quantiles calculated from medians after bootstrapping. these thus differ from those provided in the correlation table obteined from corr.transcripts."
      )
      quantile.borders <- NULL
      quantile.borders$x <-
        quantile(corrtbl$median, na.rm = T, probs = probs)
      quantile.borders$name <-
        paste0(as.character(probs * 100), "th percentile")
      quantile.borders <- as.data.frame(quantile.borders)
    }



    top <- max(table(cut(corrtbl$median, bins)), na.rm = T)

    hist <- ggplot(corrtbl, aes(x = median)) +
      geom_vline(
        data = quantile.borders,
        inherit.aes = F,
        aes(xintercept = x),
        linetype = 2,
        colour = "grey"
      ) +
      geom_histogram(
        position = "dodge",
        aes(y = ..count..),
        alpha = alpha,
        bins = bins
      ) +
      geom_text(
        data = quantile.borders,
        inherit.aes = F,
        aes(x = x, label = name),
        y = top * 0.5,
        angle = 90,
        text = element_text(size = font.size)
      ) +
      theme_scRNACorr(font.size = font.size) +
      theme(legend.position = "none") +
      xlim(c(low, up)) +
      ylab("Count")


    if (!is.null(GOI)) {
      hist <- hist +
        theme(
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()
        )

      bp <- ggplot(
        corrtbl[rownames(corrtbl) %in% GOI, ],
        aes(
          x = gn,
          middle = median,
          ymin = min,
          ymax = max,
          lower = p25,
          upper = p75
        )
      ) +
        geom_hline(
          data = quantile.borders,
          inherit.aes = F,
          aes(yintercept = x),
          linetype = 2,
          colour = "grey"
        ) +
        geom_boxplot(stat = "identity") +
        ylab("Correlation coefficient") +
        theme(legend.position = "none") +
        theme_scRNACorr(font.size = font.size) +
        ylim(low, up) +
        coord_flip() +
        theme(axis.title.y = element_blank(),
              axis.ticks.y = element_blank())

      strecher <- sum(rownames(corrtbl) %in% GOI) / 6

      out <- plot_grid(
        hist,
        bp,
        ncol = 1,
        align = "v",
        axis = "lr",
        rel_heights = c(1, strecher)
      )
    } else{
      out <- hist
    }


    return(out)
  }
