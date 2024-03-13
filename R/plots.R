

#' @title Plot 2 histograms
#' @description Plot 2 histograms on the same graph.
#' xxxxx
#'
#' @param d1 vector of values for the first histogram
#' @param d2 vector of values for the first histogram
#' @param name1 Label for first histogram
#' @param name2 Label for 2nd histogram
#' @param titlename Title of figure
#' @param xlab X-axis label
#' @param freq If True, bins heights correspond to raw counts, otherwise bins 
#' are normalized.
#'
#' @import grDevices
#' 
#' @examples
#' NULL
#' 
#' @return A plot
#' 
#'
plot2hists <- function(d1,
  d2,
  name1,
  name2,
  titlename,
  xlab = "",
  freq = TRUE) {
  c1 <- rgb(173,216,230, maxColorValue = 255, alpha = 100, names = "lt.blue")
  c2 <- rgb(255,192,203, maxColorValue = 255, alpha = 100, names = "lt.pink")
  b <- min(c(d1,d2), na.rm = TRUE) - 0.01
  e <- max(c(d1,d2), na.rm = TRUE)*(1 + 0.1)
  ax <- pretty(b:e, n = 50)
  if (freq) {
    hgA <- hist(d1, breaks = ax, plot = FALSE) # Save first histogram data
    hgB <- hist(d2, breaks = ax, plot = FALSE) # Save 2nd histogram data
    plot(hgA, col = c1, main = titlename, freq = freq,
      ylim = c(0, max(c(hgA$counts, hgB$counts))))
    plot(hgB, col = c2, add = TRUE, freq = freq)
  } else {
    hgA <- hist(d1, breaks = ax, plot = FALSE) # Save first histogram data
    hgA$density = hgA$counts/sum(hgA$counts)*100
    hgB <- hist(d2, breaks = ax, plot = FALSE) # Save first histogram data
    hgB$density = hgB$counts/sum(hgB$counts)*100
    plot(hgA, col = c1, main = titlename, xlab = xlab, freq = freq,
      ylim = c(0, max(c(hgA$density, hgB$density))))
    plot(hgB, col = c2, add = TRUE, freq = freq)
  }
  legend("topleft", c(name1, name2), fill = c(c1, c2))
}

#' @title Empirical density of peptide correlations
#' @description Plot empirical densities of correlations between peptides within
#'  PG and at random, estimated by gaussian kernel. Note that only correlations 
#'  between fully observed peptides are considered here.
#'
#' @param pep.data List representing dataset
#' @param titlename Title of the graph displayed
#' @param xlabel Label of x-axis
#' 
#' @import ggplot2
#'
#' @return The ggplot2 graph
#' @export
#'
#' @examples
#' data(subbouyssie)
#' plot_pep_correlations(subbouyssie, 'test')
#'
plot_pep_correlations <- function(pep.data, 
  titlename = NULL, 
  xlabel = "Correlations") {
  allcors = list()
  for (i in 1:ncol(pep.data$adj)) {
    pep.idx = which(pep.data$adj[, i, drop = FALSE] == 1)
    if (length(pep.idx) != 1) {
      pep_abs_pg = pep.data$peptides_ab[, pep.idx]
      cor_pg = cor(pep_abs_pg)
      mask_tri_low = lower.tri(cor_pg)
      allcors[[i]] = cor_pg[mask_tri_low]
    }
  }
  all_cors_PG_vec = unlist(allcors)
  allcors = list()
  for (i in 1:ncol(pep.data$adj)) {
    pep.idx = sample(nrow(pep.data$adj), sum(pep.data$adj[, i, drop = FALSE]))
    if (length(pep.idx) != 1) {
      pep_abs_pg = pep.data$peptides_ab[, pep.idx]
      cor_pg = cor(pep_abs_pg)
      mask_tri_low = lower.tri(cor_pg)
      allcors[[i]] = cor_pg[mask_tri_low]
    }
  }
  all_cors_rand_vec = unlist(allcors)
  data.hist = data.frame(values = c(all_cors_PG_vec, all_cors_rand_vec),
    group = factor(
      c(rep("Within PG", length(all_cors_PG_vec)),
        rep("Random", length(all_cors_rand_vec)))))
  g <- ggplot2::ggplot(data.hist, 
    ggplot2::aes(x = values, fill = group)) + xlab(xlabel) +
    # geom_histogram(position = "identity", alpha = 0.2) +
    geom_density(alpha=.2, na.rm = TRUE) + xlim(c(-1, 1)) +
    theme(legend.title=element_blank(),
      # legend.position = c(0.8, 0.9),
      panel.background = element_blank(),
      aspect.ratio = 1,
      legend.key = element_rect(fill = "white"),
      plot.title = element_text(hjust = 0.5))
  if (!is.null(titlename)) {
    g <- g + ggtitle(titlename)
  }
  g
}
