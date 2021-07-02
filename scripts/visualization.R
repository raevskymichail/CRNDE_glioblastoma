#!/usr/bin/env Rscript

#' Set a publication-ready theme for plots
#'
#' @import grid
#' @import ggthemes
#' @export
theme_Publication <- function(base_size = 14, base_family = "helvetica") {
  (theme_foundation(base_size = base_size, base_family = base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold", size = rel(1)),
               axis.title.y = element_text(angle = 90, vjust = 2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(),
               axis.line = element_line(colour = "black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour = "#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size = unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face = "italic"),
               plot.margin = unit(c(10, 5, 5, 5), "mm"),
               strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
               strip.text = element_text(face = "bold")
          ))

}

#' @import scales
#' @export
scale_fill_Publication <- function(...) {
  discrete_scale("fill", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ...)
}

#' @import scales
#' @export
scale_colour_Publication <- function(...) {
  discrete_scale("colour", "Publication", manual_pal(values = c("#386cb0", "#fdb462", "#7fc97f", "#ef3b2c", "#662506", "#a6cee3", "#fb9a99", "#984ea3", "#ffff33")), ...)
}

#' @import scales
#' @export
scale_four_colour_Publication <- function(...) {
  discrete_scale("colour", "Publication", manual_pal(values = c("#fb746c", "#eac94c", "#04bb3b", "#619cff")))
}

#' Save plot of provided analysis step
#'
#' Export ggplot2 object to desired location in vector (`.pdf`) and raster (`.tiff`) formats
#' @import ggplot2
#' @import ggpubr
#' @import extrafont
#' @export
export_analysis_plot <- function(filename, plot, path,
                                 scale, width, height, units, cairo_pdf_device = TRUE, load_fonts = FALSE) {
  if (load_fonts == TRUE) {
    extrafont::font_import(prompt = FALSE)
  }

  extrafont::loadfonts(device = "pdf")

  if (cairo_pdf_device == TRUE) {
    ggsave(
      filename = paste0(filename, ".pdf"),
      plot = plot,
      device = cairo_pdf,
      path = path,
      scale = scale,
      width = width,
      height = height,
      units = units,
      dpi = 300)
  } else {
    ggsave(
      filename = paste0(filename, ".pdf"),
      plot = plot,
      path = path,
      device = "pdf",
      scale = scale,
      width = width,
      height = height,
      units = units,
      dpi = 300)
  }

  ggsave(
    filename = paste0(filename, ".tiff"),
    plot = plot,
    device = "tiff",
    path = path,
    scale = scale,
    width = width,
    height = height,
    units = units,
    compression = "lzw",
    dpi = 1200)
}

#' PCA plot
#'
#' Creates a PCA plot (ggplot2 object) with specified PC variance
#'
#' @param data_for_pca
#'
#' @import FactoMineR
#' @import factoextra
#' @import ggplot2
#' @import ggpubr
#'
#' @export
plot_pca <- function(data_for_pca) {
  data.pca <- PCA(data_for_pca, graph = FALSE)
  pcs_variance <- as.data.frame(data.pca$eig)$`percentage of variance`
  pc1_lab <- paste0("PC1 (", round(pcs_variance[1], 2), "%)")
  pc2_lab <- paste0("PC2 (", round(pcs_variance[2], 2), "%)")

  p <- fviz_pca_ind(data.pca,
                  geom.ind = "point", # show points only (nbut not "text")
                  col.ind = metadata$dataset, # color by groups
                  addEllipses = FALSE, # Concentration ellipses
                  legend.title = NULL
                  ) +
                  labs(x = pc1_lab, y = pc2_lab) +
                  theme_Publication() +
                  theme(legend.title = element_blank(),
                        plot.title = element_blank())
  return(p)
}