# Extract ggplot legend out of plot --------------------------------------------

my.ggplot.legend <- function(a.gplot) {
  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
  
}


# Arrange plots in regular grid ------------------------------------------------

my.arrange_panels <- function (plotlist,
                               add_tag             = TRUE,
                               remove_panel_legend = TRUE,
                               common_legend       = "none",
                               common_legend_panel = 1,
                               standardise_width   = FALSE,
                               offset_lettering    = 0) {
  
  # get legend
  if (common_legend != "none")
    p.L <- my.ggplot.legend(plotlist[[common_legend_panel]] + theme(legend.position = common_legend))
  
  # add tag and remove legends
  ii <- seq_along(plotlist)
  
  plotlist <-
    lapply(ii, function(i) {
      p <- plotlist[[i]]
      
      if (add_tag)
        p <- p + labs(tag = LETTERS[i + offset_lettering])
      
      if (remove_panel_legend)
        p <- p + theme(legend.position = "none")
      
      return(p)
    })
  
  # convert to gtables
  grobs <- lapply(plotlist, ggplotGrob)
  
  # match the width of plots in a column
  if (standardise_width) {
    
    # get width
    maxWidth <- lapply(grobs, function(g) g$widths[2:5])
    maxWidth <- do.call(grid::unit.pmax, maxWidth)
    
    # apply width
    for(i in 1:length(grobs))
      grobs[[i]]$widths[2:5] <- maxWidth
    
  }
  
  # add legend
  if (common_legend != "none")
    grobs[[length(grobs) + 1]] <- p.L
  
  return(grobs)
  
}


# Simplified ggsave function ---------------------------------------------------

my.ggsave <- function (plot, filename, path = NULL, width, height,
                       units = "in", device = "tiff", dpi = 600,
                       type = "cairo", compression = "lzw") {
  
  # create output folder
  if (!dir.exists(path))
    dir.create(path, recursive = TRUE)
  
  # save plot
  for (i in seq_along(device)) {
    
    out_name <- paste(filename, device[i], sep = ".")
    
    if (device[i] == "tiff") {
      ggsave(plot = plot,
             path = path, filename = out_name,
             width = width, height = height, units = units,
             device = device[i], dpi = dpi,
             type = type, compression = compression)
    } else {
      ggsave(plot = plot,
             path = path, filename = out_name,
             width = width, height = height, units = units,
             device = device[i], dpi = dpi)
    }
  }
  
}
