#' Logo plot for subset of sites
#'
#' @param mat matrix of data to plot
#' @param plotfile name of created plot file
#' @param width width of plot
#' @param height height of plot
#' @param xlabels vector of x-tick labels
#' @param vertlines vector of vertical dotted line locations
#' @param yname y-axis label name
#' @param chars vector of characters in logo
#' @param char_colors vector of colors for each character
#' @param ylimits vector giving y minimum and maximum
#' @param title string giving title to put above plot
siteSubsetGGSeqLogo <- function(mat,
                                plotfile,
                                width,
                                height,
                                xlabels,
                                vertlines,
                                yname,
                                chars,
                                char_colors,
                                ylimits=NULL,
                                title='') {
  p <- ggseqlogo(mat, method='custom', ncol=length(xlabels),
         col_scheme=make_col_scheme(chars=chars, cols=char_colors)) +
       scale_x_continuous('', breaks=1:length(xlabels),
         labels=xlabels)

  if (nchar(trimws(yname))) {
    p <- p + scale_y_continuous(yname, limits=ylimits)
    axis.line.y = element_line(color='black')
    axis.text.y = element_text()
  } else {
    p <- p + scale_y_continuous(limits=ylimits)
    axis.text.y <- element_blank()
    axis.line.y <- element_blank()
  }

  if (length(vertlines) > 0) {
    p <- p + geom_vline(xintercept=vertlines, linetype='dotted',
                        color='black')
  }

  p <- p + theme(axis.text.x=element_text(angle=90, hjust=1),
         axis.line.x=element_line(color='black'),
         axis.text.y=axis.text.y,
         axis.line.y=axis.line.y,
         axis.text=element_text(size=12),
         axis.title=element_text(size=13),
         plot.title=element_text(hjust=0.5)
         )

  if (nchar(trimws(title))) {
    p <- p + ggtitle(title)
  }
  
  ggsave(plotfile, plot=p, width=width, height=height)
}


#' Faceted sequence logo
#'
#' @param matrices Named vector of matrices
#' @param plotfile Name of created plot file
#' @param ncol Number of columns in plot
#' @param width Width of plot
#' @param height Height of plot
#' @param xname x-axis label name
#' @param xlabels Vector of labels for x-axis ticks.
#' @param xlabelsrotate Rotate the x-axis tick labels?
#' @param yname y-axis label name
#' @param chars Vector of characters in logo
#' @param char_colors Vector of colors for each character
facetedGGSeqLogo <- function(matrices,
                             plotfile,
                             ncol,
                             width,
                             height,
                             xname,
                             xlabels,
                             xlabelsrotate,
                             xline,
                             yname,
                             chars,
                             char_colors) {
  p <- ggseqlogo(matrices, method='custom', ncol=ncol,
         col_scheme=make_col_scheme(chars=chars, cols=char_colors)) +
       scale_x_continuous(xname, breaks=1:length(xlabels),
         labels=xlabels) 

  if (xlabelsrotate) {
    axis.text.x = element_text(angle=90, hjust=1)
  } else {
    axis.text.x = element_text()
  } 

  if (xline) {
    axis.line.x <- element_line(color='black')
  } else {
    axis.line.x <- element_blank()
  }

  if (nchar(trimws(yname))) {
    p <- p + scale_y_continuous(yname)
    axis.line.y = element_line(color='black')
    axis.text.y = element_text()
  } else {
    axis.text.y <- element_blank()
    axis.line.y <- element_blank()
  }

  p <- p + theme(axis.text.x=axis.text.x,
         axis.line.x=axis.line.x,
         axis.text.y=axis.text.y,
         axis.line.y=axis.line.y,
         axis.text=element_text(size=12),
         strip.text=element_text(size=13),
         axis.title=element_text(size=13),
         panel.spacing=unit(1.75, 'lines')
         )

  ggsave(plotfile, plot=p, width=width, height=height)
}
