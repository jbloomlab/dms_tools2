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
         col_scheme=make_col_scheme(chars=letters,
           cols=letter_colors)
         ) +
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
    p <- p + scale_y_continuous(yname) +
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
