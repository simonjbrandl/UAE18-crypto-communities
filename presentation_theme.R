theme_presentation <- function (base_size = 11, base_family = "") {
  theme_bw() %+replace% 
    theme(panel.grid.major  = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "white"),
          axis.ticks = element_line(color = "white"),
          axis.text = element_text(color = "white"),
          axis.title = element_text(color = "white"))
}


theme_sjb <- function (base_size = 8, base_family = "") {
  theme_bw() %+replace% 
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 8),
          axis.title = element_text(color = "black", size = 8),
          legend.title=element_text(size=8), 
          legend.text=element_text(size=8))
}
