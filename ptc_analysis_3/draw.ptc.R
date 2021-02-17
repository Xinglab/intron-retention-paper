theme_Publication <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=margin(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

theme_empty <- function(base_size=14) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size)
    + theme(rect = element_rect(fill = "transparent"),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_text(color = 'transparent'),
            axis.ticks = element_blank(),
            panel.border=element_blank(),
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            axis.line = element_line(colour="transparent"),
            axis.title = element_text(colour="transparent")
    ))
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#ef3b2c","#7fc97f","#fdb462","#fb9a99","#a6cee3","#984ea3","#662506", "#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#ef3b2c","#7fc97f", "#fdb462","#fb9a99","#a6cee3","#984ea3", "#662506", "#ffff33")), ...)
  
}

scale_fill_compartment <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#fdb462","#fb9a99","#a6cee3","#984ea3","#662506", "#ffff33")), ...)
  
}

library('ggplot2')
library(ggpubr)
data <- read.table('ptc_response_xmeans.txt', header=T, sep='\t')
data$cell <- factor(data$cell, levels= c('ESC', 'NPC', 'Ctx'))
p <- ggplot(data, aes(x = category, y = proportion, fill  = type)) + geom_bar(stat = "identity") + facet_wrap(~cell) + theme_Publication() + scale_fill_Publication() + 
  theme(axis.title.x = element_blank(),strip.background = element_blank(), panel.grid.major.x = element_blank(), plot.title = element_text(hjust = 0.5), legend.title = element_blank(), strip.text = element_text(size=15), text = element_text(size=15), panel.spacing = unit(2, 'lines')) +
  ylab('Proportion') + scale_y_continuous(expand = c(0, 0))
print(p)
ggsave(p, file="PTC_figure.png", width=6, height=3, dpi = 300)
