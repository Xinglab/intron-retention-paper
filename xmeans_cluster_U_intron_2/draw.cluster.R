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
data.combine <- data.frame()
data <- read.table('ESC.cluster.polyA.txt', header =T, sep = '\t')
data.combine <- rbind(data.combine, data.frame(Cell_Compartment = data$Cell_Compartment, FI = data$FI, Cluster = data$Cluster, cell = 'mESC'))
data <- read.table('NPC.cluster.polyA.txt', header =T, sep = '\t')
data.combine <- rbind(data.combine, data.frame(Cell_Compartment = data$Cell_Compartment, FI = data$FI, Cluster = data$Cluster, cell = 'mNPC'))
data <- read.table('Ctx.cluster.polyA.txt', header =T, sep = '\t')
data.combine <- rbind(data.combine, data.frame(Cell_Compartment = data$Cell_Compartment, FI = data$FI, Cluster = data$Cluster, cell = 'mCtx'))
data.combine$Cell_Compartment <- factor(data.combine$Cell_Compartment, levels = c('Chr', 'Nuc', 'Cyto'))
data.combine$cell <- factor(data.combine$cell, levels = c('mESC', 'mNPC', 'mCtx'))
p <- ggplot(data.combine, aes(x = Cluster, y = FI, fill = Cell_Compartment)) + geom_boxplot(outlier.alpha = 0.001) + facet_wrap(~cell, scales = "free_x") + theme_Publication() + scale_fill_Publication() +
  theme(legend.position = "right", legend.title = element_blank(), strip.text = element_text(size=15), legend.direction = "vertical") + scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  xlab('Cluster') + ylab('PI value')

ggsave('cluster.png', width=10, height=5, dpi = 300)