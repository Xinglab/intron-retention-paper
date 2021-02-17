library(ggplot2)
library(Rtsne)
data <- read.table('tsne.input_A_D.txt', header=TRUE, sep = '\t')
## obtain data matrix
data.matrix <- as.matrix(data[,2:(dim(data)[2] - 2)])
## perform T-SNE analysis
set.seed(1890)
tsne <- Rtsne(data.matrix, check_duplicates = FALSE, pca = TRUE, 
              perplexity=5, theta=0.5, dims=2)
embedding <- as.data.frame(tsne$Y)
embedding$Class <- data[,1]
## get embedded values
embedding$strength <- data[,dim(data)[2]-1]
embedding$Category <- data[, dim(data)[2]]
p <- ggplot(embedding, aes(x=V1, y=V2, color=Category, label = Class)) + geom_point(aes(size=strength, alpha = 0.8)) + geom_text( hjust=0.3, vjust=-0.8, size = 4.5) +
  scale_alpha(guide = 'none') + scale_size_continuous(guide=FALSE) +
  theme_bw() + theme(panel.border = element_rect(colour='black',fill=NA),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     text = element_text(size=15),
                     legend.title = element_blank(),
                     legend.text=element_text(size=13),
                     strip.text = element_text(size=13),
                     legend.key = element_rect(fill = "transparent", colour = "transparent"),
                     legend.background = element_rect(fill = 'transparent', colour = 'transparent')) + 
  scale_x_continuous(expand = c(0, 0),limits = c(-150,230)) + 
  scale_color_manual(values=c("#F8766D", "#619CFF", "#56B4E9")) + 
  xlab('t-SNE 1') + ylab('t-SNE 2')
ggsave(p, file="tsne_plot.pdf", width=8, height=5.2)