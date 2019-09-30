#! /usr/bin/Rscript
library("ggdendro")
library("ggplot2")
library("vegan")
library("argparser")

# Accept input from the command line
parser <- arg_parser("Plot a dendrogram from hierarchical clustering of minhash signature similarity values.")
parser <- add_argument(parser, "--sep", help="character that separates sample names from the other components of the filenames in the header.", default=".")
parser <- add_argument(parser, "--pos", help="position of the sample name in the filenames of header.", default=1, type="integer")
parser <- add_argument(parser, "--out", help="filename of the output figure", default="minhash_dist.dendro.svg'")
parser <- add_argument(parser, "--in", help="matrix of Jaccard similarities in CSV format. The header should consist of the filenames the similarity matrix was generated from.", nargs=Inf)
args <- parse_args(parser)

hashdist <- read.table(args$in, sep=",", header=TRUE, check.names=FALSE)
samples <- colnames(hashdist)
sample.sep <- paste0("[", args$sep, "]") 
samples <- sapply(strsplit(samples, sample.sep), '[', args$pos)
colnames(hashdist) <- samples
hashdist <- as.dist(hashdist)
hashdist <- 1 - hashdist
hashdist.clust <- hclust(hashdist, method='average')
hashdist.dendro <- dendro_data(as.dendrogram(hashdist.clust), type='rectangle')

svg(args$out, width=6, height=6)
ggplot(data=hashdist.dendro$segment) + geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) + geom_text(data=hashdist.dendro$label, aes(x=x, y=y, label=label, hjust=1.1), angle=90, lineheight=10) + scale_y_continuous(expand=c(0, 0.2)) + ylab("1 - Jaccard Similarity") + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()

quit(save="no", status=0, runLast=FALSE)
