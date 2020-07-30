#! /usr/bin/Rscript
library("ggplot2")
library("argparser")
library("mclust")

theme_set(theme_bw())
my_theme <- theme(plot.title=element_text(face="bold", size=rel(1.8), hjust=0.5), text=element_text(), panel.background=element_rect(colour=NA), plot.background=element_rect(colour=NA), axis.title=element_text(face="bold", size=rel(1)), axis.title.y=element_text(angle=90, vjust=2, size=18), axis.title.x=element_text(angle=0, vjust=-2, size=18), axis.text=element_text(size=16), axis.line=element_line(colour="black"), axis.ticks=element_line(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_rect(colour=NA), legend.position="bottom", legend.direction="horizontal", legend.key.size= unit(0.2, "cm"), legend.title=element_blank(), legend.text=element_text(size=16), plot.margin=unit(c(10,5,5,5), "mm"), strip.background=element_rect(colour="#f0f0f0", fill="#f0f0f0"), strip.text=element_text(face="bold", size=12))

# Accept input from the command line
parser <- arg_parser("Generate histogram of sequence CORE scores.")
parser <- add_argument(parser, "--infile", help="file containing distribution of CORE scores")
parser <- add_argument(parser, "--png", help="filename of the output figure", default="sCORE_dist.png")
args <- parse_args(parser)

dat <- read.csv(args$infile, sep="\t", col.names=c("SeqID","CORE"))

# Use a gaussian mixture model to partition the distribution into two groups
dat.gmm <- Mclust(dat$CORE, G=2)
dat <- cbind(dat, Group=dat.gmm$classification)

g2.mean <- dat.gmm$parameters$mean[2]
g2.var <- dat.gmm$parameters$variance$sigmasq[2]
g2.sd <- sqrt(g2.var)

g2.bw <- ceiling(2 * IQR(subset(dat, Group=="2")$CORE) / length(subset(dat, Group=="2")$CORE)^(1/3))

png(args$png, height=600, width=800)
ggplot(data=dat, aes(x=CORE)) + geom_histogram(aes(y=..density.., group=Group, fill=Group), binwidth=g2.bw) + geom_density(fill="grey", alpha=0.3) + xlab("Sequence CORE index") + ylab("Frequency") + geom_vline(xintercept=g2.mean, linetype="solid", color="black", size=1) + geom_vline(xintercept=c(g2.mean - 2 * g2.sd, g2.mean + 2 * g2.sd), linetype="dashed", color="black", size=1) + my_theme + theme(legend.position="none") + scale_x_continuous(breaks=round(c(seq(0, 100, by=10), unname(g2.mean))))
dev.off()

quit(save="no", status=0, runLast=FALSE)
