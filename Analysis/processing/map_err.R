#! /usr/bin/Rscript
library(ggplot2)
library(argparser)

# Accept input from the command line
parser <- arg_parser("Generate a plot of per-base position mapping error rates.")
parser <- add_argument(parser, "--out", help="filename of the output figure", default="mapping_errors.svg")
parser <- add_argument(parser, "--in", help="histogram of mapping error rates from BBmap's mhist argument.")
args <- parse_args(parser)

mhist <- read.table(args$in, sep="\t")
colnames(mhist) <- c("Base", "MatchForward", "SubForward", "DelForward", "InsForward", "NForward", "OtherForward", "MatchReverse", "SubReverse", "DelReverse","InsReverse", "NReverse", "OtherReverse")
mhist <- data.frame(Base=mhist[,"Base"], stack(mhist, select=c("SubForward", "DelForward", "InsForward", "SubReverse", "DelReverse","InsReverse")))

svg(args$out, height=6, width=8)
ggplot(mhist, aes(x=Base, y=values, color=ind)) + geom_line() + xlab("Base Position") + ylab("Mapping Error Rate") + theme(legend.title=element_blank()) + theme_bw()
dev.off()

quit(save="no", status=0, runLast=FALSE)
