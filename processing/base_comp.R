#! /usr/bin/Rscript
library(ggplot2)
library(argparser)

# Accept input from the command line
parser <- arg_parser("Plot the distribution of nucleotide frequencies by position in the reads.")
parser <- add_argument(parser, "--out", help="filename of the output figure", default="base_frequencies.svg")
parser <- add_argument(parser, "--in", help="histogram of nucleotide frequencies from BBduk's bhist argument.")
parser <- add_argument(parser, "--maxl", help="maximum length of a read in the dataset", default=125, type="integer")
args <- parse_args(parser)

bhist <- data.frame(read.table(args$in, sep="\t", row.names=1), strand=rep(c("Forward", "Reverse"), times=c(args$maxl, args$maxl)), base=rep(0:(args$maxl-1), times=2))
colnames(bhist) <- c("A", "C", "G", "T", "N", "strand", "base")
bhist <- data.frame(bhist[,c("base", "strand")], stack(bhist, select=c("A", "C", "G", "T", "N")))

svg(args$out, height=6, width=9)
ggplot(bhist, aes(x=base, y=values, color=ind)) + geom_line() + facet_grid(~strand) + xlab("Base Position") + ylab("Frequency") + theme(legend.title=element_blank()) + theme_bw()
dev.off()

quit(save="no", status=0, runLast=FALSE)
