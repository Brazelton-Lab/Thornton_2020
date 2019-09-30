#! /usr/bin/Rscript
library("ggplot2")
library("argparser")

# Accept input from the command line
parser <- arg_parser("Plot a k-mer abundance distribution from Khmer's abundance-dist.py output.")
parser <- add_argument(parser, "--xmax", help="maximum k-mer abundance to display on plot [default: first kmer abundance with frequency zero]", type="integer")
parser <- add_argument(parser, "--out", help="filename of the output figure", default="kmer_spectrum.svg")
parser <- add_argument(parser, "--in", help="k-mer abundance distribution file from abundance-dist.py. kmer abundances should be given in the first column and their respective counts in the second. A column header is required.")
args <- parse_args(parser)

sample.name <- sapply(strsplit(basename(args$in), '[.]'), '[', 1)
kmers <- read.table(args$in, sep=",", header=TRUE)

xmax <- ifelse(is.na(args$xmax), match(0, kmers$count[2:length(kmers$count)]), args$xmax)
remaining <- sum(kmers$count[xmax: length(kmers$count)]) / sum(kmers$count)  #fraction kmers with abundance greater than xmax 
write(paste("Percentage of kmers left unplotted:", remaining * 100), stdout())

svg(args$out, height=6, width=8)
ggplot(data=kmers[2:xmax,], aes(x=abundance, y=count)) + geom_density(stat="identity", fill="#4271AE", color="#1F3552", alpha=0.8) + xlab("K-mer abundance") + ylab("log10 frequency") + scale_y_log10() + theme_bw()
dev.off()

quit(save="no", status=0, runLast=FALSE)
