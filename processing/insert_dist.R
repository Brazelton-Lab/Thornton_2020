#! /usr/bin/Rscript
library("ggplot2")
library("argparser")

# Accept input from the command line
parser <- arg_parser("generate a histogram of insert sizes from BBmap's ihist output or any tab-separated file containing a list of insert sizes as the first column and their counts as the second")
parser <- add_argument(parser, "--outplot", help="filename of the output figure", default="insert_size.svg")
parser <- add_argument(parser, "--infile", help="file of insert sizes and their counts")
parser <- add_argument(parser, "--bins", help="number of bins to divide the data among", default=30)
parser <- add_argument(parser, "--xmax", help="maximum insert size to display on the histogram", default=1500)
args <- parse_args(parser)

ishist <- read.table(args$infile, sep="\t")

svg(args$outplot, width=6, height=6)
ggplot(data.frame(insert=rep(ishist$V1, times=ishist$V2)), aes(insert)) + geom_histogram(stat="bin", bins=args$bins) + xlim(0, args$xmax) + xlab("Insert Size") + ylab("Frequency")
dev.off()

quit(save="no", status=0, runLast=FALSE)
