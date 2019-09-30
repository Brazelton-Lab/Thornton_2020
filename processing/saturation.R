#! /usr/bin/Rscript
library("ggplot2")
library("argparser")

# Accept input from the command line
parser <- arg_parser("Plot a saturation curve from Khmer's saturate-by-median report output. Tells how much new information is likely to be provided from additional sequencing effort.")
parser <- add_argument(parser, "--threshold", help="coverage threshold provided to saturate-by-median", default=5, type="integer")
parser <- add_argument(parser, "--out", help="filename of the output figure", default="read_saturation.svg")
parser <- add_argument(parser, "--infiles", help="one or more report files from saturate-by-median", nargs=Inf)
args <- parse_args(parser)

accum.curve <- data.frame()
read.totals <- c()
for (infile in args$infiles) {
  sample.name <- sapply(strsplit(basename(infile), '[.]'), '[', 1)
  sample.dat <- read.table(infile, sep=" ", header=FALSE)
  colnames(sample.dat) <- c("x", "y", "z")
  read.totals <- c(read.totals, max(sample.dat$x))
  sample.dat$sample <- sample.name
  accum.curve <- rbind(accum.curve, sample.dat)
}

svg(args$out, height=7, width=7)
ggplot(data=accum.curve, aes(x=x, y=y, color=sample)) + geom_line() + xlab("Number of reads") + ylab(paste0("Reads with coverage less than threshold (=", args$threshold, ")")) + geom_vline(xintercept=read.totals, linetype="dashed", alpha=0.6)
dev.off()

quit(save="no", status=0, runLast=FALSE)
