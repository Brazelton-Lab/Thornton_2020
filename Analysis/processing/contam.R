#! /usr/bin/Rscript
library("ggplot2")
library("argparser")

# Accept input from the command line
parser <- arg_parser("")
parser <- add_argument(parser, "--outpref", help="filename prefix of the output figures", default="contam")
parser <- add_argument(parser, "--infile", help="comma-separated file (CSV) with three columns: InputDNA, FractionAdapter, and LibraryConcentration")
args <- parse_args(parser)

dat <- read.table(args$infile, sep=",", header=TRUE)

out1 <-  paste0(args$outpref, ".adapter_contam.svg")
svg(out1, height=6, width=8)
ggplot(dat, aes(x=InputDNA, y=FractionAdapter)) + geom_point(size=4) + geom_smooth(method="loess") + xlab("Amount DNA input into library preparation") + ylab("Percent 5' adapter contaminants")
dev.off()

out2 <-  paste0(args$outpref, ".library_concentration.svg")
svg(out2, height=6, width=8)
ggplot(dat, aes(x=InputDNA, y=LibraryConcentration)) + geom_point(size=4) + geom_smooth(method="loess") + xlab("Amount DNA input into library preparation") + ylab("Final library DNA concentration")
dev.off()

quit(save="no", status=0, runLast=FALSE)
