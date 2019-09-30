#! /usr/bin/Rscript
library("Nonpareil")
library("argparser")

# Accept input from the command line
parser <- arg_parser("Plot a Nonpareil curve from one or more Nonpareil npo files. Use to predict the sequencing effort required to achieve a desired level of average coverage.")
parser <- add_argument(parser, "--star", help="desired level of average coverage", default=95, type="integer")
parser <- add_argument(parser, "--out", help="filename of the output figure", default="np_curve.svg")
parser <- add_argument(parser, "--infiles", help="one or more Nonpareil npo files", nargs=Inf)
args <- parse_args(parser)

sample.names <- sapply(strsplit(basename(args$infiles), '[.]'), '[', 1)

svg(args$out, height=6, width=8)
if (length(args$infiles > 1)) {
  ncurve <- Nonpareil.curve.batch(args$infiles, libname=sample.names, star=args$star)
  write("Sequencing effort required to reach an average coverage of 95%:", stdout())
  write(ncurve[, 'LRstar'], stdout())
} else {
  ncurve <- Nonpareil.curve(args$infiles, libname=sample.names, star=args$star)
  write(paste("Sequencing effort required to reach an average coverage of 95%:", ncurve$LRstar), stdout())
}
Nonpareil.legend('bottomright')
dev.off()

quit(save="no", status=0, runLast=FALSE)
