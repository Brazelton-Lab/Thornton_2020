#! /usr/bin/Rscript
library("ggplot2")
library("ggrepel")
library("argparser")
library("reshape2")

# Accept input from the command line
parser <- arg_parser("Determine reasonable threshold values for a custom profile HMM")
parser <- add_argument(parser, "in1", help="table of sequence and domain bitscores from custom HMM")
parser <- add_argument(parser, "in2", help="table of sequence and domain bitscores from Pfam HMMs")
parser <- add_argument(parser, "--png", help="output figure of bitscore changes between matches to custom HMM and matches to Pfam", default="bitscores.png")
parser <- add_argument(parser, "--out", help="output file of sequences suggested for incorporation into final alignment")
args <- parse_args(parser)

# Functions
find_thresholds <- function(x, fp) {
  x.sorted <- sort(x)
  fp.pos <- which(x.sorted == fp)[1]
  ga <- NA
  for(i in fp.pos:length(x.sorted)) {
    ga.prov <- x.sorted[i]
    if(!is.na(ga.prov) & ga.prov > fp) {
      ga <- ga.prov
      break
    }
  }
  ga.pos <- which(x.sorted == ga)[1]
  tc <- NA
  if(! is.na(ga.pos)) {
    for(j in ga.pos:length(x.sorted)) {
      tc.prov <- x.sorted[j]
      if(!is.na(tc.prov) & tc.prov > ga) {
        tc <- tc.prov
        break
      }
    }
  }
  return(c(ga, tc))
}

abDistance <- function(a) {
  v1 <- c(0,0) - c(1000,1000)
  v2 <- a - c(0,0)
  m <- cbind(v1, v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  return(d)
}

theme_set(theme_bw())
my_theme <- theme(plot.title=element_text(face="bold", size=rel(1.8), hjust=0.5), text=element_text(), panel.background=element_rect(colour=NA), plot.background=element_rect(colour=NA), axis.title=element_text(face="bold", size=rel(1)), axis.title.y=element_text(angle=90, vjust=2, size=18), axis.title.x=element_text(angle=0, vjust=-2, size=18), axis.text=element_text(size=16), axis.line=element_line(colour="black"), axis.ticks=element_line(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.key=element_rect(colour=NA), legend.position="bottom", legend.direction="horizontal", legend.key.size= unit(0.2, "cm"), legend.title=element_blank(), legend.text=element_text(size=16), plot.margin=unit(c(10,5,5,5), "mm"), strip.background=element_rect(colour="#f0f0f0", fill="#f0f0f0"), strip.text=element_text(face="bold", size=12))

# Dataset Formatting
dat1 <- read.csv(args$in1, sep="\t", col.names=c("Target", "Query", "sScore", "dScore"))
dat1$DB <- "New"
dat2 <- read.csv(args$in2, sep="\t", col.names=c("Target", "Query", "sScore", "dScore"))
dat2$DB <- "Pfam"

merged <- rbind(dat1, dat2)

merged.l <- melt(merged, id.vars=c("Target", "Query", "DB"), measure.vars=c("sScore", "dScore"), value.name="Score", variable.name="Region")

merged.w <- merged.w <- dcast(merged.l, Target + Region ~ DB, value.var='Score')
merged.w <- na.omit(merged.w)

merged.w$Dist <- apply(merged.w[, c("New", "Pfam")], 1, FUN=abDistance)

abrange <- 2 * sd(merged.w$Dist)

# Statistics
nmatch.new <- length(setdiff(dat1$Target, dat2$Target))
nmatch.pfam <- length(setdiff(dat2$Target, dat1$Target))
nmatch.both <- length(unique(merged.w$Target))
nmatch.total <- length(union(dat1$Target, dat2$Target))

fp.seq <- max(subset(merged.w, Region == "sScore" & New <= Pfam)$New)  #highest False Positive sequence bitscore
fp.dom <- max(subset(merged.w, Region == "dScore" & New <= Pfam)$New)  #highest False Positive domain bitscore
thresh.seq <- find_thresholds(dat1$sScore, fp.seq)  #suggested sequence thresholds
thresh.dom <- find_thresholds(dat1$dScore, fp.dom)  #suggested domain thresholds

# Graphing
png(args$png, height=600, width=1000)
ggplot(merged.w, aes(x=New, y=Pfam)) + geom_point(data=subset(merged.w, Dist >= abrange), size=4) + geom_point(data=subset(merged.w, Dist < abrange), color="darkgrey", size=4) + geom_abline(slope=1, intercept=0, color="red") + facet_wrap(~Region, 1, scales="free", labeller=as_labeller(c("sScore"="Sequence", "dScore"="Domain"))) + labs(x="Match to custom profile (bits)", y="Match to Pfam database (bits)") + my_theme
dev.off()

# Output
cat(paste("Found", nmatch.total, "total hits, of which -\n"))
cat(paste("\t", nmatch.both, "were hits to both\n"))
cat(paste("\t", nmatch.new, "were hits unique to custom profile\n"))
cat(paste("\t", nmatch.pfam, "were hits unique to Pfam database\n"))
cat("Suggested scoring thresholds are -\n")
cat(paste("\tGA:", thresh.seq[1], thresh.dom[1], "\n"))
cat(paste("\tTC:", thresh.seq[2], thresh.dom[2], "\n"))
cat(paste("\tNC:", fp.seq, fp.dom, "\n"))
cat("\n")

if(!is.na(args$out)) {
  good_seqs <- unique(subset(merge(dat1, merged.w[,c("Target", "Dist")], by="Target"), Dist >= abrange & dScore > thresh.dom[1] & sScore > thresh.seq[1])$Target)
  write(as.character(good_seqs), file=args$out, sep="\n")
}

quit(save="no", status=0, runLast=FALSE)
