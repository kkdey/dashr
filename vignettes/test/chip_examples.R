

##################  dashr on chipexo and chipseq data  #######################

tab <-  read.table("chipexo_examples/example_region1_wgEncodeSydhHistoneK562H3k27me3bUcdAlnRep1_counts.txt", header = F, stringsAsFactors = F)


extract_counts_CTCF <- function(filename){
  bed_counts <- read.table(filename, header = F, stringsAsFactors = F)
  colnames(bed_counts) <- c("chr", "start", "end", "name", "width", "counts")

  counts <- strsplit(bed_counts$counts, split = ",")[[1]]
  counts[counts == "NA"] <- 0
  counts <- as.numeric(counts)

  return(counts.l = list(chr = bed_counts$chr, start = bed_counts$start, end = bed_counts$end, counts = counts))
}

extract_counts_hg <- function(filename){
  bed_counts <- read.table(filename, header = F, stringsAsFactors = F)
  colnames(bed_counts) <- c("chr", "start", "end", "width", "counts")

  counts <- strsplit(bed_counts$counts, split = ",")[[1]]
  counts[counts == "NA"] <- 0
  counts <- as.numeric(counts)

  return(counts.l = list(chr = bed_counts$chr, start = bed_counts$start, end = bed_counts$end, counts = counts))
}

chipexo1 <- extract_counts_hg("chipexo_examples/example_region1_wgEncodeSydhHistoneK562H3k9acbUcdAlnRep1_counts.txt")
plot(chipexo1$counts, type = "h", col = "blue", ylab = "rep1 forward", xlab = "")

out <- dashr::dash_smooth(chipexo1$counts, dash_control = list(Inf_weight = 1))
lines(out$estimate, col = "red")
out2 <- smashr::smash.poiss(chipexo1$counts)
lines(out2, col = "green")

plot(chipexo1$counts, col = "gray80", type = "h", ylab = "rep1 forward", xlab = "", main = "CTCF-rep1-forward")
lines(out2, col = "blue", lwd = 2)
lines(out$estimate, col = "red", lwd = 4)
legend("topright", # places a legend at the appropriate place
       c("truth","dash-m", "smash-poiss"), # puts text in the legend
       lty=c(1,1), # gives the legend appropriate symbols (lines)
       lwd=c(2.5,2.5),
       cex = 0.5,
       col=c("gray80","red", "blue"))

ll <- list("counts" = chipexo1$counts, "dash" = out$estimate, "smash" = out2)
save(ll, file = "chipexo_wgEncodeSydhHistoneK562H3k9acbUcdAlnRep1.rda")


