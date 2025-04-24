library(optparse)
library(data.table)
option_list <- list(
  make_option(c("--target_snps", type = "character")),
  make_option(c("--flare_snps", type = "character")),
  make_option(c("--subset_pvar", type = "character", default = NULL))
)
opt <- parse_args(OptionParser(option_list = option_list))
snps_flare <- fread(opt$flare_snps)
snps_target <- fread(opt$target_snps)
colnames(snps_flare) <- c("CHROM", "POS", "ID", "REF", "ALT")
colnames(snps_target) <- c("CHROM", "POS", "ID", "REF", "ALT")

if (!is.null(opt$subset_pvar)) {
  snps_subset <- fread(opt$subset_pvar)
  snps_target <- snps_target[paste(CHROM, POS) %in% paste(snps_subset$CHROM, snps_subset$POS), ]
}

snps_target$source <- "target"
snps_flare$source <- "flare"

dt <- rbind(snps_target, snps_flare)
dt <- dt[!duplicated(dt[, paste(CHROM, POS)]), ]
dt <- setorder(dt, CHROM, POS)

prev_and_cur_flare <- (dt$source == "flare")[2:nrow(dt)] & (dt$source == "flare")[1:(nrow(dt) - 1)]
keep <- c(TRUE, !prev_and_cur_flare)
dt <- dt[keep, ]
dt <- dt[, .(CHROM, POS)]
fwrite(dt, file = "subset_regions.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
