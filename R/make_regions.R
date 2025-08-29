library(optparse)
library(data.table)
option_list <- list(
  make_option(c("--output", type = "character")),
  make_option(c("--target_snps", type = "character")),
  make_option(c("--flare_snps", type = "character")),
  make_option(c("--subset_pvar", type = "character", default = NULL))
)
opt <- parse_args(OptionParser(option_list = option_list))
snps_flare <- fread(opt$flare_snps)
snps_target <- fread(opt$target_snps)
colnames(snps_flare) <- c("#CHROM", "POS", "ID", "REF", "ALT")
colnames(snps_target) <- c("#CHROM", "POS", "ID", "REF", "ALT")

chr_flare <- gsub("chr", "", snps_flare$`#CHROM`)
chr_target <- gsub("chr", "", snps_target$`#CHROM`)

if (!is.null(opt$subset_pvar)) {
  snps_subset <- fread(opt$subset_pvar)
  chr_subset <- gsub("chr", "", snps_subset$`#CHROM`)
  snps_target <- snps_target[paste(chr_target, POS) %in% paste(chr_subset, snps_subset$POS), ]
}


snps_target$source <- "target"
snps_flare$source <- "flare"

dt <- rbind(snps_target, snps_flare)
dt <- dt[!duplicated(dt[, paste(`#CHROM`, POS)]), ]
dt <- setorder(dt, `#CHROM`, POS)

prev_and_cur_flare <- (dt$source == "flare")[2:nrow(dt)] & (dt$source == "flare")[1:(nrow(dt) - 1)]
keep <- c(TRUE, !prev_and_cur_flare)
dt <- dt[keep, ]
dt <- dt[, .(`#CHROM`, POS)]
fwrite(dt, file = opt$output, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
