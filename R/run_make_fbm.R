library(optparse)
library(data.table)

option_list <- list(
  make_option(c("--vcf_file"), type = "character"),
  make_option(c("--fbm_pref"), type = "character"),
  make_option(c("--geno_format"), type = "character"),
  make_option(c("--anc_names"), type = "character"),
  make_option(c("--chunk_size"), type = "integer"),
  make_option(c("--min_ac"), type = "integer")
)
opt <- parse_args(OptionParser(option_list = option_list))

anc_names <- unlist(strsplit(opt$anc_names, ","))

result <- HAUDI::make_fbm(
  vcf_file = opt$vcf_file,
  fbm_pref = opt$fbm_pref,
  geno_format = opt$geno_format,
  anc_names = anc_names,
  chunk_size = opt$chunk_size,
  min_ac = opt$min_ac
)

dt_bk_dims <- data.table(
  ncol = ncol(result$FBM),
  nrow = nrow(result$FBM)
)

fwrite(dt_bk_dims, file = paste0(opt$fbm_pref, "_dims.txt"))
fwrite(result$info, file = paste0(opt$fbm_pref, "_info.txt"))
writeLines(
  text = attr(x = result$FBM, which = "samples"),
  con = paste0(opt$fbm_pref, "_samples.txt")
)
