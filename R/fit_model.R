library(HAUDI)
library(bigstatsr)
library(data.table)
library(optparse)

option_list <- list(
  make_option(c("--bk_file"), type = "character"),
  make_option(c("--info_file"), type = "character"),
  make_option(c("--dims_file"), type = "character"),
  make_option(c("--gamma"), type = "integer"),
  make_option(c("--phenotype_file"), type = "character"),
  make_option(c("--model_file"), type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

code_dosage <- rep(NA_real_, 256)
code_dosage[1:201] <- seq(0, 2, length.out = 201)

fbm_info <- fread(opt$info_file)
dt_bk_dims <- fread(opt$dims_file)

bk_file <- strsplit(opt$bk_file, ".bk")[[1]]
fbm_obj <- FBM.code256(
  nrow = dt_bk_dims$nrow,
  ncol = dt_bk_dims$ncol,
  code = code_dosage,
  backingfile = bk_file,
  create_bk = FALSE,
)

y <- as.numeric(readLines(opt$phenotype_file))

mod <- HAUDI::gaudi(
  fbm_obj = fbm_obj,
  fbm_info = fbm_info,
  y = y,
  gamma_vec = opt$gamma,
  k = 5
)

saveRDS(mod, file = opt$model_file)
