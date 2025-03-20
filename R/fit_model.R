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
  make_option(c("--phenotype"), type = "character"),
  make_option(c("--fbm_samples_file"), type = "character", default = NULL),
  make_option(c("--training_samples_file"), type = "character", default = NULL),
  make_option(c("--model_file"), type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))

## Specify coding scheme for FBM
code_dosage <- rep(NA_real_, 256)
code_dosage[1:201] <- seq(0, 2, length.out = 201)

## Read in data
fbm_info <- fread(opt$info_file)
dt_bk_dims <- fread(opt$dims_file)

## Create FBM
bk_file <- strsplit(opt$bk_file, ".bk")[[1]]
fbm_obj <- FBM.code256(
  nrow = dt_bk_dims$nrow,
  ncol = dt_bk_dims$ncol,
  code = code_dosage,
  backingfile = bk_file,
  create_bk = FALSE,
)
fbm_samples <- readLines(opt$fbm_samples_file)

## Create phenotype
pheno <- fread(opt$phenotype_file)
pheno <- pheno[match(pheno[["IID"]], fbm_samples), ]
y <- pheno[[opt$phenotype]]
if (!is.null(opt$training_samples_file)) {
  training_samples <- readLines(opt$training_samples_file)
  y[!fbm_samples %in% training_samples] <- NA
}
ind_train <- which(!is.na(y))

## Fit model
mod <- HAUDI::gaudi(
  fbm_obj = fbm_obj,
  fbm_info = fbm_info,
  y = y,
  ind_train = ind_train,
  gamma_vec = opt$gamma,
  k = 5
)

saveRDS(mod, file = opt$model_file)
