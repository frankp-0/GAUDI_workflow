FROM bioconductor/bioconductor:RELEASE_3_19@sha256:c74fd973f50818c4044225c8cfcf318fa11087d5df1b6a08a5a604abbe91195a

RUN R -e 'BiocManager::install("VariantAnnotation");\
    devtools::install_github("https://github.com/frankp-0/HAUDI.git");\
    install.packages("optparse")'

COPY R/run_make_fbm.R /scripts/run_make_fbm.R
