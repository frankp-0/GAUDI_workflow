FROM bioconductor/bioconductor:RELEASE_3_19@sha256:c74fd973f50818c4044225c8cfcf318fa11087d5df1b6a08a5a604abbe91195a

RUN R -e 'BiocManager::install("VariantAnnotation");\
    devtools::install_github("frankp-0/HAUDI@f5a060806267e1097cc9bb787534c1219e3136e9");\
    install.packages("optparse")'

RUN git clone --recurse-submodules https://github.com/samtools/htslib.git \ 
    && git clone https://github.com/samtools/bcftools.git \
    && cd bcftools \
    && make

COPY R/* /scripts/
COPY test_data /test_data
