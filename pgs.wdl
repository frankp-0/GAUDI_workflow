version 1.0

workflow pgs {
    input {
        File bk_file
        File info_file
        File dims_file
        File fbm_samples_file
        Int gamma
        File phenotype_file
        String phenotype
        String output_prefix
        Int fit_gaudi_mem_gb = 20
        File? snps_file
    }

    call fit_gaudi {
        input:
            bk_file=bk_file,
            info_file=info_file,
            dims_file=dims_file,
            fbm_samples_file=fbm_samples_file,
            gamma=gamma,
            phenotype_file=phenotype_file,
            phenotype=phenotype,
            output_prefix=output_prefix,
            mem_gb=fit_gaudi_mem_gb,
            snps_file=snps_file
    }

    output {
        File model_file = fit_gaudi.model_file
        File pgs_file = fit_gaudi.pgs_file
        File effects_file = fit_gaudi.effects_file
    }

    meta {
        author: "Frank Ockerman"
        email: "frankpo@unc.edu"
    }
}

task fit_gaudi {
    input {
        File bk_file
        File info_file
        File dims_file
        File fbm_samples_file
        Float gamma
        File phenotype_file
        String phenotype
        String output_prefix
        Int mem_gb
        File? snps_file
    }

    Int disk_size = ceil(size(bk_file, "GB") + size(info_file, "GB") + size(dims_file, "GB") + size(fbm_samples_file, "GB") + size(phenotype_file, "GB") + 2)

    command <<<
        Rscript /scripts/fit_model.R \
          --bk_file ~{bk_file} \
          --info_file ~{info_file} \
          --dims_file ~{dims_file} \
          --fbm_samples_file ~{fbm_samples_file} \
          --gamma ~{gamma} \
          --phenotype_file ~{phenotype_file} \
          --phenotype ~{phenotype} \
          --output_prefix ~{output_prefix} \
          ~{if defined(snps_file) then "--snps_file " + snps_file else ""}
    >>>

    output {
        File model_file = "${output_prefix}_model.rds"
        File effects_file = "${output_prefix}_effects.txt"
        File pgs_file = "${output_prefix}_pgs.txt"
    }


    runtime {
        docker: "frankpo/run_gaudi:0.0.4"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}
