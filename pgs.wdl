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
            output_prefix=output_prefix
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
    }

    command <<<
        Rscript /scripts/fit_model.R \
          --bk_file ~{bk_file} \
          --info_file ~{info_file} \
          --dims_file ~{dims_file} \
          --fbm_samples_file ~{fbm_samples_file} \
          --gamma ~{gamma} \
          --phenotype_file ~{phenotype_file} \
          --phenotype ~{phenotype} \
          --output_prefix ~{output_prefix}
    >>>

    output {
        File model_file = "${output_prefix}_model.rds"
        File effects_file = "${output_prefix}_effects.txt"
        File pgs_file = "${output_prefix}_pgs.txt"
    }

    runtime {
        docker: "frankpo/run_gaudi:0.0.4"
    }
}
