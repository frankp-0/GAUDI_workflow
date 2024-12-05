version 1.0

workflow run_gaudi {
  input {
      File vcf_file
      File vcf_index_file
      String fbm_pref
      String geno_format
      Array[String] anc_names
      Int chunk_size
      Int min_ac
      Int gamma
      File phenotype_file
    }

    call make_fbm {
      input:  vcf_file=vcf_file,
              vcf_index_file=vcf_index_file,
              fbm_pref=fbm_pref,
              geno_format=geno_format,
              anc_names=anc_names,
              chunk_size=chunk_size,
              min_ac=min_ac
    }
 
    call fit_gaudi {
      input: fbm_pref=fbm_pref,
             bk_file=make_fbm.bk_file,
             info_file=make_fbm.info_file,
             dims_file=make_fbm.dims_file,
             gamma=gamma,
             phenotype_file=phenotype_file
      }

    output {
      File bk_file = make_fbm.bk_file
      File info_file = make_fbm.info_file
      File dims_file = make_fbm.dims_file
      File model_file = fit_gaudi.model_file
    }

    meta {
        author: "Frank Ockerman"
        email: "frankpo@unc.edu"
      }
  }

task make_fbm {
  input {
      File vcf_file
      File vcf_index_file
      String fbm_pref
      String geno_format
      Array[String] anc_names
      Int chunk_size
      Int min_ac
  }
  command <<<
    Rscript /scripts/run_make_fbm.R \
      --vcf_file ~{vcf_file} \
      --fbm_pref ~{fbm_pref} \
      --geno_format ~{geno_format} \
      --anc_names ~{sep=',' anc_names} \
      --chunk_size ~{chunk_size} \
      --min_ac ~{min_ac}
  >>>

  output {
    File bk_file = "${fbm_pref}.bk"
    File dims_file = "${fbm_pref}_dims.txt"
    File info_file = "${fbm_pref}_info.txt"
  }

  runtime {
      docker: "frankpo/run_gaudi:0.0.1"
    }
}

task fit_gaudi {
  input {
      String fbm_pref
      File bk_file
      File info_file
      File dims_file
      Float gamma
      File phenotype_file
    }

    command <<<
    Rscript /scripts/fit_model.R \
      --bk_file ~{bk_file} \
      --info_file ~{info_file} \
      --dims_file ~{dims_file} \
      --gamma ~{gamma} \
      --phenotype_file ~{phenotype_file} \
      --model_file ~{fbm_pref}_model.rds
    >>>

    output {
        File model_file = "${fbm_pref}_model.rds"
      }

  runtime {
      docker: "frankpo/run_gaudi:0.0.3"
    }
  }
