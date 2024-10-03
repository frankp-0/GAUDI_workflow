version 1.0

workflow run_gaudi {
  input {
      String vcf_file
      String fbm_pref
      String geno_format
      Array[String] anc_names
      Int chunk_size
      Int min_ac
    }

    call make_fbm {
      input:  vcf_file=vcf_file,
              fbm_pref=fbm_pref,
              geno_format=geno_format,
              anc_names=anc_names,
              chunk_size=chunk_size,
              min_ac=min_ac
    }

    output {
      File rds_file = make_fbm.rds_file
      File bk_file = make_fbm.bk_file
      File info_file = make_fbm.info_file
    }

    meta {
        author: "Frank Ockerman"
        email: "frankpo@unc.edu"
      }
  }

task make_fbm {
  input {
    String vcf_file
    String fbm_pref
    String geno_format
    Int chunk_size
  }
  command <<<
    Rscript run_make_fbm.R \
      --vcf_file=~{vcf_file} \
      --fbm_pref=~{fbm_pref} \
      --geno_format=~{geno_format} \
      --anc_names=~sep(",", ~{anc_names}) \
      --chunk_size=~{chunk_size} \
      --min_ac=~{min_ac}
  >>>
  output {
    File rds_file = "${fbm_pref}.rds"
    File bk_file = "${fbm_pref}.bk"
    File info_file = "${fbm_pref}.info"
  }
  runtime {
    docker: "frankpo/run_gaudi:0.0.1"
  }
}
