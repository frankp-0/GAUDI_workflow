version 1.0

workflow run_gaudi {
  input {
      File vcf_file
      File vcf_index_file
      File flare_vcf_file
      File flare_vcf_index_file
      File? snp_list
      String fbm_pref
      String geno_format
      Array[String] anc_names
      Int chunk_size
      Int min_ac
      Int gamma
      File phenotype_file
      String phenotype
      String output_prefix
    }

    call subset_target {
      input:
        target_vcf=vcf_file,
        target_vcf_index=vcf_index_file,
        flare_vcf=flare_vcf_file,
        flare_vcf_index=flare_vcf_index_file,
        snp_list=snp_list
      }

    call make_fbm {
      input:  vcf_file=subset_target.subset_vcf,
              vcf_index_file=subset_target.subset_vcf_index,
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
             fbm_samples_file=make_fbm.samples_file,
             gamma=gamma,
             phenotype_file=phenotype_file,
             phenotype=phenotype,
             output_prefix=output_prefix
      }

    output {
      File bk_file = make_fbm.bk_file
      File info_file = make_fbm.info_file
      File dims_file = make_fbm.dims_file
      File fbm_samples_file = make_fbm.samples_file
      File model_file = fit_gaudi.model_file
      File pgs_file = fit_gaudi.pgs_file
      File effects_file = fit_gaudi.effects_file
    }

    meta {
        author: "Frank Ockerman"
        email: "frankpo@unc.edu"
      }
  }

task subset_target {
    input {
        File target_vcf
        File target_vcf_index
        File flare_vcf
        File flare_vcf_index
        File? snp_list
  }
  command <<<
    bcftools query -f '%ID\n' ~{target_vcf} > temp_target_snps.txt
    bcftools query -f '%ID\n' ~{flare_vcf} > temp_flare_snps.txt
    if [[ "~{snp_list}" = "" ]]
    then
      awk 'NR==FNR{a[$0]; next} !($0 in a){print $0} END{for (x in a) print x}' temp_target_snps.txt temp_flare_snps.txt > subset_snps.txt
    else
      awk 'NR==FNR{a[$0]; next} $0 in a' temp_target_snps.txt ~{snp_list} | awk 'NR==FNR{c[$0]; next} !($0 in c){print $0} END{for (x in c) print x}' temp_flare_snps.txt - > subset_snps.txt
    fi
    bcftools view --include ID==@subset_snps.txt ~{target_vcf} -Oz -o subset_target.vcf.gz
    bcftools index -t subset_target.vcf.gz
    bcftools annotate -c FORMAT -a ~{flare_vcf} subset_target.vcf.gz -Oz -o subset.anc.vcf.gz
    bcftools index -t subset.anc.vcf.gz
  >>>

  output {
      File subset_vcf = "subset.anc.vcf.gz"
      File subset_vcf_index = "subset.anc.vcf.gz.tbi"
    }
  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
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
    File samples_file = "${fbm_pref}_samples.txt"
  }

  runtime {
      docker: "frankpo/run_gaudi:0.0.4"
    }
}

task fit_gaudi {
  input {
      String fbm_pref
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
