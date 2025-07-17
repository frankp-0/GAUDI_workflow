version 1.0

import "prepare_input.wdl" as prep
import "pgs.wdl" as pgs

workflow run_gaudi {
    input {
        Array[File] vcf_files
        Array[File] vcf_index_files
        Array[File] flare_vcf_files
        Int min_ac
        File? fbm_subset_pvar
        String fbm_prefix
        String geno_format
        Array[String] anc_names
        Int chunk_size
        Int subset_target_mem_gb = 8
        Int make_fbm_mem_gb = 16
        Int gamma
        File phenotype_file
        String phenotype
        String output_prefix
        Int fit_gaudi_mem_gb = 20
        File? snps_file
    }

    call prep.prepare_input as prepare_input {
        input:
            vcf_files = vcf_files,
            vcf_index_files = vcf_index_files,
            flare_vcf_files = flare_vcf_files,
            fbm_subset_pvar = fbm_subset_pvar,
            subset_target_mem_gb = subset_target_mem_gb,
            fbm_prefix = fbm_prefix,
            geno_format = geno_format,
            anc_names = anc_names,
            chunk_size = chunk_size,
            min_ac = min_ac,
            make_fbm_mem_gb = make_fbm_mem_gb
    }

    call pgs.fit_gaudi {
        input:
            bk_file=prepare_input.bk_file,
            info_file=prepare_input.info_file,
            dims_file=prepare_input.dims_file,
            fbm_samples_file=prepare_input.fbm_samples_file,
            gamma=gamma,
            phenotype_file=phenotype_file,
            phenotype=phenotype,
            output_prefix=output_prefix,
            mem_gb=fit_gaudi_mem_gb,
            snps_file=snps_file
    }

    output {
      File bk_file = prepare_input.bk_file
      File info_file = prepare_input.info_file
      File dims_file = prepare_input.dims_file
      File fbm_samples_file = prepare_input.fbm_samples_file
      File model_file = fit_gaudi.model_file
      File pgs_file = fit_gaudi.pgs_file
      File effects_file = fit_gaudi.effects_file
    }

    meta {
        author: "Frank Ockerman"
        email: "frankpo@unc.edu"
    }
}
