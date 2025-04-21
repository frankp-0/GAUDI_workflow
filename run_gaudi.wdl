version 1.0

import "prepare_input.wdl" as prep
import "pgs.wdl" as pgs

workflow run_gaudi {
    input {
        File vcf_file
        File vcf_index_file
        File flare_vcf_file
        File flare_vcf_index_file
        Int min_ac
        File? snp_list
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
    }

    call prep.subset_target {
        input:
          target_vcf=vcf_file,
          target_vcf_index=vcf_index_file,
          flare_vcf=flare_vcf_file,
          flare_vcf_index=flare_vcf_index_file,
          snp_list=snp_list,
          fbm_pref=fbm_prefix,
          mem_gb=subset_target_mem_gb
    }

    call prep.make_fbm {
        input:
            vcf_file=subset_target.subset_vcf,
            vcf_index_file=subset_target.subset_vcf_index,
            fbm_pref=fbm_prefix,
            geno_format=geno_format,
            anc_names=anc_names,
            chunk_size=chunk_size,
            min_ac=min_ac,
            mem_gb=make_fbm_mem_gb
    }
 
    call pgs.fit_gaudi {
        input:
            bk_file=make_fbm.bk_file,
            info_file=make_fbm.info_file,
            dims_file=make_fbm.dims_file,
            fbm_samples_file=make_fbm.samples_file,
            gamma=gamma,
            phenotype_file=phenotype_file,
            phenotype=phenotype,
            output_prefix=output_prefix,
            mem_gb=fit_gaudi_mem_gb
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
