version 1.0

workflow prepare_input {
    input {
        File vcf_file
        File vcf_index_file
        File flare_vcf_file
        File flare_vcf_index_file
        File? snp_list
        String fbm_prefix
        String geno_format
        Array[String] anc_names
        Int chunk_size = 1000
        Int min_ac = 10
        Int subset_target_mem_gb = 8
        Int make_fbm_mem_gb = 16
    }

    call subset_target {
        input:
          target_vcf=vcf_file,
          target_vcf_index=vcf_index_file,
          flare_vcf=flare_vcf_file,
          flare_vcf_index=flare_vcf_index_file,
          snp_list=snp_list,
          fbm_pref=fbm_prefix,
          mem_gb=subset_target_mem_gb
        }

    call make_fbm {
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
 
    output {
        File bk_file = make_fbm.bk_file
        File info_file = make_fbm.info_file
        File dims_file = make_fbm.dims_file
        File fbm_samples_file = make_fbm.samples_file
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
        String fbm_pref
        Int mem_gb
    }

    Int disk_size = ceil(2 * (size(target_vcf, "GB") + size(target_vcf_index, "GB") + size(flare_vcf, "GB") + size(flare_vcf_index, "GB")))

    command <<<
        bcftools query -f '%ID\n' ~{target_vcf} > temp_target_snps.txt
        bcftools query -f '%ID\n' ~{flare_vcf} > temp_flare_snps.txt
        if [[ "~{snp_list}" = "" ]]
        then
          awk 'NR==FNR{a[$0]; next} !($0 in a){print $0} END{for (x in a) print x}' temp_target_snps.txt temp_flare_snps.txt > subset_snps.txt
        else
          awk 'NR==FNR{a[$0]; next} $0 in a' temp_target_snps.txt ~{snp_list} | awk 'NR==FNR{c[$0]; next} !($0 in c){print $0} END{for (x in c) print x}' temp_flare_snps.txt - > subset_snps.txt
        fi
        bcftools view --include ID==@subset_snps.txt ~{target_vcf} -Oz -o ~{fbm_pref}.vcf.gz
        bcftools index -t ~{fbm_pref}.vcf.gz
        bcftools annotate -c FORMAT -a ~{flare_vcf} ~{fbm_pref}.vcf.gz -Oz -o ~{fbm_pref}.anc.vcf.gz
        bcftools index -t ~{fbm_pref}.anc.vcf.gz
    >>>

    output {
        File subset_vcf = "~{fbm_pref}.anc.vcf.gz"
        File subset_vcf_index = "~{fbm_pref}.anc.vcf.gz.tbi"
    }
    runtime {
        docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
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
        Int mem_gb
    }

    Int disk_size = ceil(4 * (size(vcf_file, "GB"))) + 1

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
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}
