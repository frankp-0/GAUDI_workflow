version 1.0

workflow prepare_input {
    input {
        Array[File] vcf_files
        Array[File] vcf_index_files
        Array[File] flare_vcf_files
        File? samples
        File? fbm_subset_pvar
        Int subset_target_mem_gb = 8
        String fbm_prefix
        String geno_format
        Array[String] anc_names
        Int chunk_size = 1000
        Int min_ac = 10
        Int make_fbm_mem_gb = 16
    }

    scatter (i in range(length(vcf_files))) {
        call get_regions as get_regions_chrom {
            input:
                target_vcf=vcf_files[i],
                flare_vcf=flare_vcf_files[i],
                fbm_subset_pvar=fbm_subset_pvar,
                mem_gb=subset_target_mem_gb
          }

        call make_anc_vcf as make_anc_vcf_chrom {
            input:
                regions=get_regions_chrom.regions,
                target_vcf=vcf_files[i],
                target_vcf_index=vcf_index_files[i],
                flare_vcf=flare_vcf_files[i],
                samples=samples
          }
      }

    call concat_vcfs {
        input:
            vcfs = make_anc_vcf_chrom.vcf_file,
            vcf_indices = make_anc_vcf_chrom.vcf_index,
            fbm_pref = fbm_prefix
    }

    call make_fbm {
        input:
            vcf_file=concat_vcfs.vcf,
            vcf_index_file=concat_vcfs.index,
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

task get_regions {
    input {
        File target_vcf
        File flare_vcf
        File? fbm_subset_pvar
        Int mem_gb
    }

    Int disk_size = ceil(4 * (size(target_vcf, "GB") + size(flare_vcf, "GB")))
    String raw_name = basename(target_vcf)
    String base_name = select_first([
        sub(raw_name, "\\.vcf\\.gz$", ""),
        sub(raw_name, "\\.vcf$", "")
    ])

     command <<<
        /bcftools/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ~{target_vcf} > target.pvar
        /bcftools/bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ~{flare_vcf} > flare.pvar

        Rscript /scripts/make_regions.R \
            --output ~{base_name}.regions.txt \
            --target_snps target.pvar \
            --flare_snps flare.pvar \
            ~{if defined(fbm_subset_pvar) then "--subset_pvar " + fbm_subset_pvar else ""}
    >>>

    output {
        File regions = "~{base_name}.regions.txt"
        File target_pvar = "target.pvar"
        File flare_pvar = "flare.pvar"
    }
    runtime {
        docker: "frankpo/run_gaudi:0.0.5"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}

task make_anc_vcf {
    input {
        File regions
        File target_vcf
        File target_vcf_index
        File flare_vcf
        File? samples
    }

    Int disk_size = ceil(4 * (size(target_vcf, "GB") + size(flare_vcf, "GB")))
    String raw_name = basename(target_vcf)
    String base_name = select_first([
        sub(raw_name, "\\.vcf\\.gz$", ""),
        sub(raw_name, "\\.vcf$", "")
    ])


    command <<<
        /bcftools/bcftools view -R ~{regions} \
          ~{if defined(samples) then "-S " + samples else ""} \
          ~{target_vcf} -Oz -o subset_target.vcf.gz
        tabix subset_target.vcf.gz
        tabix ~{flare_vcf}
        /bcftools/bcftools view -R ~{regions} \
          ~{if defined(samples) then "-S " + samples else ""} \
          ~{flare_vcf} -Oz -o subset_flare.vcf.gz
        tabix subset_flare.vcf.gz
        /bcftools/bcftools annotate -c FORMAT -a subset_flare.vcf.gz subset_target.vcf.gz -Oz -o ~{base_name}.lanc.vcf.gz
        tabix ~{base_name}.lanc.vcf.gz
    >>>

    output {
        File vcf_file = "~{base_name}.lanc.vcf.gz"
        File vcf_index = "~{base_name}.lanc.vcf.gz.tbi"
    }

    runtime {
        docker: "frankpo/run_gaudi:0.0.5"
        disks: "local-disk ~{disk_size} SSD"
        memory: "4G"
    }
}

task concat_vcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_indices
        String fbm_pref
    }

    Int disk_size = ceil(2.5 * size(vcfs, "GB"))

    command <<<
        /bcftools/bcftools concat ~{sep=' ' vcfs} -Oz -o ~{fbm_pref}.vcf.gz
        tabix ~{fbm_pref}.vcf.gz
    >>>

    output {
        File vcf = "~{fbm_pref}.vcf.gz"
        File index = "~{fbm_pref}.vcf.gz.tbi"
    }

    runtime {
        docker: "frankpo/run_gaudi:0.0.5"
        disks: "local-disk ~{disk_size} SSD"
        memory: "4G"
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
        docker: "frankpo/run_gaudi:0.0.5"
        disks: "local-disk ~{disk_size} SSD"
        memory: "~{mem_gb}G"
    }
}
