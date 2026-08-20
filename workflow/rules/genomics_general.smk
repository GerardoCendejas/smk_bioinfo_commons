import os

rule genomics_general_geno:
    """
    This rule is a general rule for generating the format geno.gz expected by genomics general.
    Its input is a bcf file and its output is a geno.gz file.
    """
    # input:
    #     bcf: "input.bcf"
    # output:
    #     geno= "results/{prefix}.geno.gz"
    # log:
    #     "logs/genomic_general_geno/{prefix}.log"
    threads: 4

    resources:
        mem_mb = lambda wildcards, threads: threads * 1000
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        ( \
                echo -ne "#CHROM\\tPOS\\t" && \
                bcftools query -l {input.bcf} | tr '\\n' '\\t' | sed 's/\\t$//' && \
                echo "" && \
                bcftools query -f '%CHROM\\t%POS[\\t%TGT]\\n' {input.bcf} \
        ) 2> {log} | gzip > {output.geno}
       """


rule genomics_general_phyml_sliding_window:
    """
    Recipe: Run PhyML in sliding windows using genomics_general script.
    
    Defaults:
      - Window: 100kb
      - Step: 100kb
      - Model: GTR
      - MinSites: 200
    """
    
    # input:
    #    geno = "input.geno.gz"
    
    output:
        # trees_g = "results/genomics_general_phyml_sliding_window/{prefix}.trees.gz",
        # The tsv is also and output, but was not included when first running the rule
        # Add it here for completeness when moving to commons
        
        # tsv = "results/genomics_general_phyml_sliding_window/{prefix}.data.tsv"
    
    # log:
    #     "logs/genomics_general_phyml_sliding_window/{prefix}.log"
        
    params:
        # --- PARAMETERS (Defaults) ---
        w_size = 100000,
        s_size = 100000,
        min_sites = 200,
        wind_type = "coordinate",
        model = "GTR",
        out_prefix = lambda w, output: output.trees_g.replace(".trees.gz", "")

    threads: 4
    resources:
        mem_mb = lambda wildcards, threads: threads * 1000
    conda:
        "envs/genomics_general.yaml"
    shell:
        """
        phyml_sliding_windows.py \
            -T {threads} \
            -g {input.geno} \
            --prefix {params.out_prefix} \
            -w {params.w_size} \
            -S {params.s_size} \
            -M {params.min_sites} \
            --windType {params.wind_type} \
            --model {params.model} \
            --phyml phyml \
        > {log} 2>&1
        """
