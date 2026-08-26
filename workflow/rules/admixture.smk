# PREFIX: admixture

# Rules for running admixture

# Similar rule in bcftools now. Redundant for the future.
rule admixture_vcf_preprocess:
    """
    Preprocess vcf/bcf to remove 'chr' string
    """
    # input:
    #     bcf = "input/{prefix}.bcf"
    # output:
    #     bcf = get_output("admixture_vcf_preprocess","_{sampe}.bcf")
    #     csi = get_output("admixture_vcf_preprocess","_{sampe}.bcf.csi")
    # log:
    #     get_log_wild("admixture_vcf_preprocess","{sampe}")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("bcftools")
    shell:
        """
        bcftools view \
        -Ou \
        {input.bcf} \
        | sed 's/^chr//' \
        | bcftools view \
        -O b \
        --write-index=csi \
        -o {output.bcf} \
        --threads {threads} \
        > {log} 2>&1
        """

## Note: This rule may be used for LD pruning in other contexts as well.
rule admixture_plink2_ld_prune:
    """
    This will do LD pruning using plink2.
    """
    # input:
    #     bcf = "input/{sample}.bcf"
    # output:
    #     bed = get_output("admixture_plink_ld_prune","_{sample}.bed"),
    #     bim = get_output("admixture_plink_ld_prune","_{sample}.bim"),
    #     fam = get_output("admixture_plink_ld_prune","_{sample}.fam")
    # log:
    #     get_log_wild("admixture_plink_ld_prune","{sample}")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("admixture")
    params:
        o_prefix = lambda w, output: output.bed.replace(".bed", ""),
        chr_num = 28, # This is important to change, this used as an example of birds
        # It is flexible for diferent organisms other than humans, but make sure to use your correct
        # number of chromosomes
        window_size = 50,
        step_size = 10,
        r2_threshold = 0.1,
        args = "" # Any additional arguments you want to give
        # E.g: --bad-ld for when using less than 50 samples, even though not recommended by plink
    shell:
        """
        plink2 \
        --bcf {input.bcf} \
        --allow-extra-chr \
        --chr-set {params.chr_num} \
        --set-all-var-ids @:# \
        --indep-pairwise {params.window_size} {params.step_size} {params.r2_threshold} \
        --threads {threads} \
        {params.args} \
        --out {params.o_prefix} \
        > {log} 2>&1

        plink2 \
        --bcf {input.bcf} \
        --allow-extra-chr \
        --chr-set {params.chr_num} \
        --set-all-var-ids @:# \
        --extract {params.o_prefix}.prune.in \
        --make-bed \
        --threads {threads} \
        --out {params.o_prefix} \
        > {log} 2>&1
        """

rule admixture_bim_preprocess:
    """
    Changes chrom name to 0 (admixture does not take non-human chroms)
    """
    # input:
    #     bim = "input/{prefix}.bim"
    #     bed = "input/{prefix}.bed"
    #     fam = "input/{prefix}.fam"
    # output:
    #     bim = get_output("admixture_bim_preprocess","_{sample}.bim")
    #     bed = get_output("admixture_bim_preprocess","_{sample}.bed")
    #     fam = get_output("admixture_bim_preprocess","_{sample}.fam")
    # log:
    #     get_log_wild("admixture_bim_preprocess","{sample}")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    shell:
        """
        awk '{{$1="0";print $0}}' \
        {input.bim} \
        > {output.bim} &&

        cp {input.bed} {output.bed} \
        > {log} 2>&1

        cp {input.fam} {output.fam} \
        > {log} 2>&1
        """

rule admixture_run:
    """
    Run admixture with n number of ancestral populations (K) and cross-validation
    """
    # input:
    #     bed = "input/{prefix}.bed"
    # output:
    #     # We expect the outputs to be in a folder, only specifying the Q file here
    #     q_file = get_output("admixture_run:","_{sample}/{k}.Q")
    # log:
    #     get_log_wild("admixture_run:","{sample}/{k}")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    params:
        # k = 
    conda:
        get_env("admixture")
    shell:
        """
        INPUT_ABS=$(readlink -f {input.bed})
        
        OUTPUT_ABS=$(readlink -f {output.q_file})
        OUTPUT_DIR=$(dirname $OUTPUT_ABS)

        LOG_ABS=$(readlink -f {log})

        mkdir -p $OUTPUT_DIR && \
        cd $OUTPUT_DIR
        
        admixture \
        -j{threads} \
        --cv $INPUT_ABS \
        {params.k} \
        > $LOG_ABS 2>&1

        mv *{params.k}.Q $OUTPUT_ABS
        """

rule admixture_plot:
    """
    Plots the admixture results by individual, cluster by population/species
    """
    # input:
    #     q_file = expand(rule.admixture_run.output, admixture_k=config["admixture_ks"]),
    #     fam = "input/{sample}.fam",
    #     # Ind Sp
    #     pop_map = "input/{sample}.tsv",
    #     # Pop order
    #     txt  = "input/{sample}.txt" # One population per line, in the order you want them to appear in the plot
    # output:
    #     png = get_output("admixture_plot","_{sample}.png")
    # log:
    #     get_log_wild("admixture_plot","{sample}")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    params:
        r_bin = "Rscript",
        script = get_script("admixture/plot_admixture.R"),
        q_dir = lambda w, input: os.path.dirname(input.q_file[0]),
        # min_k = min(config["admixture_ks"]),
        # max_k = max(config["admixture_ks"]),
    conda:
        get_env("r_plots")
    shell:
        """
        {params.r_bin} \
        {params.script} \
        --dir {params.q_dir} \
        --fam {input.fam} \
        --pop_map {input.pop_map} \
        --sort {input.txt} \
        --min  {params.min_k} \
        --max  {params.max_k} \
        --out {output.png} \
        > {log} 2>&1
        """


