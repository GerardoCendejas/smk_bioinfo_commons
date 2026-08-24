# PREFIX: dsuite

# Rule definition for running the dsuite
# In this context .tre refers to a single tree file in newick format

rule dsuite_tree_preprocess:
    """
    This rule will preprocess the .tre file for dsuite:
    - Removes node labels
    - Leaves one Outgroup tip (needed for f-branch)
    """
    # input:
    #     tre = "input/{prefix}.tre"
    # output:
    #     tre = get_output("dsuite_tree_preprocess","_{sample}.tre")
    # log:
    #     get_log_wild("dsuite_tree_preprocess","{sample}")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    params:
        r_bin = "Rscript",
        script = "workflow/scripts/dsuite/tree_preprocess.R",
        # outgroup_names = "Outgroup1,Outgroup2,Outgroup3" # Comma separated list of outgroup names
    conda:
        get_env("r_plots")
    shell:
        """
        {params.r_bin} \
        {params.script} \
        -i {input.tre} \
        -o {output.tre} \
        -g {params.outgroup_names} \
        > {log} 2>&1
        """

rule dsuite_dtrios_parallel_tree:
    """
    This rule will run the dsuite Dtrios parallel command.
    """
    # input:
    #     tre= "input/{sample}.tre",
    #     vcf= "input/{sample}.vcf.gz",
    #     pop_map= "input/{sample}.groups.tsv" # Pop map file with outgroup named as Outgroup
    # output:
    #     bbaa = get_output("dsuite_dtrios_parallel","/{sample}_BBAA.txt"),
    #     dmin = get_output("dsuite_dtrios_parallel","/{sample}_Dmin.txt"),
    #     f_branch = get_output("dsuite_dtrios_parallel","/{sample}_tree.txt")
    # log:
    #     get_log_wild("dsuite_dtrios_parallel","{sample}")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("dsuite")
    params:
        jack_rep = 100,
    shell:
        """
        TRE_ABS=$(readlink -f {input.tre})
        VCF_ABS=$(readlink -f {input.vcf})
        POPMAP_ABS=$(readlink -f {input.pop_map})
        LOG_ABS=$(readlink -f {log})

        OUT_BBAA=$(readlink -f {output.bbaa})
        OUT_DMIN=$(readlink -f {output.dmin})
        OUT_FBRANCH=$(readlink -f {output.f_branch})

        WORK=$(mktemp -d -p "{resources.tmpdir}" dtp.XXXXXXXX)
        trap 'rm -rf "$WORK"' EXIT INT TERM

        cd "$WORK"
        cp "$POPMAP_ABS" local_groups.tsv
        
        DtriosParallel \
        -t "$TRE_ABS" \
        -k {params.jack_rep} \
        --cores {threads} \
        local_groups.tsv \
        "$VCF_ABS" \
        > $LOG_ABS 2>&1


        mv ./*combined_BBAA.txt "$OUT_BBAA"
        mv ./*combined_Dmin.txt "$OUT_DMIN"
        mv ./*combined_tree.txt "$OUT_FBRANCH"

        """


rule dsuite_fbranch:
    """
    This rule will run the Fbranch command of the Dsuite tool
    """
    # input:
    #     tre = "input/{sample}.tre" # From dsuite_tree_preprocess
    #     f_branch = "input/{sample}_tree.txt" # From dsuite_dtrios_parallel_tree
    # output:
    #     f_branch_o = get_output("dsuite_fbranch","_{sample}.txt") # Dsuite Fbranch output
    # log:
    #     get_log_wild("dsuite_fbranch","{sample}")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda: 
        get_env("dsuite")
    shell:
        """
        Dsuite \
        Fbranch \
        --Zb-matrix \
        --Pb-matrix \
        {input.tre} \
        {input.f_branch} \
        > {output} \
        2> {log}
        """
rule dsuite_fbranch_plot:
    """
    This rule will plot the results of running dsuite_dtrios
    """
    # input:
    #     f_branch_o = "input/{sample}.txt",
    #     tre = "input/{sample}.tre"
    # output:
    #     png = get_output("dsuite_fbranch_plot","_{sample}.png")
    # log:
    #     get_log_wild("dsuite_fbranch_plot","{sample}")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("dsuite")
    shell:
        """
        FBRANCH_ABS=$(readlink -f {input.f_branch_o})
        TREE_ABS=$(readlink -f {input.tre})

        OUTPUT_ABS=$(readlink -f {output.png})
        OUTPUT_DIR=$(dirname "$OUTPUT_ABS"")

        LOG_ABS=$(readlink -f {log})
        
        WORK=$(mktemp -d -p "{resources.tmpdir}" dfp.XXXXXXXX)
        trap 'rm -rf "$WORK"' EXIT INT TERM

        cd "$WORK"

        sed -n '/^# Z-scores:/q;p' "$FBRANCH_ABS" > temp_filtered_fbranch.txt
        
        dtools.py \
        temp_filtered_fbranch.txt \
        "$TREE_ABS" \
        > "$LOG_ABS" 2>&1

        mv fbranch.png "$OUTPUT_ABS"

        """
rule dsuite_fbranch_plot_net:
    """
    This rule plots the fbranch statistics of the terminal tips as a network (signifficant p-values only)
    """
    # input:
    #     f_branch_o = "input/{prefix}.txt"
    # output:
    #     png = get_output("dsuite_fbranch_plot_net","_{sample}.png")
    # log:
    #     get_log_wild("dsuite_fbranch_plot_net","{sample}")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("r_plots")
    params:
        r_bin = "Rscript",
        script = "workflow/scripts/dsuite/plot_fbranch_graph.R",
        igraph_source = "workflow/scripts/dsuite/igraphplot2.R",
        p_thresh = 0.01
    shell:
        """
        {params.r_bin} \
        {params.script} \
        -i {input.f_branch_o} \
        -o {output.png} \
        -p {params.p_thresh} \
        -g {params.igraph_source} \
        > {log} 2>&1
        """
rule dsuite_dinvestigate:
    """
    This rule will run the Dinvestigate on different trios
    """
    # input:
    #     vcf = "input/{sample}.vcf.gz",
    #     pop_map= "input/{sample}.tsv"      
    # output:
    #     dinv = get_output("dsuite_dinvestigate","_{sample}.txt") # To be used in plotting
    # log:
    #     get_log_wild("dsuite_dinvestigate","{sample}")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("dsuite")
    params:
        window_size = 50,
        step_size = 25,
        # p1 = trio["p1"],
        # p2 = trio["p2"],
        # p3 = trio["p3"]
    shell:
        """
        VCF_ABS=$(readlink -f {input.vcf})
        POPMAP_ABS=$(readlink -f {input.pop_map})

        OUTPUT_ABS=$(readlink -f {output.dinv})
        OUTPUT_DIR=$(dirname $OUTPUT_ABS)

        LOG_ABS=$(readlink -f {log})

        TMP_OUT=$(basename {output.dinv}).tmp
        
        WORK=$(mktemp -d -p "{resources.tmpdir}" dfp.XXXXXXXX)
        trap 'rm -rf "$WORK"' EXIT INT TERM

        cd "$WORK"

        TRIOS_FILE="temp_trio_${{TMP_OUT}}.txt"
        echo -e "{params.p1}\t{params.p2}\t{params.p3}" > "$TRIOS_FILE"
        
        Dsuite \
        Dinvestigate \
        -w {params.window_size},{params.step_size} \
        --run-name $TMP_OUT \
        "$VCF_ABS" \
        "$POPMAP_ABS" \
        "$TRIOS_FILE" \
        > "$LOG_ABS" 2>&1

        mv *localFstats*${{TMP_OUT}}* "$OUTPUT_ABS"

        """

rule dsuite_dinvestigate_plot:
    """
    This rule plots the results of the dinvestigate by chr in an output folder
    """
    # input:
    #     dinv = "input/{sample}.txt"
    #     txt = "input/{sample}.txt" # Chromosome plotting order
    # output:
    #     # Must be success.log inside the dir where you want the plots, do not use the same for different samples, gives a lot of plots, and could cause conflict if running samples in parallel (might not run if finds already success.log in dir)
    #     log = get_output("dsuite_dinvestigate_plot","_{sample}/success.log")
    # log:
    #     get_log_wild("dsuite_dinvestigate_plot","{sample}")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    params:
        r_bin = "Rscript",
        script = "workflow/scripts/dsuite/plot_dinvestigate.R",
        o_prefix = lambda w, output: output.log.replace("success.log", "")
    conda:
        get_env("r_plots")
    shell:
        """
        {params.r_bin} \
        {params.script} \
        -i {input.dinv} \
        -o {params.o_prefix} \
        -c {input.txt} \
        > {log} 2>&1
        """


rule dsuite_fbranch_plot_map:
    """
    This rule plots the fbranch statistics of the terminal tips on top of the galapagos map
    """
    # input:
    #     f_branch_o = "input/{prefix}.txt" # From dsuite_fbranch
    #     tsv = "config/island_coords.tsv" # tsv file: location lon lat
    #     It assumes that ur ids are genus_species_location
    # output:
    #     png = get_output("dsuite_fbranch_plot_map","_{sample}.png")
    # log:
    #     get_log_wild("dsuite_fbranch_plot_map","{sample}")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("dsuite")
    params:
        r_bin = "Rscript",
        script = "workflow/scripts/dsuite/plot_fbranch_map.R",
        p_thresh = 0.01,
        o_prefix = lambda w, output: output.png.replace(".png", "")
    shell:
        """
        {params.r_bin} \
        {params.script} \
        -i {input.f_branch_o} \
        -c {input.tsv} \
        -o {params.o_prefix} \
        -p {params.p_thresh} \
        > {log} 2>&1
        """
