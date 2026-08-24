## PREFIX: nwku

## These rules are general rules for plotting phylogenetic trees in Newick format.


rule nwku_plot_tree:
    """
    This is a generic rule to plot a phylogenetic tree from a Newick format file.
nn    """
    # input:
    #     tre = "input/{prefix}.tre"
    # output:
    #     png = get_output("plot_tree",".png")
    # log:
    #     get_log("plot_tree")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    params:
        r_bin = "Rscript",
        script = "workflow/scripts/newick_utils/plot_tree.R",
        rooted = "-r", # "-r" Plot rooted tree
        node_labels = "", # "-n" Plot node labels
        size = 800, # 800 Size of the output image in pixels
    conda:
        get_env("r_plots")
    shell:
        """
        {params.r_bin} \
        {params.script} \
        -i {input.tre} \
        -o {output.png} \
        -s {params.size} \
        {params.rooted} \
        {params.node_labels} \
        2> {log}
        """

rule nwku_plot_tree_color_tip:
    """
    This rule will plot a tree with color tips
    """
    # input:
    #     tre="input/{sample}.tre",
    #     pop_map="input/{sample}.csv" # Use pm_hap (popmap hap file if using haplotypes), but the name is still pop_map
    # output:
    #     png = get_output("plot_tree_color_tip","{sample}.png")
    # log:
    #     get_log("plot_tree_color_tip")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    params:
        r_bin = "Rscript",
        script = "workflow/scripts/newick_utils/plot_tree_color_tips.R",
        rooted = "-r", # "-r" Plot rooted tree
        node_labels = "", # "-n" Plot node labels
        size = 10, # Size of the output image in inches
        dpi = 300, # DPI of the output image
    conda:
        get_env("r_plots")
    shell:
        """
        {params.r_bin} \
        {params.script} \
        -i {input.tre} \
        -o {output.png} \
        -m {input.pop_map} \
        -s {params.size} \
        -d {params.dpi} \
        {params.rooted} \
        {params.node_labels} \
        2> {log}
        """
        
rule nwku_reroot_tree:
    """
    Reroot the tree given the names in a txt file or a comma separted string
    """
    # input:
    #     tre="input/{prefix}.tre"
    #     samples="input/{prefix}_outgroup.txt" # Samples file (one name per line) of comma separated string of names to use as outgroup
    # output:
    #     tre=get_output("reroot_tree",".tre")
    # log:
    #     get_log("reroot_tree")
    conda:
        get_env("r_plots")
    params:
        r_bin = "Rscript",
        script = "workflow/scripts/newick_utils/reroot_tree.R",
    shell:
        """
        {params.r_bin} \
        {params.script} \
        -i {input.tre} \
        -o {output.tre} \
        -g {input.samples} \
        2> {log}
        """


