# Rules for running ASTRAL to infer species trees from gene trees using the coalescent model.
# This can be window trees obtained with genomics_general module.

rule astral_coalescent:
    """
    This rule will run ASTRAL to infer a species tree from a set of gene trees using the coalescent model.
    """
    # input:
    #     trees = "trees.trees"
    # output:
    #     tre = "results/astral_coalescent/{prefix}.tre"
    # log: 
    #    "results/logs/astral_coalescent/{prefix}.log"
    threads: 4
    resources:
        mem_mb = lambda wildcards, threads: threads * 1000
    conda: 
        get_env("astral4")
    shell:
        """
        astral4 \
        -t {threads} \
        -o {output.tre} \
        {input.trees} \
        2> {log}
        """
        
rule astral_coalescent_map:
    """
    This rule will generate a coalescent tree from .trees file.
    It will use the population/species mapping in the .map file.
    """
    input:
        trees = "input/{prefix}.trees",
        pm_hap = "input/{prefix}.map"
    # output:
    #     tre = get_output("astral_coalescent_map",".tre")
    # log:
    #     get_log("astral_coalescent_map")
    threads: 4
    conda:
        get_env("astral4")
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    shell:
        """
        astral4 \
        -t {threads} \
        -a {input.pm_hap} \
        -i {input.trees} \
        -o {output.tre} \
        2>{log}
        """

