# Rules for visualizing ancestry hmm results

rule ahmm_plot:
    """
    Plots general results for ahmm from .posterior files
    """
    # input:
    #     i_dir = directory(ahmmm_output_dir),
    # output:
    #     png = get_output("ahmm_plot","_{samples}/overall_density.png") # This is the name given to the final plot, follow this or
    # chromosome_densities.png, genomic_ancestry_manhattan.png, ternary_density.png
    # log:
    #     get_log_wild("ahmm_plot","{sample}_success")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    params:
        r_bin = "Rscript",
        script = "workflow/scripts/ahmm/plot_ahmm.R",
        o_prefix = lambda w, output: os.path.dirname(output.png),
        title = "Genome"
    conda:
        get_env("r_plots")
    shell:
        """
        
        export OMP_NUM_THREADS=1
        export OPENBLAS_NUM_THREADS=1
        export MKL_NUM_THREADS=1
        export VECLIB_MAXIMUM_THREADS=1
        export NUMEXPR_NUM_THREADS=1
        
        {params.r_bin} \
        {params.script} \
        -i {input.i_dir} \
        -o {params.o_prefix} \
        -l {params.title} \
        > {log} 2>&1
        """
        
