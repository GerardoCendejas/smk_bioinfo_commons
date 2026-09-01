# PREFIX: phlash

# Set of rules for demographic inference using the Bayesian method of Phlash (Terhorst, 2025).

rule phlash_run:
    """
    Runs the phlash MCMC inference on the complete set of chromosome and samples given, gets the output in a .tsv file
    This will contain the time in generations as first column and effective population size in the other columns, (500 MCMC samples)
    """
    # input:
    #     bcf = "input/{prefix}.bcf.gz",
    #     txt_chr = "input/{prefix}.chroms.txt",
    #     samples = "input/{prefix}.samples.txt"
    # output:
    #     tsv = get_output("phlash_run","_{prefix}.tsv")
    # log:
    #     get_log_wild("phlash_run","{prefix}")
    threads: 4
    params:
        python_bin = "python",
        script = get_script("phlash/phlash_run.py"),
        m = 2.04e-9 # WARNING: This is a mutation rate per site per generation, for birds, change according to your species.
    resources:
        mem_mb=lambda wildcards, threads, params: threads * 1000 
    conda:
        get_env("phlash")
    shell:
        """
        export XLA_PYTHON_CLIENT_PREALLOCATE=false

        # Phlash uses Jax library, does not care about Snakemake core management
        # Will use systemd limits to force the use of only N cores
        systemd-run --user --scope -p CPUQuota={threads}00% \
	{params.python_bin} \
        {params.script} \
        -i {input.bcf} \
        -c {input.txt_chr} \
        -s {input.samples} \
        -m {params.m} \
        -o {output.tsv} \
        > {log} 2>&1
	"""


rule phlash_plot:
    """
    Simple plot of the results of the phlash run, with time is given in years, just provide the generation time.
    """
    # input:
    #     tsv = get_output("phlash_run","_{prefix}.tsv") # From the phlash_run rule
    # output:
    #     png = get_output("phlash_plot","_{prefix}.png")
    # log:
    #     get_log_wild("phlash_run","{prefix}")
    threads: 1
    params:
        r_bin = "Rscript",
        script = get_script("phlash/plot_phlash.R"),
        g = 5, # WARNING: This is the generation time in years, change according to your species.
        s = 100 # Number of generations ago to start the plot from, This is a good idea to set,  since there is little confidence in the estimates for the most recent generations, and it is better to focus on the older generations.
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("r_plots")
    shell:
        """
        {params.r_bin} \
        {params.script} \
        -i {input.tsv} \
        -g {params.g} \
        -s {params.s} \
        -o {output.png} \
        > {log} 2>&1
	"""

rule phlash_plot_pop:
    """
    Similar to the previous one, but with the possibility to plot the results of multiple runs on the same plot, to compare different populations for example.
    This will color according to a greater inclusivity level, in this example runs are populations colored by species.
    This mapping is defined in the input.txt file.
    """
    # input:
    #     tsv = get_output("phlash_run","_{prefix}.tsv") # list of all the tsv files to plot, comma separated
    #     txt = "input/{prefix}.txt" # FORMAT: file with the mapping between population and species, with two columns, the first one is the population name and the second one the species name
    # output:
    #     png = get_output("phlash_plot","_{prefix}.png")
    # log:
    #     get_log_wild("phlash_run","{prefix}")
    threads: 1
    params:
        r_bin = "Rscript",
        script = get_script("phlash/phlash_plot_pop.R"),
        # p = "population_1,population_2,population_3", # Comma separated list of populations to plot, if not provided, all populations will be plotted
        g = 5, # WARNING: This is the generation time in years, change according to your species.
        s = 100, # Number of generations ago to start the plot from
        tsv = lambda w, input: ",".join(input.tsv) # list of all the tsv files to plot, comma separated
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("r_plots")
    shell:
        """
        {params.r_bin} \
        {params.script} \
        -i {params.tsv} \
        -p {params.p} \
        -m {input.txt} \
        -g {params.g} \
        -s {params.s} \
        -o {output.png} \
        > {log} 2>&1
	"""

rule phlash_plot_ccr_three:
    """
    Plots the cross coalescent rate for three population history...
    This rule DOES NOT really plot the CCR, but rather the implementation introduced in the phlash paper (Terhorst, 2025):
    
    n_combined/(n_pop1 + n_pop2)

    Where n represents the instantaneous coalescent rate (eta), and combined is when two populations are considered as one in the phlash_run rule.

    The rate is calculated pairwise, so for the three populations, we will have three pairwise comparisons, and the plot will show the three curves for each pair of populations.
    """
    # input:
    #     tsv_1 = "input/{prefix}_ind1.tsv", # tsv file with the results of the phlash_run for the first population
    #     tsv_2 = "input/{prefix}_ind2.tsv", # tsv
    #     tsv_3 = "input/{prefix}_ind3.tsv", # tsv file with the results of the phlash_run for the third population
    #     tsv_12 = "input/{prefix}_comb12.tsv", # tsv file with the results of the phlash_run for the combination of the first and second population
    #     tsv_13 = "input/{prefix}_comb13.tsv", # tsv
    #     tsv_23 = "input/{prefix}_comb23.tsv", # tsv file with the results of the phlash_run for the combination of the second and third population
    # output:
    #     png = get_output("phlash_plot","_{prefix}.png")
    #     gif = get_output("phlash_plot","_{prefix}.gif") # A gif of the plot, showing the change in the CCR over time
    #     gif3d = get_output("phlash_plot","_{prefix}_3d.gif") # A gif of the plot, showing the change in the CCR over time, but in 3D
    # log:
    #     get_log_wild("phlash_run","{prefix}")
    threads: 1
    params:
        r_bin = "Rscript",
        script = get_script("phlash/plot_ccr_three.R"),
        g = 5, # WARNING: Generation time in years. If 1, then the x-axis will be in generations. Adjust to your species.
        s = 100, # Number of generations ago to start the plot from
        c = 0.95, # confidence interval to plot (quantile of the mcmc samples from phlash if run in different ARG inferences, like with singer)
        # sp_1 = "Species_1", # Species name to be plot
        # sp_2 = "Species_2",
        # sp_3 = "Species_3"
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("r_plots")
    shell:
        """
        {params.r_bin} \
        {params.script} \
        --ind1 {input.tsv_1} \
        --ind2 {input.tsv_2} \
        --ind3 {input.tsv_3} \
        --comb12 {input.tsv_12} \
        --comb13 {input.tsv_13} \
        --comb23 {input.tsv_23} \
        --sp1 {params.sp_1} \
        --sp2 {params.sp_2} \
        --sp3 {params.sp_3} \
        -g {params.g} \
        -s {params.s} \
        -c {params.c} \
        -o {output.png} \
        --gif_output {output.gif} \
        --gif3d_output {output.gif3d} \
        > {log} 2>&1
	"""

rule phlash_get_theta_ne:
    """
    Gets the theta Wattersons estimate and historical Ne from the log of phlash_run
    """
    # input:
    #     rlog = "input/{prefix}.log" # from rules.phlash_run.log # this can be multiple files, it will output everything into one tsv file.
    # output:
    #     tsv = get_output("phlash_get_theta_ne","_{prefix}.tsv")
    # log:
    #     get_log_wild("phlash_get_theta_ne","{prefix}")
    threads: 4
    params:
        r_bin = "Rscript",
        script = get_script("phlash/phlash_get_theta_ne.R"),
        m = 2.04e-9, # WARNING: This is a mutation rate per site per generation, for birds, change according to your species.
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("r_plots")
    shell:
        """
        {params.r_bin} \
        {params.script} \
        -m {params.m} \
        -o {output.tsv} \
        {input.rlog} \
        > {log} 2>&1
        """
