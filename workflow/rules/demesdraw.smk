# Rules to plot demographic models using the `demesdraw` library in Python.

rule demesdraw_plot_model:
    """
    Plots a demographic model from a .yaml file
    """
    # input:
    #     demes = "input/{prefix}.yaml" # a demes formated yaml file
    # output:
    #     png_1 = get_output("demesdraw_plot_model","_{prefix}_tubes.png"),
    #     png_2 = get_output("demesdraw_plot_model","_{prefix}_size.png")
    # log:
    #     get_log_wild("demesdraw_plot_model","{prefix}")
    threads: 1
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    params:
        python_bin = "python",
        script = "workflow/scripts/demesdraw/plot_model.py",
    conda:
        "envs/demesdraw.yaml"
    shell:
        """
        {params.python_bin} \
        {params.script} \
        -i {input.demes} \
        -t {output.png_1} \
        -s {output.png_2} \
        > {log} 2>&1
	"""


