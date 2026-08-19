rule seqkit_stats:
    """
    Stats of  a fastq file using seqkit stats.
    """
    # input:
    #     fastq = "input/{prefix}.fastq"
    # output:
    #     txt = get_output("seqkit_stats","_prefix.txt")
    # log:
    #     get_log_wild("seqkit_stats","prefix")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        "envs/fastq_tools.yaml"
    shell:
        """
        seqkit stats \
        {input.fastq} \
        -j {threads} \
        > {output.txt} 2>{log}
	"""


