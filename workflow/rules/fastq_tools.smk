# PREFIX: fqt

rule fqt_seqkit_stats:
    """
    Stats of  a fastq file using seqkit stats.
    """
    # input:
    #     fastq = "input/{prefix}.fastq"
    # output:
    #     txt = get_output("seqkit_stats","_prefix.txt") # txt file, one path to .fastq per line
    # log:
    #     get_log_wild("seqkit_stats","prefix")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        get_env("fastq_tools")
    shell:
        """
        seqkit stats \
        {input.fastq} \
        -j {threads} \
        > {output.txt} 2>{log}
	"""


