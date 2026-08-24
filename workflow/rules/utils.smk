# PREFIX: utils

# Utils tools for Snakemake workflows
# File processing and helper functions

rule utils_gunzip:
    """
    Generic rule to decompress any .gz file.
    Usage: Request the file without .gz extension.
    """
    # input:
    #     gz = "{filepath}.gz"
    # output:
    #     uzip = "{filepath}"
    # log:
    #     "logs/gunzip/{filepath}.log"
    wildcard_constraints:
        filepath=r".*(?<!\.gz)"
    shell:
        # -c: Writes to stdout, so we can redirect to a file
        # Keep the original file, so we can use it again if needed
        "gunzip -c {input.gz} > {output.uzip} 2> {log}"

rule utils_touch:
    """
    touches a file
    """
    # output:
    #     txt = get_output("touch","_prefix.txt")
    # log:
    #     get_log_wild("touch","prefix")
    threads: 1
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    shell:
        """
        touch {output.txt} > {log} 2>&1
	"""
rule utils_csv2tsv:
    """
    Transforms a csv file into a tsv file.
    """
    # input:
    #     csv = "input/{prefix}.csv"
    # output:
    #     tsv = get_output("utils_csv2tsv","_prefix.tsv")
    # log:
    #     get_log_wild("utils_csv2tsv","prefix")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    shell:
        """
        awk -v FS="," -v OFS="\t" '{ $1 = $1; print $0 }' \
        {input.csv} > {output.tsv} \
        2> {log}
	"""



        


        
