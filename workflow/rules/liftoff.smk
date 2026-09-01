# PREFIX: liftoff 

# Liftoff for lifting over anotations

rule liftoff_run:
    """
    Runs liftoff to lift over annotations from one reference to another.
    """
    # input:
    #     fa_ref = "input/{prefix}.fa",
    #     gff_ref = "input/{prefix}.gff",
    #     fa_target = "input/{prefix}_target.fa"
    # output:
    #     gff = get_output("liftoff_run","_{sample}.gff")
    # log:
    #     get_log_wild("liftoff_run","{sample}")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda: 
        getenv("liftoff")
    shell:
        """
        OUT_ABS=$(readlink -f {output})
        OUT_DIR=$(dirname $OUT_ABS)
        
        liftoff \
        -g {input.gff_ref} \
        -o {output.gff} \
        -dir $OUT_DIR \
        {input.fa_target} \
        {input.fa_ref} \
        -p {threads} \
        > {log} 2>&1
        """


