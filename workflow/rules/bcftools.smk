# Rules for bcftools

rule bcftools_filter_samples:
    """
    Generic rule to filter a VCF/BCF file to keep only a subset of samples.
    The samples to keep must be provided in a text file (one sample name per line).
    """
    # We do not define input/output files here, as this is a generic rule.

    # input:
    #    bcf = "path/to/input.bcf", # Can take vcf as input, this is flexible, but try to prefer bcf for storage reasons
    #    samples = "path/to/samples.txt"
    # output:
    #     bcf = "results/genotypes/filtered/{prefix}.bcf",
    #     csi = "results/genotypes/filtered/{prefix}.bcf.csi" # writes csi index by default
    # log:
    #     "logs/bcftools_filter_samples/{prefix}.log"
    threads: 4
    resources:
        mem_mb = lambda wildcards, threads: threads * 1000
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools view \
            --samples-file {input.samples} \
            --output-type b \
            --threads {threads} \
            --write-index=csi \
            --output {output.bcf} \
            {input.bcf} \
        2> {log}
        """

rule bcftools_filter_samples_flexible:
    """
    Filters vcf, accepts the id\t species and comma separated to keep
    """
    # input:
    #     bcf = "path/to/input.bcf",
    #     pop_map = "data/metadata.tsv"
    # output:
    #     bcf = "results/genotypes/filtered/{prefix}.vcf.gz",
    #     csi = "results/genotypes/filtered/{prefix}.vcf.gz.csi"
    # log:
    #     "logs/bcftools_filter/{prefix}.log"
    threads: 4
    params:
        # species_list = "conirostris,magnirostris",
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        REGEX=$(echo "{params.species_list}" | sed 's/,/|/g' | sed 's/^/(/' | sed 's/$/)/')

        awk -v rel="$REGEX" '$2 ~ rel {{print $1}}' {input.pop_map} > {output.bcf}.tmp_list

        bcftools view \
        --threads {threads} \
        --samples-file {output.bcf}.tmp_list \
        -O b \
        -o {output.bcf} \
        --write-index=csi \
        {input.bcf} \
        > {log} 2>&1

        rm {output.bcf}.tmp_list
        """


rule bcftools_filter_informative:
    """
    Recipe: Filter variants to keep only informative ones.
    Default criteria:
      - No singletons (AC>1)
      - Allele frequency (Ref>=2, Alt>=2)
      - Missing data < 20%
    """

    # input:
    #     bcf = "results/genotypes/filtered/{prefix}.bcf"
    # output:
    #     bcf = "results/genotypes/informative/{prefix}.bcf",
    #     csi = "results/genotypes/informative/{prefix}.bcf.csi"
    # log:
    #     "logs/bcftools_filter_informative/{prefix}.log"
    params:
        filter_expr = 'AC>1 && COUNT(GT!="ref")>=2 && COUNT(GT="ref")>=2 && F_MISSING<0.2'
    threads: 4
    resources:
        mem_mb = lambda wildcards, threads: threads * 1000
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools view \
            --include '{params.filter_expr}' \
            --output-type b \
            --threads {threads} \
            --write-index=csi \
            --output {output.bcf} \
            {input.bcf} \
        2> {log}
        """

rule bcftools_stats:
    """
    This rule gets the stats of a vcf file
    """
    # input:
    #     bcf = "input/{prefix}.bcf"
    # output:
    #     txt = get_output("bcftools_stats",".txt")
    # log:
    #     get_log("vcf_stats")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools stats --threads {threads} {input.bcf} > {output.txt}

        PREFIX=$(echo "{output.txt}" | sed 's/\\.[^.]*$//')

        plot-vcfstats -p ${{PREFIX}} {output.txt}
        """


rule bcf2vcf:
    """
    Converts a bcf into vcf.gz format
    """
    # input:
    #     bcf = "input/{prefix}.bcf"
    # output:
    #     vcf = get_output("bcf2vcf","_{sample}.vcf.gz") # never unzip vcfs, not good
    # log:
    #     get_log_wild("bcf2vcf","{sample}")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools view \
        {input.bcf} \
        -O z \
        -o {output.vcf} \
        --threads {threads} \
        2> {log}

        tabix {output.vcf} 2> {log}
        """

rule vcf2bcf:
    """
    Converts a vcf.gz into bcf format
    """
    # input:
    #     vcf = "input/{prefix}.vcf.gz"
    # output:
    #     bcf = get_output("bcf2vcf","_{sample}.bcf") # never unzip vcfs, not good
    # log:
    #     get_log_wild("bcf2vcf","{sample}")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools view \
        {input.vcf} \
        -O b \
        -o {output.bcf} \
        --write-index=csi \
        --threads {threads} \
        2> {log}
        """


rule bcftools_no_chr:
    """
    Removes the 'chr' string from chromosome names
    This is need in some places, like plink
    """
    # input:
    #     vcf = "input/{prefix}.vcf.gz"
    # output:
    #     vcf = get_output("bcftools_no_chr","_{sample}.vcf.gz")
    # log:
    #     get_log_wild("bcftools_no_chr","{sample}")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools \
        view \
        --threads {threads} \
        {input.vcf} \
        | sed 's/^chr//' \
        | bgzip \
        --threads {threads} \
        > {output.vcf} \
        2> {log}
        """


rule bcftools_extract_locus:
    """
    Extract a region of a bcf/vcf file and outputs a new vcf.gz file
    """
    # input:
    #     bcf = "input/{prefix}.vcf.gz"
    # output:
    #     bcf = get_output("bcftools_extract_locus","_{locus}.vcf.gz")
    # log:
    #     get_log_wild("bcftools_extract_locus","{locus}")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    params:
        # chrom = "chr1",
        # start = 0,
        # end = 0
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools \
        view \
        --threads {threads} \
        -r {params.chrom}:{params.start}-{params.end} \
        {input.bcf} \
        -O b \
        -o {output.bcf} \
        > {log} 2>&1

        """

rule tabix:
    """
    Tabindex a vcf.gz file
    """
    # input:
    #     vcf = "input/{prefix}.vcf.gz"
    # output:
    #     tbi = get_output("tabix","_{sample}.vcf.gz.tbi")
    # log:
    #     get_log_wild("tabix","{sample}")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        tabix \
        {input.vcf} \
        2> {log}
        """


rule bcftools_get_chr:
    """
    This rule will split the vcf file into individual vcf files given the {chr} provided
    """
    # input:
    #     # Note that this should have the respective bcf.csi
    #     bcf = "input/{sample}.bcf"
    # output:
    #     bcf = get_output("bcftools_get_chr","_{sample}_{chr}.bcf")
    # log:
    #     get_log("bcftools_get_chr")
    params:
        # region= "{chrom}"
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools view \
        --output-type b \
        --threads {threads} \
        --write-index=csi \
        -r {params.region}\
        --output {output.bcf} \
        {input.bcf} \
        > {log} 2>&1
        """


rule bcftools_fill_tags:
    """
    Fills the tags for the bcf file (needed for some next steps)
    """
    # input:
    #     bcf = "input/{prefix}.bcf"
    # output:
    #     bcf = get_output("bcftools_fill_tags","_{sample}.bcf")
    # log:
    #     get_log_wild("bcftools_fill_tags","{sample}")
    threads: 4
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools +fill-tags \
        --output-type b \
        --threads {threads} \
        --write-index=csi \
        --output {output.bcf} \
        {input.bcf} \
        -- -t all \
        > {log} 2>&1
	"""

rule vcftools_get_allele_freq:
    """
    Gets the allele frequency using vcftools
    """
    # input:
    #     bcf = "input/{prefix}.bcf"
    # output:
    #     tsv = get_output("vcftools_get_allele_freq","_{sample}.tsv")
    # log:
    #     get_log_wild("vcftools_get_allele_freq","{sample}")
    threads: 2
    resources:
        mem_mb=lambda wildcards, threads: threads * 1000
    params:
        o_prefix = lambda w, output: str(output.tsv).replace(".tsv", ""),
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        vcftools \
        --bcf {input.bcf} \
        --counts2 \
        --out {params.o_prefix} \
        --max-alleles 2 \
        > {log} 2>&1

        mv {params.o_prefix}.frq.count {output.tsv}

        sed -i '1s/.*/chrom\tpos\tn_alleles\tn_chr\tref_count\talt_count/' {output.tsv}

        """
