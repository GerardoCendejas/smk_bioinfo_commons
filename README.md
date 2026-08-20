

# Introduction

This is a repository with all the common functionalities and snakemake rules for my standard bioinformatic pipelines. It is intended to be used as a submodule in other repositories, so that I can easily run common instruction for different projects.


## Architecture

The architecture of this repository and that or the projects using it is as follows:

-   Project
    -   config
        -   config.yaml
    -   workflow
        -   Snakefile
        -   rules
            -   rule1.smk
            -   rule2.smk
        -   scripts
            -   rule1
                -   script1.py
            -   rule2
                -   script2.py
        -   envs
            -   rule1.yaml
            -   rule2.yaml
    -   results
    -   tests
    -   data
    -   results
    -   resources


## Understanding the workflow

The idea of this repo is pretty much inspired by the idea of Biobricks in Synthetic Biology, in which each standard part is meant to be reused and interchangeable with others of the same function/type. Given this, the rules in this repo are meant to be used as standard parts that can be reused in different projects. The rules are designed to be as generic as possible, so that they can be used in different contexts. The rules are also designed to be as modular as possible, so that they can be easily combined with other rules.

The idea is that we can just built our bioinformatic pipelines just by piping (literally just connecting parts) that are compatible and already tested. In this way we can build complex pipelines without having to worry about the details of each part, and we can also easily replace parts with others that are compatible and already tested.

As an example, we can have a general rule to extract samples from a `.bcf` file, and a rule to do any analysis on these type of files, if we have many different rules that do `.bcf` to `.bcf` files, any of these can be use in the pipeline, and this makes easy to change existing files and still be able to pipe them into the next step of the pipeline. This is a very powerful concept, and it is the basis of the idea of Biobricks in Synthetic Biology.


## Basic usage


### Defining output

In order to have standard output creation for tractability, the functions in the [commons module](#orgcc0c2db) are built so you don't have to worry about messy directories, the output that you define in the files will **ALWAYS** be created in the `results` directory, and the log files will be created in the `logs` directory, same for benchmarks. The idea is that you don't have to worry about creating directories, and you can just focus on the analysis.

The usage of this module is pretty simple, in the creation of `input:` and `output:` directives when defining a Snakemake rule, you can use the functions in the [commons module](#orgcc0c2db) to define the proper path, never writing outside of the `results` directory.

The most basic usage would be to use the `get_output` function to define the output file, in this case you need to define the rule name as used in the [Snakefile](#org8b14021) and the [config file](#org059f9e4), as well as the file name per see (don't worry if you don't understand this right now, see the [config file](#org059f9e4) for an explanation on the reason why we work this way). 

    
    get_output("rule_name", "file_name.txt")

For log and benchmark files, you can use the `get_log` and `get_benchmark` functions, respectively, and the file name is not necessary.

1.  Wildcarding

    There is a special case in which you may want to use wildcards for the output and log files, the `get_output` rule works just as well with that, but for the log and benchmark files, you would need to pass the wildcard as a second argument to `get_log` and `get_benchmark`, respectively. In this case you need the filename prefix, with no `.log` extension, as in this example:
    
        
        get_output("rule_name","{wildcard}.txt")
        
        get_log("rule_name", "{wildcard}")
        
        get_benchmark("rule_name", "{wildcard}")
    
    This way is highly encouraged, as it allows for more flexibility in the output files, and it is also more efficient, as it allows for the creation of multiple output files with a single rule, and it will also keep the log and benchmark files organized and one for each rule run.


### Config file

The configuration file present in `config/config.yaml` is the main configuration file for the project, the main purpose of this file is define local variables to use in the Snakefile and specially to define some rule attributes necessary for the proper directory organization of the output files.

The basic structure of a rule description in the configuration file is as follows:

    
    rule_name:
      dir: "output_directory"
      prefix: "descriptor"

There reason for this architecture is that the helper functions in [commons module](#orgcc0c2db) use this configuration tags to name output and log files. merging the prefix and filename with a `_`. Let's look at an example:

    
    get_output("rule_name","test.txt")

`/path/to/project/results/output_directory/descriptor_test.txt`

This organization is very useful, all output will always be in results, in a subdirectory defined by `dir:` and with a descriptor prefix defined by `prefix:`. This makes extremely easy to keep things organized, and your `config.yaml` file will be just a place to tell the rules output where to go. Since the helper functions just considers this things as strings, you can easily define the `dir:` as a path (eg: `dir: "subdir1/subdir2"`), still living inside `results/`.

Other things that would live in your config file would be any variables that you want to use as wildcards or config in your Snakefile.


### Snakefile

TBD: Importing with git submodule to other project.

The basic structure of a rule usage in your `workflow/Snakefile` is as follows:

    
    # Import rule from smk module
    use rule rule_name as local_rule_name with:
        input:
            txt = rules.previous_rule.output.txt
        output:
            csv = get_output("local_rule_name","{wildcard}.csv")
        log:
            get_log("local_rule_name","{wildcard}")
        benchmark:
            get_benchmark("local_rule_name","{wildcard}")
        threads: config["local_rule_name"]["threads"]

The usage of the rule is intended to be local always, so the rule can be called as many times as needed, they are basic building blocks of your pipeline.

Input and output files are always defined with a name (`txt` and `csv` in this example), see specific rules in [modules](#orge1b0b24) for the specific input and output files that are defined for each rule, these names are descriptors of the file type, and are used to define the input and output files in the rule.

They help undestanding also which rules could be used in a pipeline, as the output of one rule could be the input of another rule, and the names of the files help to understand which rules could be used together. All rules have descriptors for input and output, none of them has unnamed input or output files, this is a design decision to make the rules more understandable and easier to use, also facilitates expansion of the rules, as new input and output files can be added without breaking existing rules.

The rule name in the [config file](#org059f9e4) should be the local name of the rule, not the general one, it will not work if you do not define the local name in the config file.


# Modules


## Commons

For a better understanding of this module please first read the [introduction](#org7d12bde) section, it will help you understand the purpose of this module and how to use it.


### get\_output

Helper function to get common output for rules.


### get\_output\_filename

If you have multiple output files in a rule, and for some reason want to use a different prefix instead of making filenames different with the filename, use this helper function instead of [get\_output](#orga1c67f2):

    
    get_output_filename("rule_name","alt_prefix","filename.txt")

What this does is look for the prefix to use in the config file, so it would need to be like this:

    
    rule_name:
      dir: "directory_name"
      prefix: "default_prefix"
      alt_prefix: "alternative_prefix"


### get\_log

Returns the path for logging the rule. Takes an optional `wildcard` argument (file prefix, no extension) to support wildcarding.


### get\_benchmark

Returns the path for benchmarking the rule. Takes an optional `wildcard` argument (file prefix, no extension) to support wildcarding.


### lines2string

Takes a text file and returns a list of each line as a separate string. This is useful when yoy have defined names (such as sample ids) in a file and want to use them as wildcards in the Snakefile.


## [ahmm](workflow/rules/ahmm.smk)

Analyzes the results of a run from [ancestry hmm](https://github.com/russcd/Ancestry_HMM).


## [bcftools](workflow/rules/bcftools.smk)

General rules for working with `.bcf` and `.vcf` files.


## [admixture](workflow/rules/admixture.smk)

General rules for running admixture on a `.bcf`.


## [demesdraw](workflow/rules/demesdraw.smk)

Rules for drawing demographic models in demesformat.


## [dsuite](workflow/rules/dsuite.smk)

Rules for running standard analysis from <https://github.com/millanek/Dsuite.git>. It include analysis like D statistics (ABBA-BABA, fd, df) in genomic windows with dinvestigate. Also performs the f-branch test with fbranch.

Includes functions to plot the f-branch results with the function provided by the dsuite package, and also a rule to plot f-branch values (only between extant populations, no internal branches) in a network graph, as well as a rule to project this same graph onto geographical space using a tsv file determinating the coordinates of each population.


## [genomics\_general](workflow/rules/genomics_general.smk)

Rules for running analysis from [genomics general](https://github.com/simonhmartin/genomics_general.git) package. includes functions to get maximum likelihood phylogenetic trees from `.vcf` files in genomic windows by using phyml.


# File formats

In this architecture, it is important to understand the file formats that are used in the rules, as they are the standard formats for bioinformatic pipelines, as well as some defined formats used for standardization of rules that are described here.

It is important to understand this because the input and output format of the rules also define which rules can be piped together, as the output of one rule can be the input of another rule, or some rules return the same output format as input, and this is important to understand when building pipelines.


## General formats


### log

This is a very important format in this suite. Some rules (usually beacuse of restrictions of programs that they use) cannot either determine the specific files output, or its easier to track the output when we generate a `.log` file.

This is not a log file from snakemake, that will go as normal to `results/logs`. This is only a log file to track files in the pipeline. Every time you see a log flag in the input or output of a rule, check the naming of it (some required a specific name as `success.log`). Some rules use log just to track that the rule ran succesfully and it is not expected to be used afterwards, while some rules that do need a log file as input will tell you from which previous rule log output to use CHECK THE RULE.smk file always to use the proper log file.


### txt

This is a general format of plain text that will be used for very basic configuration or input files, or to output some non standard formatted files from certain programs.


### csv

Comma separated values, common and general format for tabular data, used for input and output of many programs, and also for some intermediate files in the pipelines.


### pop\_map

This is a special format used for population mapping, it is a tab separated file with two columns, the first column is the sample id, and the second column is the population id. This format is used for input of some programs that require population information, and also for some intermediate files in the pipelines.

We use this format as a master helper, this helps keep more readable the function of some rules, and also keeps it easy to understand the input and output of the rules, as well as the pipelines that can be built with them.


### dir

A directory, not a file.


### png

A common format for images, usually means that it will be the final output of a rule.


## Bioinformatic formats


### fastq

The raw reads from NGS.


### vcf

This is the standard format for variant call files, used for input and output of many programs, and also for some intermediate files in the pipelines.


### bcf

This is the binary version of the vcf format, used for input and output of many programs, and also for some intermediate files in the pipelines. It is more efficient than vcf, as it is smaller and faster to read and write. Many rules output this format, preferred to vcf because of the efficiency.


### csi

An index for a `.bcf` files, many rules in this repo already produce a `.csi` index file for the `.bcf` output as common practice, as it is required for many programs that use `.bcf` files as input.


### geno

A format for genotype info, needed for genomics\_general pipelines.


### tre

Newick format for phylogenetic trees.


## Local formats

Some formats named only for certain programs pipelines


### dinv

The result format of d investigate from dsuite package.


### f\_branch

A format needed to run fbranch from dsuited] package.


### f\_branch\_o

Output from fbranch from dsuite package, used for plotting fbranch graphs.

