# rule msa:
#     input:
#         sequences = expand("results/fasta/{sample}.fa", sample=IDS) # added "results/" before each path
#     output:
#         alignment = "results/msa/alignment.fasta"
#     log:
#         "results/logs/msa/msa.log"
#     threads: 6
#     conda:
#         "../envs/env.yaml"
#     shell:
#         "python3 -m augur align --sequences {input.sequences} -o {output.alignment} --threads {threads} &> {log}"

rule filter_msa:
    """
    Leaves out samples specified by the user. Runs a designated python script.
    """
    input:
        msa="results/msa/core_gene_alignment.aln",
        black_list=config["blacklist"]
    output:
        filtered = "results/msa/core_gene_alignment_filtered.aln"
    log:
        "results/logs/blacklist.log"
    threads:
        1
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/filter_samples.py"


rule tree:
    """
    Creates a phylogenetic tree using the augur software of nextstrain. The method can be set in the configuration file.
    """
    input: # takes filtered file if blacklist is specified by the user
        alignment = "results/msa/core_gene_alignment.aln" if config["blacklist"] != 'True' else "results/msa/core_gene_alignment_filtered.aln"   
    output:
        tree = "results/tree/tree.nwk"
    log:
        "results/logs/tree/tree.log"
    threads: 6
    params:
        method = config["tree"]["method"],
        extra = config["tree"]["extra"],
    conda:
        "../envs/phylogenetics.yaml"
    shell:
        "augur tree --method {params.method} {params.extra} --alignment {input.alignment} --output {output.tree} --threads {threads} &> {log}"

rule visualize_tree:
    """
    Visualizes the newick format of the phylogenetic tree using a python script. 
    """
    input:
        tree = "results/tree/tree.nwk"
    output:
        png = 'results/tree/tree.png'
    conda:
        "../envs/phylogenetics.yaml"
    threads: 1
    script:
        "../scripts/treevisual.py"

