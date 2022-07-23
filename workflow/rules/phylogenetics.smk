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
        msa="results/msa/core_gene_alignment.fa",
        black_list=config["blacklist"]
    output:
        filtered = "results/msa/core_gene_alignment_filtered.fa"
    log:
        "results/logs/blacklist.log"
    threads:
        1
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/filter_samples.py"


# rule tree:
#     """
#     Creates a phylogenetic tree using the augur software of nextstrain. The method can be set in the configuration file.
#     """
#     input: # takes filtered file if blacklist is specified by the user
#         alignment = "results/msa/core_gene_alignment.fa" if config["blacklist"] == '' else "results/msa/core_gene_alignment_filtered.fa"   
#     output:
#         tree = "results/tree/tree.nwk"
#     log:
#         "results/logs/tree/tree.log"
#     threads: 6
#     params:
#         method = config["tree"]["method"],
#         extra = config["tree"]["extra"],
#     conda:
#         "../envs/phylogenetics.yaml"
#     shell:
#         "augur tree --method {params.method} {params.extra} --alignment {input.alignment} --output {output.tree} --nthreads {threads} &> {log}"

rule roary_tree:
    """
    create phylogenetic tree out of core genome alignment using FastTree script, provided by roary
    https://sanger-pathogens.github.io/Roary/
    """
    input:
        alignment =  "results/pangenome/core_gene_alignment.aln"
    output:
        tree = "results/tree/tree.newick"
    log:
        "results/logs/tree/roary_tree.log"
    threads:
        1  # FastTree doesn't seem to have thread specification flag
    params:
        extra = config["tree"]["extra"]
    conda:
        "../envs/pangenome.yaml"
    shell:
        "FastTree -nt -gtr {input.alignment} {params.extra} > {output.tree} 2> {log}"

rule install_plot_roary:
    """
    install script for nice visualization of core genome phlogenetic tree, designed for roary output
    https://github.com/sanger-pathogens/Roary/tree/master/contrib/roary_plots
    """
    output:
        script = "results/installations/roaryplots/roary_plots.py"
    log:
        "results/logs/installation_roaryplots.log"
    threads: 1
    shell:
        "wget https://raw.githubusercontent.com/sanger-pathogens/Roary/master/contrib/roary_plots/roary_plots.py -P results/installations/roaryplots/ 2> {log}"


rule plot_tree:
    """
    Plot tree output of core genome alignmetn by roary, using roary's additional designated script 
    https://github.com/sanger-pathogens/Roary/tree/master/contrib/roary_plots
    """
    input:
        tree = "results/tree/tree.newick",
        gene_presence_csv = "results/pangenome/gene_presence_absence.csv",
        script = "results/installations/roaryplots/roary_plots.py"
    output:
        out_dir = directory("results/pangenome/plots/"),
        freq = "results/pangenome/plots/pangenome_frequency.png",
        mat = "results/pangenome/plots/pangenome_matrix.png", 
        pie = "results/pangenome/plots/pangenome_pie.png"
    conda:  
        "../envs/roaryplot.yaml"
    threads: 1
    params:
        extra = config["plot_roary"]["extra"]
    shell: # run roary_plots.py script and move output from main directory to results
        """
        python {input.script} {input.tree} {input.gene_presence_csv} {params.extra}
        mv pangenome_frequency.png pangenome_matrix.png pangenome_pie.png {output.out_dir}
        """


rule visualize_tree:
    """
    Visualizes the newick format of the phylogenetic tree using a python script. 
    """
    input:
        tree = "results/tree/tree.newick"
    output:
        png = 'results/tree/tree.png'
    conda:
        "../envs/phylogenetics.yaml"
    log:
        "results/logs/tree/visualize_tree.log"
    threads: 1
    script:
        "../scripts/treevisual.py"

