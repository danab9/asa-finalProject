rule filter_msa:
    """
    Leaves out samples specified by the user. Runs a designated python script.
    """
    input:
        msa="results/pangenome/core_gene_alignment.aln",
        gpa = "results/pangenome/gene_presence_absence.csv",
        black_list=config["blacklist"]
    output:
        filtered_fa = "results/pangenome/core_gene_alignment_filtered.aln",
        filtered_gpa = "results/pangenome/gene_presence_absence_filtered.csv"
    log:
        "results/logs/blacklist.log"
    threads:
        1
    conda:
        "../envs/roaryplot.yaml"
    script:
        "../scripts/filter_samples.py"


rule roary_tree:
    """
    create phylogenetic tree out of core genome alignment using FastTree script, provided by roary
    https://sanger-pathogens.github.io/Roary/
    """
    input:
        alignment =  "results/pangenome/core_gene_alignment.aln" if config["blacklist"] == '' else "results/pangenome/core_gene_alignment_filtered.aln"
    output:
        tree = "results/tree/tree.newick"
    log:
        "results/logs/tree/roary_tree.log"
    threads:
        1  # FastTree doesn't seem to have thread specification option
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
        gene_presence_csv = "results/pangenome/gene_presence_absence.csv" if config["blacklist"] == '' else "results/pangenome/gene_presence_absence_filtered.csv",
        script = "results/installations/roaryplots/roary_plots.py"
    output:
        out_dir = directory("results/pangenome/plots/"),
        freq = "results/pangenome/plots/pangenome_frequency.png",
        mat = "results/pangenome/plots/pangenome_matrix.png", 
        pie = "results/pangenome/plots/pangenome_pie.png"
    conda:  
        "../envs/roaryplot.yaml"
    threads: 1  # there is no multithreads option to roary_plots.py script 
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

