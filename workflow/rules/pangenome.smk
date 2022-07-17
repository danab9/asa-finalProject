rule pangenome: 
    """
    Creates a core genome using roary https://sanger-pathogens.github.io/Roary/, through a multiple sequence alignment.  
    """
    input:
        #genome = "results/genomes/{sample}.fasta",
        annotation = expand("results/annotations/{sample}_genes.gff", sample=IDS)
    output:
        #dir = directory("results/pangenome"),
        #msa = "results/msa/core_gene_alignment.fa"
        msa = "results/pangenome/core_gene_alignment.aln",
        gene_presence_csv = "results/pangenome/gene_presence_absence.csv"
    log:
        "results/logs/pangenome/roary.log"
    threads: 8
    params:
        extra = config["roary"]["extra"],
        percentage = config["roary"]["percentage_threshold"],
    conda:
        "../envs/pangenome.yaml"
    shell:  # run roary, then extract content from randomally names folder to the results/pangenome/ folder.  
        """
        roary -f results/pangenome/  -e --mafft -p {threads} -cd {params.percentage} {params.extra} {input.annotation} &> {log}
        mv results/pangenome/*/* results/pangenome/
        rmdir results/pangenome/*/
        """  
        # cp results/pangenome/*/core_gene_alignment.aln {output.msa}
#removed {wildcards.sample}  from -o. 
#-o results/pangenome/

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

# rule plot_tree:
#     """
#     Plot tree output of core genome alignmetn by roary, using roary's designated script
#     """
#     input:
#         tree = "results/tree/tree.newick",
#         gene_presence_csv = "results/pangenome/gene_presence_absence.csv"
#     output:
#         touch("treeplot.done")
#     conda:
#         "../envs/pangenome.yaml"
#     shell:
#         "roary_plots.py {input.tree} {input.gene_presence_csv}"
