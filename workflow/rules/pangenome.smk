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

