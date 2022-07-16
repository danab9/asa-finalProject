rule pangenome: 
    """
    Creates a core genome using roary https://sanger-pathogens.github.io/Roary/, through a multiple sequence alignment.  
    """
    input:
        #genome = "results/genomes/{sample}.fasta",
        annotation = expand("results/annotations/{sample}_genes.gff", sample=IDS)
    output:
        #dir = directory("results/pangenome"),
        msa = "results/pangenome/core_gene_alignment.aln",
    log:
        "results/logs/pangenome/roary.log"
    threads: 8
    params:
        extra = config["roary"]["extra"],
        percentage = config["roary"]["percentage_threshold"],
    conda:
        "../envs/pangenome.yaml"
    shell:
        """
        roary -f "results/pangenome/"  -e --mafft -p {threads} -cd {params.percentage} {params.extra} {input.annotation} &> {log}
        cp results/pangenome/*/core_gene_alignment.aln results/pangenome/
        """  
#removed {wildcards.sample}  from -o. 
#-o results/pangenome/