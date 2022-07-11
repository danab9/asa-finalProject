rule pangenome: 
    """
    Creates a core genome using roary https://sanger-pathogens.github.io/Roary/, through a multiple sequence alignment.  
    """
    input:
        #genome = "results/genomes/{sample}.fasta",
        annotation = expand("results/annotations/{sample}_genes.gff", sample=IDS)
    output:
        #dir = directory("results/pangenome"),
        table = "results/pangenome/{sample}.txt",
    log:
        "results/logs/pangenome/{sample}.log"
    threads: 8
    params:
        extra = config["roary"]["extra"],
        percentage = config["roary"]["percentage_threshold"],
    conda:
        "../envs/pangenome.yaml"
    shell:
        """roary -f "results/pangenome/" -o results/pangenome/{wildcards.sample} -e --mafft -p {threads} -cd {params.percentage} {params.extra} {input.annotation} 2> {log}"""  
