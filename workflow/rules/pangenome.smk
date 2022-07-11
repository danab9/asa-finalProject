
# % Roary: https://sanger-pathogens.github.io/Roary/ , https://anaconda.org/bioconda/roary
# % only seems to use genes. 
# % `-cd FLOAT percentage of isolates a gene must be in to be core [99]`

# % plus it can make an MAFFT alignment directly:
# % `-n        fast core gene alignment with MAFFT, use with -e 
# % -e        create a multiFASTA alignment of core genes using PRANK `
# sample gff https://github.com/AlgoLab/MALVIRUS-tutorial-data/blob/master/sars-cov-2.gff

rule pangenome:  
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
