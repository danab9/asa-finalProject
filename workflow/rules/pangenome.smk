
# % Roary: https://sanger-pathogens.github.io/Roary/ , https://anaconda.org/bioconda/roary
# % only seems to use genes. 
# % `-cd FLOAT percentage of isolates a gene must be in to be core [99]`

# % plus it can make an MAFFT alignment directly:
# % `-n        fast core gene alignment with MAFFT, use with -e 
# % -e        create a multiFASTA alignment of core genes using PRANK `
# sample gff https://github.com/AlgoLab/MALVIRUS-tutorial-data/blob/master/sars-cov-2.gff

rule pangenome:  
    input:
        genome = "results/genomes/{sample}.fasta",
        annotation = "results/annotations/{sample}_genes.gff"
    output:
        dir = directory("results/pangenome"),
        table = "results/pangenome.txt",
    log:
        "results/logs/pangenome.log"
    threads: 8
    params:
        extra = config["roary"]["extra"],
        percentage = config["roary"]["percentage_threshold"],
    conda:
        "../envs/pangenome.yaml"
    shell:
        """roary -e --mafft -p {threads} -cd {params.percentage} {params.extra} {input.annotation} 2> {log}"""  
