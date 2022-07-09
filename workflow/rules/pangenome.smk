
# % Roary: https://sanger-pathogens.github.io/Roary/ , https://anaconda.org/bioconda/roary
# % only seems to use genes. 
# % `-cd FLOAT percentage of isolates a gene must be in to be core [99]`

# % plus it can make an MAFFT alignment directly:
# % `-n        fast core gene alignment with MAFFT, use with -e 
# % -e        create a multiFASTA alignment of core genes using PRANK `

rule pangenome: 
    input:
        genome = "results/genomes/{sample}.fasta",
        annotation = 
    output:
        dir = directory("results/pangenome/{sample}"),
        table = "results/pangenome/{sample}.txt",
    log:
        "results/logs/pangenome/{sample}.log"
    threads: 8
    params:
        extra = config["roary"]["extra"],
        percentage = config["roary"]["percentage_threshold"]
    conda:
        "../envs/pangenome.yaml"
    shell:
        """roary -e --mafft -p {threads} –f {output.dir} -cd {params.percentage} {params.extra} *.gff 2> {log}"""  