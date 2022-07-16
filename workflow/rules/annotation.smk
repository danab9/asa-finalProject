# rule prodigal:
#     """
#     Generates gene annotation in a .GFF format for each sample using prodigal - https://www.biocode.ltd/catalog-tooldata/prodigal 
#     """
#     input:
#         genome = "results/genomes/{sample}.fasta"
#     output:
#         annotation = "results/annotations/{sample}_genes.gff"
#     log:
#         "results/logs/prodigal/{sample}.log"
#     conda:
#         "../envs/annotation.yaml"
#     threads:
#         1 # Currently prodigal does not handle threads. 
#     params:
#         extra = config["prodigal"]["extra"]
#     shell:
#         "prodigal  -i {input.genome} -o {output.annotation} -f gff {params.extra} &> {log}"


rule prokka:
    """
    Generates gene annotation in a .GFF format for each sample using prokka - https://github.com/tseemann/prokka
    """
    input:
        genome = "results/genomes/{sample}.fasta"
    output:
        annotation = "results/annotations/{sample}_genes.gff"
    log:
        "results/logs/prokka/{sample}.log"
    conda:
        "../envs/prokka.yaml"
    threads:
        8 
    params:
        extra = config["prokka"]["extra"],
        dir = directory("results/annotations/")
    shell:  # using '--force' because snakemkae automatically creates the output directory prior to the run, which gives an error 
        "prokka --outdir {params.dir} --prefix {wildcards.sample}_genes --force {params.extra} {input.genome} &> {log}"
