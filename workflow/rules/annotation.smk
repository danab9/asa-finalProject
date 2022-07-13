rule prodigal:
    """
    Generates gene annotation in a .GFF format for each sample using prodigal - https://www.biocode.ltd/catalog-tooldata/prodigal 
    """
    input:
        genome = "results/genomes/{sample}.fasta"
    output:
        annotation = "results/annotations/{sample}_genes.gff"
    log:
        "results/logs/prodigal/{sample}.log"
    conda:
        "../envs/annotation.yaml"
    threads:
        1 # Currently prodigal does not handle threads. 
    params:
        extra = config["prodigal"]["extra"]
    shell:
        "prodigal  -i {input.genome} -o {output.annotation} -f gff {params.extra} &> {log}"
