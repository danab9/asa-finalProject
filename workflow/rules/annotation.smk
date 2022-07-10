# prodigal - https://www.biocode.ltd/catalog-tooldata/prodigal (conda, unsupervised, gff)
# All GFF3 files created by Prokka are valid with Roary and this is the recommended way of generating the input files.
# 
rule prodigal:
    input:
        "results/assembly/{sample}/assembly.fasta"
    output:
        "results/annotations/{sample}_genes.gff"
    log:
        "results/logs/prodigal/{sample}.log"
    conda:
        "../envs/annotation.yaml"
    threads:
        1 # currently prodigal does not handle threads 
    params:
        extra = config["prodigal"]["extra"]
    shell:
        "prodigal  -i {input} -o {output} -f gff {params.extra} &> {log}"