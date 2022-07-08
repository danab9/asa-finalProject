rule mlst:
    input:
        genomes = expand("results/genomes/{sample}.fasta", sample=IDS)
    output:
        table = "results/mlst.csv"
    log:
        "results/logs/mlst.log"
    threads: 4
    conda:
        "../envs/mlst.yaml"
    shell:
        "mlst --csv {input.genomes} > {output.table} 2> {log}" 

rule resistance: #per genome; https://github.com/arpcard/rgi#using-rgi-main-genomes-genome-assemblies-metagenomic-contigs-or-proteomes
    input:
        genome = "results/genomes/{sample}.fasta"
    output:
        dir = directory("results/resistance/{sample}"),
        table = "results/resistance/{sample}.tsv",
    log:
        "results/logs/resistance_{sample}.log"
    threads: 10
    conda:
        "../envs/resistance.yaml"
    shell:
        "rgi main --i {input.genome} --output_file {output.dir} --input_type contig --local --clean --threads {threads} 2> {log}"
        

        
        #"rgi {input.genomes} 2> {log}" #genomes/* | 
#  % conda: rgi https://anaconda.org/bioconda/rgi
# https://github.com/arpcard/rgi#using-rgi-main-genomes-genome-assemblies-metagenomic-contigs-or-proteomes
# %Alternatively use blast. 


# rule aggregate_tables:
#     input:
#         mlst = "results/mlst.tsv" if config["mlst"] == "True" else [],
#         resistance = "results/resistance.tsv"