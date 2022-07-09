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

rule resistance: 
    input:
        genome = "results/genomes/{sample}.fasta"
    output:
        dir = directory("results/resistance/{sample}"),
        table = "results/resistance/{sample}.txt",
    log:
        "results/logs/resistance/{sample}.log"
    threads: 8
    params:
        extra = config["rgi"]["extra"],
    conda:
        "../envs/resistance.yaml"
    shell:
        "rgi main -i {input.genome} --output_file {output.dir} {params.extra} --input_type contig -n {threads} 2> {log}" 
#"rgi {input.genomes} 2> {log}" #genomes/* | 
#  % conda: rgi https://anaconda.org/bioconda/rgi
# https://github.com/arpcard/rgi#using-rgi-main-genomes-genome-assemblies-metagenomic-contigs-or-proteomes
# %Alternatively use blast. 

# rule plasmids:
#     input:
#         genome = "results/genomes/{sample}.fasta"
#     output:
#         plasmids = 
# % - plsdb https://ccb-microbe.cs.uni-saarland.de/plsdb
# % conda: https://anaconda.org/ccb-sb/plsdbapi
# % download from: https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/plsdb.fna.bz2

rule download_plasmids_db:
    output:
        database_sequences = "results/plasmids_database/plsdb.fna" 
    log:
        "results/logs/plasmids/download_database_download.log"
    shell: 
        """
        wget -P results/plasmids_database https://ccb-microbe.cs.uni-saarland.de/plsdb/plasmids/download/plsdb.fna.bz2 &> {log}
        bzip2 -dk results/plasmids_database/plsdb.fna.bz2 &> {log}
        """

rule make_plasmids_db:
    input:
        sequences = "results/plasmids_database/plsdb.fna"
    output:
        database = multiext(
            "results/plasmids_database/plsdb.fna",
            ".nhr",
            ".nin",
            ".nog",
            ".nsd",
            ".nsi",
            ".nsq"
        )
    conda:
        "../envs/virulence.yaml"
    log:
        "results/logs/plasmids/download_database_make.log"
    shell:
        """
        makeblastdb -in results/plasmids_database//plsdb.fna -dbtype nucl -parse_seqids -logfile {log}
        """

rule plasmids:
    input:
        genome = "results/genomes/{sample}.fasta",
        database = multiext(
            "results/plasmids_database/plsdb.fna",
            ".nhr",
            ".nin",
            ".nog",
            ".nsd",
            ".nsi",
            ".nsq"
        ) if config["plasmids_database"] == "" else multiext(
            config["plasmids_database"],
            ".nhr",
            ".nin",
            ".nog",
            ".nsd",
            ".nsi",
            ".nsq"
        ),
        database_sequences = "results/plasmids_database/plsdb.fna"
    output:
        table = "results/plasmids/{sample}.tsv",
    log:
        "results/logs/plasmids/{sample}.log"
    threads: 8
    params:
        extra = config["blastn"]["extra"],
    conda:
        "../envs/virulence.yaml"
    shell:
       "blastn -query {input.genome} -db {input.database_sequences} -outfmt 6 -out {output} -num_threads {threads} 2> {log}" 


rule download_virulence_db:
    output:
        database_sequences = "results/virulence_database/VFDB_setA_nt.fas" #config["virulence_database"] #
    log:
        "results/logs/virulence/download_database_download.log"
    shell: 
        """
        wget -P results/virulence_database http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz &> {log}
        gunzip results/virulence_database/VFDB_setA_nt.fas.gz 2> {log}
        """

rule make_virulence_db:
    input:
        sequences = "results/virulence_database/VFDB_setA_nt.fas"
    output:
        database = multiext(
            "results/virulence_database/VFDB_setA_nt.fas",
            ".nhr",
            ".nin",
            ".nog",
            ".nsd",
            ".nsi",
            ".nsq"
        )
    conda:
        "../envs/virulence.yaml"
    log:
        "results/logs/virulence/download_database_make.log"
    shell:
        """
        makeblastdb -in results/virulence_database//VFDB_setA_nt.fas -dbtype nucl -parse_seqids -logfile {log}
        """

rule virulence:
    input:
        genome = "results/genomes/{sample}.fasta",
        database = multiext(
            "results/virulence_database/VFDB_setA_nt.fas",
            ".nhr",
            ".nin",
            ".nog",
            ".nsd",
            ".nsi",
            ".nsq"
        ) if config["virulence_database"] == "" else multiext(
            config["virulence_database"],
            ".nhr",
            ".nin",
            ".nog",
            ".nsd",
            ".nsi",
            ".nsq"
        ),
        database_sequences = "results/virulence_database/VFDB_setA_nt.fas"
    output:
        table = "results/virulence/{sample}.tsv",
    log:
        "results/logs/virulence/{sample}.log"
    threads: 8
    params:
        extra = config["blastn"]["extra"],
    conda:
        "../envs/virulence.yaml"
    shell:
       "blastn -query {input.genome} -db {input.database_sequences} -outfmt 6 -out {output} -num_threads {threads} 2> {log}" 




# blast using wget DNA db from http://www.mgc.ac.cn/VFs/download.htm (why DNA and not protein)
# PathoFact https://git-r3lab.uni.lu/laura.denies/PathoFact (difficult to use)
# multiple fastas at ones, by concatanating them in a multi fasta https://www.researchgate.net/post/How_can_I_create_a_local_BLAST_database_using_multiple_FASTA_files

# rule aggregate_tables:
#     input:
#         mlst = "results/mlst.tsv" if config["mlst"] == "True" else [],
# 1 row per sample 
#         resistance = expand("results/resistance/{sample}.txt",sample=IDS) if config["resistance"] == "True" else [],
# multiple rows (genes/ORFs) per sample  
#         plasmids = ,
#  rows (genes/ORFs) per sample 
#         virulence = expand("results/virulence/{sample}.tsv",sample=IDS) if config["resistance"] == "True" else [], 
# multiple rows different hits) per sample 


