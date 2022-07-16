rule mlst:
    """
    Multi-locus sequence typing using the mlst tool, ran for all samples at once. 
    Optional: if mlst is set to 'True' in the configuration file. 
    """ 
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


# rule download_resistance_db:
#     """
#     Downloads a plasmid reference database from https://ccb-microbe.cs.uni-saarland.de/plsdb
#     in case no database was provided by the user under "resistance_database" in the configuration file. 
#     """
#     output:
#         database_sequences = "results/resistance_database/nucleotide_fasta_protein_homolog_model.fasta"
#     log:
#         "results/logs/resistance/download_database_download.log"
#     shell: 
#         """
#         wget -P results/resistance_database https://card.mcmaster.ca/latest/data &> {log}
#         mv results/resistance_database/data results/resistance_database/card-data.tar.bz2
#         bzip2 -dk results/resistance_database/card-data.tar.bz2 &> {log}
#         tar -xvf results/resistance_database/card-data.tar -C results/resistance_database &> {log}
#         """

# rule make_resistance_db:
#     """
#     Makes the downloaded plamid database ready for search with blast 
#     in case no database was provided by the user under "resistance_database" in the configuration file. 
#     """
#     input:
#         sequences = "results/resistance_database/nucleotide_fasta_protein_homolog_model.fasta"
#     output:
#         database = multiext(
#             "results/resistance_database/nucleotide_fasta_protein_homolog_model.fasta",
#             ".nhr",".nin",".nog",".nsd",".nsi",".nsq")
#     conda:
#         "../envs/virulence.yaml"
#     log:
#         "results/logs/resistance/download_database_make.log"
#     shell:
#         "makeblastdb -in results/resistance_database/nucleotide_fasta_protein_homolog_model.fasta -dbtype nucl -parse_seqids -logfile {log}"

# rule resistance:
#     """
#     Screens for resistance genes in each sample using the PLSD database and blast. 
#     Optional: if resistance is set to 'True' in the configuration file. 
#     """ 
#     input:
#         genome = "results/genomes/{sample}.fasta",
#         database = multiext(
#             "results/resistance_database/nucleotide_fasta_protein_homolog_model.fasta",
#             ".nhr",".nin",".nog",".nsd",".nsi",".nsq"
#         ) if config["resistance_database"] == "" else multiext(
#             config["resistance_database"],
#             ".nhr",".nin",".nog",".nsd",".nsi",".nsq"),
#         database_sequences = "results/resistance_database/nucleotide_fasta_protein_homolog_model.fasta"
#     output:
#         table = "results/resistance/{sample}.tsv",
#     log:
#         "results/logs/resistance/{sample}.log"
#     threads: 8
#     params:
#         extra = config["blastn"]["extra"],
#     conda:
#         "../envs/virulence.yaml"
#     shell:
#        "blastn -query {input.genome} -db {input.database_sequences} -outfmt 6 -out {output} -num_threads {threads} 2> {log}" 

rule resistance: 
    """
    Screens for antibiotic resistance genes in each sample using the CARD database 
    and the rgi tool, https://github.com/arpcard/rgi#using-rgi-main-genomes-genome-assemblies-metagenomic-contigs-or-proteomes
    Optional: if resistance is set to 'True' in the configuration file. 
    """ 
    input:
        genome = "results/genomes/{sample}.fasta"
    output:
        table = "results/resistance/{sample}.txt" #"results/resistance/{sample}.json",
    log:
        "results/logs/resistance/{sample}.log"
    threads: 8
    params:
        extra = config["rgi"]["extra"],
        dir = "results/resistance/{sample}",
    conda:
        "../envs/resistance.yaml"
    shell:
        """
        rgi main -i {input.genome} --output_file {output.dir} {params.extra} --input_type contig -n {threads} 2> {log}
        #rgi heatmap --input ../storage/mi/danab93/asa-finalProject-myrthe/asa-finalProject/results/resistance
        """ 
# TODO: check with 1 core. 

rule download_plasmids_db:
    """
    Downloads a plasmid reference database from https://ccb-microbe.cs.uni-saarland.de/plsdb
    in case no database was provided by the user under "plasmids_database" in the configuration file. 
    """
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
    """
    Makes the downloaded plamid database ready for search with blast 
    in case no database was provided by the user under "plasmids_database" in the configuration file. 
    """
    input:
        sequences = "results/plasmids_database/plsdb.fna"
    output:
        database = multiext(
            "results/plasmids_database/plsdb.fna",
            ".nhr",".nin",".nog",".nsd",".nsi",".nsq")
    conda:
        "../envs/virulence.yaml"
    log:
        "results/logs/plasmids/download_database_make.log"
    shell:
        "makeblastdb -in results/plasmids_database//plsdb.fna -dbtype nucl -parse_seqids -logfile {log}"

rule plasmids:
    """
    Screens for plasmids genes in each sample using the PLSD database and blast. 
    Optional: if plasmids is set to 'True' in the configuration file. 
    """ 
    input:
        genome = "results/genomes/{sample}.fasta",
        database = multiext(
            "results/plasmids_database/plsdb.fna",
            ".nhr",".nin",".nog",".nsd",".nsi",".nsq"
        ) if config["plasmids_database"] == "" else multiext(
            config["plasmids_database"],
            ".nhr",".nin",".nog",".nsd",".nsi",".nsq"),
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
    """
    Downloads the virulence factors (VF) DNA core database from http://www.mgc.ac.cn/VFs/download.htm
    in case no database was provided by the user under "virulence_database" in the configuration file. 
    """
    output:
        database_sequences = "results/virulence_database/VFDB_setA_nt.fas" 
    log:
        "results/logs/virulence/download_database_download.log"
    shell: 
        """
        wget -P results/virulence_database http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz &> {log}
        gunzip results/virulence_database/VFDB_setA_nt.fas.gz 2> {log}
        """

rule make_virulence_db:
    """
    Makes the downloaded VF database ready for search with blast 
    in case no database was provided by the user under "virulence_database" in the configuration file. 
    """
    input:
        sequences = "results/virulence_database/VFDB_setA_nt.fas"
    output:
        database = multiext(
            "results/virulence_database/VFDB_setA_nt.fas",
            ".nhr",".nin",".nog",".nsd",".nsi",".nsq")
    conda:
        "../envs/virulence.yaml"
    log:
        "results/logs/virulence/download_database_make.log"
    shell:
        """
        makeblastdb -in results/virulence_database//VFDB_setA_nt.fas -dbtype nucl -parse_seqids -logfile {log}
        """

rule virulence:
    """
    Screens for virulence factors in each sample using the VF database and blast. 
    Optional: if virulence is set to 'True' in the configuration file. 
    """ 
    input:
        genome = "results/genomes/{sample}.fasta",
        database = multiext(
            "results/virulence_database/VFDB_setA_nt.fas",
            ".nhr",".nin",".nog",".nsd",".nsi",".nsq"
            ) if config["virulence_database"] == "" else multiext(
            config["virulence_database"],
            ".nhr",".nin",".nog",".nsd",".nsi",".nsq"),
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


rule aggregate_plasmids:
    """
    Aggregates the results from the blast plasmids search. 
    """ 
    input:
        plasmids = expand("results/plasmids/{sample}.tsv",sample=IDS),
    params:
        eval = config["blastn"]["eval"],
    log:
        "results/logs/plasmids/aggregate.log"
    output: 
        csv = "results/plasmids.csv"
    script:
        "../scripts/aggregate_plasmids.py"


rule aggregate_virulence:
    """
    Aggregates the results from from the blast virulence search. 
    """ 
    input:
        virulence = expand("results/virulence/{sample}.tsv",sample=IDS),
    params:
        eval = config["blastn"]["eval"],
    log:
        "results/logs/virulence/aggregate.log"
    output: 
        csv = "results/virulence.csv"
    script:
        "../scripts/aggregate_virulence.py"


rule aggregate_resistance:
    """
    Aggregates the results from the rgi resistance search. 
    """ 
    input:
        resistance = expand("results/resistance/{sample}.txt",sample=IDS),
    params:
        eval = config["blastn"]["eval"],
    log:
        "results/logs/resistance/aggregate.log"
    output: 
        csv = "results/resistance.csv"
    script:
        "../scripts/aggregate_resistance.py"

rule downstream_analyses:
    """
    Triggers the various analyses steps, including MLST, plasmid, virulence factor and antbiotic resistance screening. 
    """ 
    input:
        mlst = "results/mlst.csv" if config["mlst"] == "True" else [],
        resistance = "results/resistance.csv" if config["resistance"] == "True" else [],
        plasmids = "results/plasmids.csv" if config["plasmids"] == "True" else [],
        virulence = "results/virulence.csv" if config["virulence"] == "True" else [], 
    output: 
        touch_file = touch(temp("results/logs/downstream_analyses.done"))


