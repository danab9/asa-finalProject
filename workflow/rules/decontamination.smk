rule short_kraken_contaminated: 
    """
    Scanning short contaminated reads for decontamination using Kraken,
    such that the reports can be analysed by the user to provide proper decontamination references. 
    Optional: if screening:short is set to 'True' in the configuration file.  
    """ 
    input:
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1']) 
            if config["trimming"]["short"]=='False' else "results/fastq/trimmed/short/{sample}_1_P.fastq.gz",
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2']) 
            if config["trimming"]["short"]=='False' else "results/fastq/trimmed/short/{sample}_2_P.fastq.gz"
    output:
        clasified_reads_1 = "results/fastq/kraken/contaminated/short/{sample}_1.fq",
        clasified_reads_2 = "results/fastq/kraken/contaminated/short/{sample}_2.fq",
        report = "results/kraken/contaminated/short/{sample}_short_contaminated.txt"
    log:
        "results/logs/kraken/contaminated/short/{sample}.log"
    threads: 10
    params:
        kraken_db = config["kraken"]["database"], 
        quick = "--quick" if config["kraken"]["quick"] == "True" else ""
    conda:
        "../envs/decontamination.yaml"
    shell:
        "kraken2 -db {params.kraken_db} {params.quick} --paired --classified-out results/fastq/kraken/contaminated/short/{wildcards.sample}#.fq {input.r1} {input.r2} --report {output.report} --threads {threads} &> {log}"

rule short_kraken_decontaminated: 
    """
    Scanning short already decontaminated reads for decontamination using Kraken,
    such that the verify that the contamination has properly been removed. 
    Optional: if screening:short is set to 'True' in the configuration file.  
    """ 
    input:
        r1="results/fastq/decontaminated/short/{sample}_1.fq", 
        r2="results/fastq/decontaminated/short/{sample}_2.fq"
    output:
        clasified_reads_1 = "results/fastq/kraken/decontaminated/short/{sample}_1.fq",
        clasified_reads_2 = "results/fastq/kraken/decontaminated/short/{sample}_2.fq",
        report = "results/kraken/decontaminated/short/{sample}_short_decontaminated.txt"
    log:
        "results/logs/kraken/short/decontaminated/{sample}.log"
    threads: 10
    params:
        kraken_db = config["kraken"]["database"], 
        quick = "--quick" if config["kraken"]["quick"] == "True" else ""
    conda:
        "../envs/decontamination.yaml"
    shell:
        "kraken2 -db {params.kraken_db} {params.quick} --paired --classified-out results/fastq/kraken/decontaminated/short/{wildcards.sample}#.fq {input.r1} {input.r2} --report {output.report} --threads {threads} &> {log}"


rule long_kraken_contaminated: 
    """
    Scanning long contaminated reads for decontamination using Kraken,
    such that the reports can be analysed by the user to provide proper decontamination references. 
    Optional: if screening:long is set to 'True' in the configuration file.  
    """ 
    input:
        long=(lambda wildcards: samples.at[wildcards.sample, 'ONT']) 
            if config["trimming"]["long"]=='False' else "results/fastq/trimmed/long/{sample}.fastq.gz",
    output:
        clasified_reads = "results/fastq/kraken/contaminated/long/{sample}.fq",
        report = "results/kraken/contaminated/long/{sample}_long_contaminated.txt"
    log:
        "results/logs/kraken/contaminated/long/{sample}.log"
    threads: 10
    params:
        kraken_db = config["kraken"]["database"], 
        quick = "--quick" if config["kraken"]["quick"] == "True" else ""
    conda:
        "../envs/decontamination.yaml"
    shell:
        "kraken2 -db {params.kraken_db} {params.quick} --classified-out {output.clasified_reads} {input.long} --report {output.report} --threads {threads} &> {log}"

rule long_kraken_decontaminated: 
    """
    Scanning long already decontaminated reads for decontamination using Kraken,
    such that the verify that the contamination has properly been removed. 
    Optional: if screening:long is set to 'True' in the configuration file.  
    """ 
    input:
        long= "results/fastq/decontaminated/long/{sample}.fq",
    output:
        clasified_reads = "results/fastq/kraken/decontaminated/long/{sample}.fq",
        report = "results/kraken/decontaminated/long/{sample}_long_decontaminated.txt"
    log:
        "results/logs/kraken/decontaminated/long/{sample}.log"
    threads: 10
    params:
        kraken_db = config["kraken"]["database"], 
        quick = "--quick" if config["kraken"]["quick"] == "True" else ""
    conda:
        "../envs/decontamination.yaml"
    shell:
        "kraken2 -db {params.kraken_db} {params.quick} --classified-out {output.clasified_reads} {input.long} --report {output.report} --threads {threads} &> {log}"

rule screen: 
    """
    Combines the screening results for short and long, and contaminated and decontaminated reads
    and creates a comprehensive multiqc report.
    Optional: if screening is set to 'True' in the configuration file.  
    """ 
    input:
        screen_short_contaminated = expand("results/kraken/contaminated/short/{sample}_short_contaminated.txt",sample=IDS) 
            if config["screening"]["short"] == "True" else [], 
        screen_long_contaminated = expand("results/kraken/contaminated/long/{sample}_long_contaminated.txt",sample=IDS) 
            if config["screening"]["long"] == "True" else [],
        screen_short_decontaminated = expand("results/kraken/decontaminated/short/{sample}_short_decontaminated.txt",sample=IDS) 
            if config["screening"]["short"] == "True" and config["decontamination"]["short"] == "True" else [], 
        screen_long_decontaminated = expand("results/kraken/decontaminated/long/{sample}_long_decontaminated.txt",sample=IDS) 
            if config["screening"]["long"] == "True" and config["decontamination"]["long"] == "True" else [],
    output:
        multiqc_report = "results/kraken/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "results/logs/multiqc/multiqc_kraken.log"
    threads: 1
    params:
        extra = config["multiqc"]["extra"]  
    shell:
        "multiqc results/kraken {params.extra} -f -o results/kraken &> {log}" #-f is used to overwrite existing reports, to prevent Snakemake from resulting in an error. 

rule short_bowtie2_build_contamination:
    """
    Builds the reference database for the decontamination of the short reads 
    Optional: if decontamination:short is set to 'True' in the configuration file.  
    """ 
    input:
        references = config["contamination_reference"]["short"]
    output:
        multiext(
            "results/references/contamination/contamination_reference_short", 
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",)
    log:
        "results/logs/bowtie2_build/build_contamination.log"
    threads: 4
    conda:
        "../envs/decontamination.yaml"
    shell:
        """
        mkdir -p results/references/contamination 
        bowtie2-build {input.references} results/references/contamination/contamination_reference_short --threads {threads} &> {log}
        """ 
        
rule short_bowtie_map_contaminations:
    """
    Maps the short reads to the reference database for decontamination 
    Optional: if decontamination:short is set to 'True' in the configuration file.  
    """ 
    input:
        indexed_ref= multiext(
            "results/references/contamination/contamination_reference_short", 
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",),
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1'])
                if config["trimming"]["short"]=='False' else "results/fastq/trimmed/short/{sample}_1_P.fastq.gz",
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2'])
                if config["trimming"]["short"]=='False' else "results/fastq/trimmed/short/{sample}_2_P.fastq.gz"
    output:
        "results/sam_contaminations/short/{sample}.sam"
    log:
        "results/logs/bowtie2/contamination_alignment/{sample}.log"
    threads: 6
    conda:
        "../envs/decontamination.yaml"
    shell:
        "bowtie2 -x results/references/contamination/contamination_reference_short -1 {input.r1} -2 {input.r2} -S {output} --threads {threads} &> {log}"

rule short_keep_unmapped:  
    """
    Removes the mapped reads, since these are contaminations. 
    Note: this is a negative filtering approach. 
    For a positive filtering approach the flag can simply be changed. 
    Optional: if decontamination:short is set to 'True' in the configuration file.  
    """ 
    input:
        "results/sam_contaminations/short/{sample}.sam"
    output:
        "results/bam_decontaminated/short/{sample}.bam"
    conda:
        "../envs/decontamination.yaml"
    threads: 4
    log:
        "results/logs/samtools/contaminations/{sample}_bam_unmapped.log"
    shell:
        "samtools view -b -f 12 {input} --threads {threads} > {output} 2> {log}" #The flag -f 12 discards all unmapped mates (paired reads) 

rule short_sam_to_fastq:
    """
    Produces a FASTQ file from the short unpapped reads. 
    Optional: if decontamination:short is set to 'True' in the configuration file.  
    """ 
    input:
         "results/bam_decontaminated/short/{sample}.bam"
    output:
        short_1="results/fastq/decontaminated/short/{sample}_1.fq", 
        short_2="results/fastq/decontaminated/short/{sample}_2.fq"
    conda:
        "../envs/decontamination.yaml"
    log:
        "results/logs/bamtofq/{sample}_decontaminated.log"
    shell:
        "bedtools bamtofastq -i {input} -fq {output.short_1} -fq2 {output.short_2} 2> {log}" 

rule long_decontamination:
    """
    Maps the long reads to the reference database for decontamination 
    Removes the mapped reads, since these are contaminations. 
    Note: this is a negative filtering approach. 
    For a positive filtering approach the flag can simply be changed, from -f to -rf
    Optional: if decontamination:short is set to 'True' in the configuration file.  
    """ 
    input:
        reference = config["contamination_reference"]["long"],
        long = "results/fastq/trimmed/long/{sample}.fastq.gz" if config["trimming"]["long"]=='True' else (lambda wildcards: samples.at[wildcards.sample, 'ONT']),
    output:
        bam = "results/bam_decontaminated/long/{sample}.bam",
        sam ="results/sam_contaminations/long/{sample}.sam",
    log:
        "results/logs/artificialreference/{sample}.log"
    threads: 4
    conda:
        "../envs/decontamination.yaml"
    shell:
        """
        minimap2 -t {threads} -a {input.reference} {input.long} -o {output.sam} 2> {log}  
        samtools view -b -f 4 {output.sam} > {output.bam} 2> {log}
        """  #The flag -f 4 discards all unmapped reads. 


rule long_sam_to_fastq:
    """
    Produces a FASTQ file from the long unpapped reads. 
    Optional: if decontamination:short is set to 'True' in the configuration file.  
    """ 
    input:
        bam = "results/bam_decontaminated/long/{sample}.bam"
    output:
        long="results/fastq/decontaminated/long/{sample}.fq"
    conda:
        "../envs/decontamination.yaml"
    log:
        "results/logs/bamtofq/{sample}_decontaminated.log"
    shell:
        "bedtools bamtofastq -i {input} -fq {output.long} 2> {log}"