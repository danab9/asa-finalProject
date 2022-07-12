rule short_kraken: 
    """

    """ #TODO: add comments 
    input:
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1']) if config["trimming"]["short"]=='False' else "results/fastq/trimmed/short/{sample}_1_P.fastq.gz",
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2']) if config["trimming"]["short"]=='False' else "results/fastq/trimmed/short/{sample}_2_P.fastq.gz"
    output:
        clasified_reads_1 = "results/fastq/kraken/{sample}_1.fq",
        clasified_reads_2 = "results/fastq/kraken/{sample}_2.fq",
        report = "results/report/kraken/short/{sample}.txt"
    log:
        "results/logs/kraken/short/{sample}.log"
    threads: 10
    params:
        kraken_db = config["kraken"]["database"], 
        quick = "--quick" if config["kraken"]["quick"] == "True" else ""
    conda:
        "../envs/decontamination.yaml"
    shell:
        "kraken2 -db {params.kraken_db} {params.quick} --paired --classified-out results/fastq/kraken/{wildcards.sample}#.fq {input.r1} {input.r2} --report {output.report} --threads {threads} &> {log}"

rule long_kraken: 
    input:
        long=(lambda wildcards: samples.at[wildcards.sample, 'ONT']) if config["trimming"]["long"]=='False' else "results/fastq/trimmed/short/{sample}_long.fastq.gz",
    output:
        clasified_reads = "results/fastq/kraken/{sample}_long.fq",
        report = "results/report/kraken/long/{sample}.txt"
    log:
        "results/logs/kraken/long/{sample}.log"
    threads: 10
    params:
        kraken_db = config["kraken"]["database"], 
        quick = "--quick" if config["kraken"]["quick"] == "True" else ""
    conda:
        "../envs/decontamination.yaml"
    shell:
        "kraken2 -db {params.kraken_db} {params.quick} --classified-out {output.clasified_reads} {input.long} --report {output.report} --threads {threads} &> {log}"

rule screen: # for short and long, #creates a multiqc report of the screened reads in the folder qc/screened
    input:
        screen_short = expand("results/report/kraken/short/{sample}.txt",sample=IDS) if config["screening"]["short"] == "True" else [], #for all samples. 
        screen_long = expand("results/report/kraken/long/{sample}.txt",sample=IDS) if config["screening"]["long"] == "True" else [],
    output:
        "results/report/kraken/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    log:
        "results/logs/multiqc/multiqc_kraken.log"
    threads: 1
    params:
        config["multiqc"]["extra"]  
    shell:
        "multiqc results/report/kraken {params} -o results/report/kraken &> {log}"

rule short_bowtie2_build_contamination:
    input:
        config["contamination_reference"]["short"]
    output:
        multiext(
            "results/references/contamination/contamination_reference_short", 
            ".fa",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    log:
        "results/logs/bowtie2_build/build_contamination.log"
    threads: 4
    conda:
        "../envs/decontamination.yaml"
    shell:
        """
        mkdir -p results/references/contamination && cp {input} results/references/contamination/contamination_reference_short.fa
        bowtie2-build {input} results/references/contamination/contamination_reference_short --threads {threads} &> {log}
        """
    # Error No output file specified!, to fix this I added {input} for a 2nd time

rule short_bowtie_map_contaminations:
    input:
        ref="results/references/contamination/contamination_reference_short.fa",
        #@Sandro: It is true that we don't use the config defined directory here, but this dir is used in the build rule above,
        # this rule copies the user defined file name to this new location, so that bowtie build can use it.
        indexed_ref= multiext(
            "results/references/contamination/contamination_reference_short",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2"),
        r1=(lambda wildcards: samples.at[wildcards.sample, 'fq1'])
                if config["trimming"]["short"]=='False' else "results/fastq/trimmed/short/{sample}_1_P.fastq.gz",
        r2=(lambda wildcards: samples.at[wildcards.sample, 'fq2'])
                if config["trimming"]["short"]=='False' else "results/fastq/trimmed/short/{sample}_2_P.fastq.gz"
    output:
        "results/sam_contaminations/{sample}.sam"
    log:
        "results/logs/bowtie2/contamination_alignment/{sample}.log"
    threads: 6
    conda:
        "../envs/decontamination.yaml"
    shell:
        "bowtie2 -x results/references/contamination/contamination_reference_short -1 {input.r1} -2 {input.r2} -S {output} --threads {threads} &> {log}"

rule short_keep_unmapped:   # TODO: see https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd
    input:
        "results/sam_contaminations/{sample}.sam"
    output:
        "results/bam_decontaminated/{sample}.bam"
    conda:
        "../envs/decontamination.yaml"
    threads: 4
    log:
        "results/logs/samtools/contaminations/{sample}_bam_unmapped.log"
    shell:
        "samtools view -b -f 4 {input} --threads {threads} > {output} 2> {log}"

rule short_sam_to_fastq:
    input:
         "results/bam_decontaminated/{sample}.bam"
    output:
        short_1="results/fastq/decontaminated/short/{sample}_1.fq", 
        short_2="results/fastq/decontaminated/short/{sample}_2.fq"
    conda:
        "../envs/decontamination.yaml"
    log:
        "results/logs/bamtofq/{sample}_decontaminated.log"
    shell:
        "bedtools bamtofastq -i {input} -fq {output.short_1} -fq2 {output.short_2} 2> {log}"

rule long_sam_to_fastq:
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

rule long_decontamination:
    input:
        reference = config["contamination_reference"]["long"],
        long = "results/fastq/trimmed/long/{sample}_1_P.fastq.gz" if config["trimming"]["long"]=='True' else (lambda wildcards: samples.at[wildcards.sample, 'ONT']),
    output:
        bam = "results/bam_decontaminated/long/{sample}.bam"
    log:
        "results/logs/artificialreference/{sample}.log"
    threads: 4
    conda:
        "../envs/decontamination.yaml"
    shell:
        "minimap2 -t {threads} -a {input.reference} {input.long} 2> {log} | samtools view -b - > {output.bam} 2> {log}"
