"""
script to filter samples prior to phylogeny,
by names appearing in a 'black list', provided by the user
"""
from Bio import SeqIO
import sys
import pandas as pd  # TO DO: make sure pandas in environment

with open(snakemake.log[0], "w") as log:  # redirect output to log file
    sys.stderr = sys.stdout = log

    # get names and add '_genes' as all headers in roary's output.
    black_list=[line.strip() + '_genes' for line in open(snakemake.input["black_list"])] # read list of unwanted genomes into black list

    # write new fasta msa file.
    with open(snakemake.output["filtered_fa"], 'w') as out_fa:
        
        for record in SeqIO.parse(snakemake.input["msa"], "fasta"): 
            if record.id in black_list:
                continue # skip if id is on black list 
            SeqIO.write(record, out_fa, "fasta")  # copy id to new file


    # write new gene presence and absence table without black list samples
    gpa = pd.read_csv(snakemake.input["gpa"])
    to_keep = [name for name in gpa.columns if name not in black_list]
    gpa_filtered = gpa.filter(to_keep, axis=1)  # filter out black list names

    gpa_filtered.to_csv(snakemake.output["filtered_gpa"], index=False)

