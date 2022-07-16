"""
script to filter samples prior to phylogeny,
by names appearing in a 'black list', provided by the user
"""
from Bio import SeqIO
import sys

with open(snakemake.log[0], "w") as log:  # redirect output to log file
    sys.stderr = sys.stdout = log


black_list=[line.strip() for line in open(snakemake.input["black_list"])] # read list of unwanted genomes into black list

with open(snakemake.output["filtered"], 'w') as out_file:
    
    for record in SeqIO.parse(snakemake.input["msa"], "fasta"): 
        if record.id.removesuffix('_genes') in black_list:
            continue # skip if id is on black list 
        SeqIO.write(record, out_file, "fasta")  # copy id to new file

