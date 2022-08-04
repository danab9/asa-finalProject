from Bio import Phylo
import matplotlib.pyplot as plt
import sys

with open(snakemake.log[0], "w") as log:  # redirect output to log file
    sys.stderr = sys.stdout = log

    tree_dir = snakemake.input["tree"]
    image_dir = snakemake.output["png"]

    tree = Phylo.read(tree_dir,"newick")
    Phylo.draw(tree)

    plt.savefig(image_dir)