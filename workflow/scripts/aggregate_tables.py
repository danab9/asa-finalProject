#aggregating tsv tables for plasmids
# read tsv with pandas. 
import pandas as pd
import sys

with open(snakemake.log[0], "w") as log:  # redirect output to log file
    sys.stderr = sys.stdout = log


plasmids_dirs = snakemake.input["plasmids"]
plasmids_csv = snakemake.output["plasmids"]
eval_threshold = snakemake.params["eval"]

# create a new csv, and write a header = ? 
new_df = ? header? 
new_df.to_csv(plasmids_csv, index=False, header=True) 

for resistance_dir in resistance_dirs:
    sample_name = os.path.basename(resistance_dir).split(".")[0] 
    old_df = pd.read_csv(resistance_dir, sep="\t", header=None)
    # make sure that there is only one entry per contig - target match. 
        # filter on E-value, using eval_threshold
            # add to new_df or remove row. 

    

new_df.to_csv(plasmids_csv, mode='a', index=False, header=False) # all append to the same output file? 
