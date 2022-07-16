"""
Aggregates the txt tables produces by the resistance search using rgi
Outputs a csv file: for each sample, duplicate hits are dropped, keeping the ID with the best hit score. 
Lengthy columns, e.g. with sequence information, are dropped in order to keep the results interpretable
"""

import pandas as pd
import os

resistance_dirs = snakemake.input["resistance"] 
resistance_csv = snakemake.output["csv"] 
eval_threshold = snakemake.params["eval"] 

col_names =  ["sample"] + list((pd.read_csv(resistance_dirs[0], sep="\t", header=0)).columns)
new_df = pd.DataFrame(columns = col_names)

for dir in resistance_dirs:
    sample_name = os.path.basename(dir).split(".")[0] 
    old_df = pd.read_csv(dir, sep="\t", header=0)
    old_df["sample"] = sample_name
    old_df = old_df.sort_values('Best_Hit_Bitscore').drop_duplicates('ID')
    new_df = new_df.append(old_df) 

new_df = new_df.drop(["ORF_ID", "Predicted_DNA","Predicted_Protein", "CARD_Protein_Sequence" ], axis=1) 
new_df.to_csv(resistance_csv, index=False, header=True) 