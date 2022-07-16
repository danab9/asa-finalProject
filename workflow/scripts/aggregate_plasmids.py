"""
Aggregates the tsv tables produces by the plasmids search using blast
Outputs a csv file: for each sample, duplicate hits are dropped, keeping, the one with the lowest E-value and highest percent identity. 
hits that don not fullfill the E-value threshold are also pruned. 
"""

import pandas as pd
import os
import sys

with open(snakemake.log[0], "w") as log:  # redirect output to log file
    sys.stderr = sys.stdout = log

plasmids_dirs = snakemake.input["plasmids"] #["/storage/mi/danab93/asa-finalProject-myrthe/asa-finalProject/results/plasmids/sample1.tsv", "/storage/mi/danab93/asa-finalProject-myrthe/asa-finalProject/results/plasmids/sample2.tsv",] #
plasmids_csv = snakemake.output["csv"] #"test.csv" #
eval_threshold = snakemake.params["eval"] #0.00005 #

col_names =  ["sample", "accession_number", "percent_identity", "length", "drop_column_1", "drop_column_2", "start_position_query", 
"end_position_query", "start_position_hit", "end_position_hit", "eval", "drop_column_3"]
new_df = pd.DataFrame(columns = col_names)

for dir in plasmids_dirs:
    sample_name = os.path.basename(dir).split(".")[0] 
    old_df = pd.read_csv(dir, sep="\t", header=None, names=col_names[1:])
    old_df["sample"] = sample_name
    old_df = old_df.sort_values('percent_identity', ascending=False)
    old_df = old_df.sort_values('eval').drop_duplicates('accession_number')
    new_df = new_df.append(old_df) 
new_df = new_df[new_df['eval'] < eval_threshold] 
new_df = new_df.drop(["drop_column_1", "drop_column_2", "drop_column_3" ], axis=1)
new_df.to_csv(plasmids_csv, index=False, header=True)