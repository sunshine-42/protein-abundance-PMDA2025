import pandas as pd
import numpy as np

#########################
# 1. Filter high-quality data
#########################
def filter_high_quality(df):
    return df[
        (df["is_integrated"] == True) &
        (df["uniprot_status"] == "reviewed") &
        (df["abundance"].notna()) &
        (df["nog_id"].notna())
    ].copy()


#########################
# 2. Log-transform abundance
#########################
def log_transform_abundance(df):
    df["log_abundance"] = np.log1p(df["abundance"])  # log1p for stability (handles 0)
    return df


#########################
# 3. Pivot for intra-species modeling
#########################
def pivot_by_protein(df, organism):
    df_org = df[df["organism_name"] == organism]
    pivot = df_org.pivot_table(index="UniprotAccession", 
                               columns="sample_organ", 
                               values="log_abundance")
    return pivot


#########################
# 4. Pivot for cross-species modeling (via nog_id)
#########################
def pivot_by_nog(df, organism):
    df_org = df[df["organism_name"] == organism]
    pivot = df_org.pivot_table(index="nog_id", 
                               columns="sample_organ", 
                               values="log_abundance")
    return pivot
