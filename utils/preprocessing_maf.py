import pandas as pd
# import dask.dataframe as dd
import glob
import itertools
from tqdm import tqdm
import numpy as np
from sys import argv
# Directory containing MAF files
maf_dir = f"{argv[1]}/"#"TCGA-CHOL/"  
maf_files = glob.glob(maf_dir + "*/*.maf.*")

# List to store individual MAF dataframes
tumor_list = []
normal_list = []

# Columns of interest
columns_needed = ["Gene", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode", "Variant_Classification", "Reference_Allele", "Tumor_Seq_Allele2", "Match_Norm_Seq_Allele2"]

mutation_types = [
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
    "In_Frame_Ins", "In_Frame_Del", "Splice_Site","Translation_Start_Site", "Nonstop_Mutation",
]

# Read and process MAF files
for i, maf_file in enumerate(tqdm(maf_files)):
    try:
        df = pd.read_csv(maf_file, sep="\t", comment="#", usecols=columns_needed, low_memory=False)
    except Exception as e:
        print(f"{str(e)}: {maf_file}")
    df['Tumor_Sample_Barcode'] = df['Tumor_Sample_Barcode'].apply(lambda x: x[:12])
    df['Matched_Norm_Sample_Barcode'] = df['Matched_Norm_Sample_Barcode'].apply(lambda x: x[:12])
    df = df[df["Variant_Classification"].isin(mutation_types)]
    # Separate tumor and normal sample mutations
    tumor_df = df[["Gene", "Tumor_Sample_Barcode", "Tumor_Seq_Allele2", "Reference_Allele"]].dropna()
    normal_df = df[["Gene", "Matched_Norm_Sample_Barcode", "Match_Norm_Seq_Allele2", "Reference_Allele"]].dropna()

    # Rename columns for consistency
    normal_df.rename(columns={"Matched_Norm_Sample_Barcode": "Tumor_Sample_Barcode"}, inplace=True)

    tumor_list.append(tumor_df)
    normal_list.append(normal_df)

def classify_tumor(row):
    # if pd.notna(row["Tumor_Sample_Barcode"]) and pd.isna(row["Matched_Norm_Sample_Barcode"]):
    #     return "TUMOR_ONLY"
    # elif pd.notna(row["Tumor_Sample_Barcode"]) and pd.notna(row["Matched_Norm_Sample_Barcode"]):
    if row["Tumor_Seq_Allele2"] != row["Reference_Allele"]:
            return 1
    else:
        return 0

def classify_normal(row):
    # if pd.notna(row["Tumor_Sample_Barcode"]) and pd.isna(row["Matched_Norm_Sample_Barcode"]):
    #     return "TUMOR_ONLY"
    # elif pd.notna(row["Tumor_Sample_Barcode"]) and pd.notna(row["Matched_Norm_Sample_Barcode"]):
    if row["Match_Norm_Seq_Allele2"] != row["Reference_Allele"]:
            return 1
    else:
        return 0
    # else:
    #     return "UNKNOWN"
# Combine all tumor and normal data
merged_tumor = pd.concat(tumor_list, ignore_index=True).drop_duplicates()
merged_tumor['Mutation_Status'] = merged_tumor.apply(classify_tumor, axis=1)
merged_tumor = merged_tumor.loc[merged_tumor['Mutation_Status'] != 0]
merged_normal = pd.concat(normal_list, ignore_index=True).drop_duplicates()
merged_normal['Mutation_Status'] = merged_normal.apply(classify_normal, axis=1)
merged_normal = merged_normal.loc[merged_normal['Mutation_Status'] != 0]
# Extract unique genes and patients
tumor_genes = merged_tumor["Gene"].unique()
tumor_patients = merged_tumor["Tumor_Sample_Barcode"].unique()

normal_genes = merged_normal["Gene"].unique()
normal_patients = merged_normal["Tumor_Sample_Barcode"].unique()
# Create binary mutation matrices
combos = pd.MultiIndex.from_product([tumor_genes,tumor_patients, ], names=['Gene','Sample'])
tumor_matrix = pd.DataFrame(index=combos).sort_values(by=['Gene','Sample']).reset_index()
mutation_set = set(zip(merged_tumor['Gene'], merged_tumor['Tumor_Sample_Barcode']))

# Create a new 'Mutation_Status' column using a vectorized approach
tumor_matrix['-1'] = tumor_matrix[['Gene', 'Sample']].apply(
    lambda row: 1 if (row['Gene'], row['Sample']) in mutation_set else 0,
    axis=1
)
tumor_matrix[f"{len(tumor_genes)}"] = (tumor_matrix['Gene'] != tumor_matrix['Gene'].shift()).cumsum()
tumor_matrix[f"{len(tumor_patients)}"] = tumor_matrix.groupby(f"{len(tumor_genes)}").cumcount() + 1
tumor_matrix[f"{len(tumor_genes)}"] -= 1  
tumor_matrix[f"{len(tumor_patients)}"] -= 1  
tumor_matrix = tumor_matrix[[f"{len(tumor_genes)}", f"{len(tumor_patients)}", '-1', 'Gene', 'Sample']]
tumor_matrix.to_csv(f"Tumor_matrix_attempt_2_{maf_dir[:-1]}.txt", sep="\t", index=False)

sample_mapping = tumor_matrix[['Sample', f"{len(tumor_patients)}"]].drop_duplicates()
merged_normal = merged_normal.merge(sample_mapping, how="left", left_on="Tumor_Sample_Barcode", right_on='Sample')

print(merged_normal)
max_sample_number = sample_mapping[f"{len(tumor_patients)}"].max()  # Get last assigned sample number
unmatched_samples = merged_normal[f"{len(tumor_patients)}"].isna()  # Find unmatched rows
merged_normal.loc[unmatched_samples, f"{len(tumor_patients)}"] = range(max_sample_number + 1, max_sample_number + 1 + unmatched_samples.sum())
merged_normal[f"{len(tumor_patients)}"] += 1
merged_normal.sort_values(by=[f"{len(tumor_patients)}"], inplace=True)
merged_normal = merged_normal[['Gene', f"{len(tumor_patients)}"]]
merged_normal.to_csv(f"Normal_list_attempt_2_{maf_dir[:-1]}.txt", sep="\t", index=False, header=False)
# Fill in the matrices

