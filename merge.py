import pandas as pd
import glob
import os

# Path where your CSV files are stored
path = "csvfiles"  # adjust if needed

# Get all CSV files inside the folder
all_files = glob.glob(os.path.join(path, "*.csv"))

# Read and merge all CSVs
dfs = []
for filename in all_files:
    df = pd.read_csv(filename)
    dfs.append(df)

# Combine all into one dataset
merged_df = pd.concat(dfs, ignore_index=True)

# Save the merged dataset
merged_df.to_csv("combined_dataset.csv", index=False)

print(f"Merged {len(all_files)} files. Output saved as combined_dataset.csv")
