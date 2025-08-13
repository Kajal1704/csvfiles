import pandas as pd
from io import StringIO

soft_file = "GDS4844_full.soft"  # Change to your file name

with open(soft_file, "r") as f:
    lines = f.readlines()

# Locate the start and end of the dataset table
start_idx = next(i for i, line in enumerate(lines) if line.startswith("!dataset_table_begin")) + 1
end_idx = next(i for i, line in enumerate(lines) if line.startswith("!dataset_table_end"))

# Extract and load into DataFrame
table_data = "".join(lines[start_idx:end_idx])
df = pd.read_csv(StringIO(table_data), sep="\t")

# Save to CSV
output_file = "GDS4844.csv"
df.to_csv(output_file, index=False)

print(f"✅ Conversion complete. CSV saved as: {output_file}")
print("🔍 Preview:")
print(df.head())
