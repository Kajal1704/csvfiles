import tabula

# Path to your PDF
pdf_path = "mmc2-2.pdf"

# Output CSV file
csv_path = "genotype_data.csv"

# Extract all tables into a DataFrame
tables = tabula.read_pdf(pdf_path, pages="all", multiple_tables=False)

# Save as CSV
tables[0].to_csv(csv_path, index=False)

print(f"✅ CSV saved as {csv_path}")
