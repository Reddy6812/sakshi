import os
import csv
import yaml
from pathlib import Path
import sys
import logging
import subprocess
import pandas as pd
# Get the current module directory
module_dir = Path(__file__).parent

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)
        
def process_molecule_variants(config):
    csv_file = config['operations']['extract_scores']['output_csv']
    df = pd.read_csv(csv_file)

    # Step 1: Split the 'Name' column to separate Mol_ID and Variant
    split_columns = df['Name'].str.rsplit('-', n=1, expand=True)

    # Check if the split resulted in two columns or not
    if len(split_columns.columns) == 1:
        # If there's no underscore, assign a default variant value of '1'
        df['Mol_ID'] = split_columns[0]
        df['Variant'] = 1
    else:
        df['Mol_ID'] = split_columns[0]
        df['Variant'] = split_columns[1].astype(int)

    # Step 2: Sort by Mol_ID, then by Docking Score (ascending), and Variant (ascending)
    df = df.sort_values(by=['Mol_ID', 'Docking Score', 'Variant'], ascending=[True, True, True])

    # Step 3: Identify and track duplicates to be deleted
    duplicates_to_delete = df[df.duplicated(subset=['Mol_ID'], keep='first')]

    # Print the Mol_IDs and their variants that will be deleted
    if not duplicates_to_delete.empty:
        print("Deleted Mol_IDs and their variants:")
        for _, row in duplicates_to_delete.iterrows():
            print(f"Mol_ID: {row['Mol_ID']}, Variant: {row['Variant']}, Docking Score: {row['Docking Score']}")
    else:
        print("No duplicates found to delete.")

    # Step 4: Remove duplicates, keeping the row with the smallest Docking Score
    df = df.drop_duplicates(subset=['Mol_ID'], keep='first')

    # Step 5: Combine Mol_ID and Variant back into the 'Name' column
    df['Name'] = df['Mol_ID'] + '-' + df['Variant'].astype(str)

    # Drop the temporary columns
    df = df.drop(columns=['Mol_ID', 'Variant'])

    # Write the processed data back to the CSV
    df.to_csv(csv_file, index=False)
    logging.info(f"Processed molecule variants and updated {csv_file}")



# Step 3: Extract docking scores and write to CSV
#extract_name_and_docking_score(config)
if __name__ == "__main__":
    config_path = "configure.yaml"  # Path to your configuration file
    config = read_config(config_path)
    # Step 4: Process molecule variants
    process_molecule_variants(config)
