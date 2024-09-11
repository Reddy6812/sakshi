import os
import csv
from pathlib import Path
import logging
import yaml
module_dir = Path(__file__).parent
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)
        
def add_fingerprints_to_csv(config):
    # Load the existing CSV file
    csv_file = config['operations']['extract_scores']['output_csv']
    fingerprints_dir = config['operations']['add_fingerprints']['fingerprints_dir']
    
    rows = []
    with open(csv_file, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            rows.append(row)
    
    # Create a new list to store rows with fingerprints
    updated_rows = []
    
    # Iterate through each row and add the fingerprint data if available
    for row in rows:
        name = row['Name']
        fingerprint_file = Path(fingerprints_dir) / f"{name}.fp"
        
        if fingerprint_file.exists():
            with open(fingerprint_file, 'r') as f:
                fingerprint = f.read().strip()
            row['Fingerprint'] = fingerprint
            logging.info(f"Added fingerprint for {name}")
            updated_rows.append(row)
        else:
            logging.warning(f"No fingerprint file found for {name}. Row will be removed.")
        
    # Write the updated rows back to the same CSV file
    fieldnames = ['Name', 'Docking Score', 'Fingerprint']
    
    with open(csv_file, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in updated_rows:
            writer.writerow(row)
    
    logging.info(f"Updated CSV file saved as {csv_file}")

if __name__ == "__main__":
    config_path = module_dir/'../../configure.yaml'  # Path to your configuration file
    #config_path = config_path.resolve()
    #print(config_path)
    config = read_config(config_path)

    # Add fingerprints to the CSV file
    add_fingerprints_to_csv(config)
