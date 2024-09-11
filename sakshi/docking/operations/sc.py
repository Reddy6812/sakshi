import os
import csv
import yaml
from pathlib import Path

def read_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def extract_name_and_docking_score(config):
    pdbqt_dir = config['operations']['extract_scores']['input_dir']
    output_csv = config['operations']['extract_scores']['output_csv']

    data = []

    for pdbqt_file in Path(pdbqt_dir).glob("*.pdbqt"):
        name = None
        docking_score = None

        with open(pdbqt_file, "r") as file:
            lines = file.readlines()
            for line in lines:
                #print(line[0:50])
                if line.startswith("REMARK VINA RESULT:"):
                    docking_score = float(line.split()[3])
                if line.startswith("REMARK  Name ="):
                    name = line.split('=')[1].strip()
                    break  # Once we have both name and score, we can stop reading this file

        if name and docking_score is not None:
            data.append({'Name': name, 'Docking Score': docking_score})
            print(f"Extracted Name: {name}, Docking Score: {docking_score}")
        else:
            print(f"Failed to extract data from {pdbqt_file.name}")

    # Write to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Name', 'Docking Score']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for row in data:
            writer.writerow(row)

    print(f"Data written to {output_csv}")

if __name__ == "__main__":
    config_path = "configure.yaml"  # Path to your configuration file
    config = read_config(config_path)

    extract_name_and_docking_score(config)
