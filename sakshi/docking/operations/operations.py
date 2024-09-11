import os
import csv
import yaml
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel
#from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import pandas as pd
import joblib
import sys
import logging
import subprocess
# Get the current module directory
module_dir = Path(__file__).parent

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def run_script(script_path, args=[]):
    result = subprocess.run(['python3', script_path] + args)
    if result.returncode != 0:
        logging.error(f"Error running {script_path}")
        print(f"Error running {script_path}")
        sys.exit(result.returncode)
    logging.info(f"Script {script_path} ran successfully")

'''
def convert_pdbqt_to_smiles(input_dir, output_dir, output_format="smi"):
    os.makedirs(output_dir, exist_ok=True)
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdbqt", output_format)

    for pdbqt_file in Path(input_dir).glob("*.pdbqt"):
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, str(pdbqt_file))
        smiles = obConversion.WriteString(mol).strip()
        
        output_file = Path(output_dir) / f"{pdbqt_file.stem}.{output_format}"
        #output_file = Path(output_dir) / f"	{pdbqt_file.name}.{output_format}"
        with open(output_file, "w") as f:
            f.write(smiles)
        logging.info(f"Converted {pdbqt_file.name} to {output_format.upper()} format.")
'''

def convert_pdbqt_to_smiles(input_dir, output_file, output_format="smi"):
    # Ensure the output directory exists
    os.makedirs(Path(output_file).parent, exist_ok=True)
    
    # Initialize OpenBabel conversion
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdbqt", output_format)

    # Open the output file in append mode
    with open(output_file, "a") as f:
        # Iterate over all PDBQT files in the input directory
        for pdbqt_file in Path(input_dir).glob("*.pdbqt"):
            mol = openbabel.OBMol()
            obConversion.ReadFile(mol, str(pdbqt_file))
            smiles = obConversion.WriteString(mol).strip()
            
            # Write the SMILES string to the output file with a newline
            f.write(smiles + "\n")
            logging.info(f"Converted {pdbqt_file.name} to {output_format.upper()} format.")
         
'''
dont need this/rough
def convert_pdbqt_to_smiles(input_dir, output_file, output_format="smi"):
    # Ensure the output directory exists
    os.makedirs(Path(output_file).parent, exist_ok=True)
    
    # Initialize OpenBabel conversion
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdbqt", output_format)

    # Open the output file in write mode
    with open(output_file, "w") as f:
        # Iterate over all PDBQT files in the input directory
        for pdbqt_file in Path(input_dir).glob("*.pdbqt"):
            mol = openbabel.OBMol()
            obConversion.ReadFile(mol, str(pdbqt_file))
            smiles = obConversion.WriteString(mol).strip()
            
            # Split the SMILES string and the name
            if smiles:
                smiles_data = smiles.split()
                if len(smiles_data) == 2:
                    smile_str = smiles_data[0]
                    name = smiles_data[1]
                    # Write the SMILES string and the name to the output file
                    f.write(f"{smile_str}\t{name}\n")
                    logging.info(f"Converted {pdbqt_file.name} to {output_format.upper()} format.")
                else:
                    logging.warning(f"SMILES format incorrect for {pdbqt_file.name}: {smiles}")
'''
'''
def generate_fingerprints(smiles_dir):
    fingerprint_dir = Path(smiles_dir).parent / "fingerprints"
    os.makedirs(fingerprint_dir, exist_ok=True)
    
    for smiles_file in Path(smiles_dir).glob("*.smi"):
        with open(smiles_file, "r") as f:
            smiles = f.read().strip()
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logging.warning(f"Failed to parse SMILES: {smiles_file.name}")
            continue
        
        try:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            fp_string = ''.join(map(str, fp.ToBitString()))
            fp_file = fingerprint_dir / f"{smiles_file.stem}.fp"
            with open(fp_file, "w") as f:
                f.write(fp_string)
            logging.info(f"Generated fingerprint for {smiles_file.name}.")
        except Exception as e:
            logging.error(f"Error generating fingerprint for {smiles_file.name}: {e}")
'''
def generate_fingerprints(smiles_file):
    # Set up the directory to save fingerprints
    fingerprint_dir = Path(smiles_file).parent / "../fingerprints1/"
    os.makedirs(fingerprint_dir, exist_ok=True)

    # Open the .smi file and process each line
    with open(smiles_file, "r") as f:
        for line in f:
            try:
                # Split the line into SMILES and name
                smiles, name = line.strip().split('\t')
                # Convert the SMILES to a molecule object
                mol = Chem.MolFromSmiles(smiles)

                if mol is None:
                    logging.warning(f"Failed to parse SMILES for {name}: {smiles}")
                    continue

                # Generate the fingerprint
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
                fp_string = ''.join(map(str, fp.ToBitString()))

                # Write the fingerprint to a file named after the corresponding compound name
                fp_file = fingerprint_dir / f"{name}.fp"
                with open(fp_file, "w") as fp_f:
                    fp_f.write(fp_string)

                logging.info(f"Generated fingerprint for {name} from SMILES.")

            except Exception as e:
                logging.error(f"Error processing line '{line.strip()}': {e}")
                
def extract_name_and_docking_score(config):
    print(module_dir)
    run_script(module_dir / f'sc.py')

    '''
    pdbqt_dir = config['stages']['extract_scores']['input_dir']
    output_csv = config['stages']['extract_scores']['output_csv']
    print(module_dir)
    data = []

    for pdbqt_file in Path(pdbqt_dir).glob("*.pdbqt"):
        name = None
        docking_score = None

        with open(pdbqt_file, "r") as file:
            lines = file.readlines()
            for line in lines:
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
    '''

def add_fingerprints_to_csv(config):
    print(module_dir)
    run_script(module_dir / f'addfing.py')
    
def prepcsv():
    run_script(module_dir / f'prepcsv.py')
    
if __name__ == "__main__":
    config_path = "configure.yaml"  # Path to your configuration file
    config = read_config(config_path)

    # Step 1: Convert PDBQT to SMILES
    
    convert_pdbqt_to_smiles(
        config['operations']['parse_pdbqt']['input_dir'],
        #config['operations']['parse_pdbqt']['output_dir'],
        config['operations']['parse_pdbqt']['output_file'],
        config['operations']['parse_pdbqt']['output_format']
    )
    '''
    
    # Step 2: Generate Fingerprints
    generate_fingerprints(
        config['operations']['parse_pdbqt']['output_dir']
    )
    '''
    
    # Step 2: Generate Fingerprints
    generate_fingerprints(
        config['operations']['parse_pdbqt']['output_file']
    )
    
    # Step 3: Extract docking scores and write to CSV
    extract_name_and_docking_score(config)
    
    #step 4: append fingerprints to csv file
    add_fingerprints_to_csv(config)
    

    #execute prepare csv
    #prepcsv()
    