import os
from openbabel import openbabel
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def convert_pdb_to_pdbqt(input_dir, output_dir):
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize OpenBabel conversion
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdbqt")

    # Iterate over all PDB files in the input directory
    for pdb_file in Path(input_dir).glob("*.pdb"):
        mol = openbabel.OBMol()
        
        # Read the PDB file
        obConversion.ReadFile(mol, str(pdb_file))
        
        # Prepare output file path
        output_file = Path(output_dir) / f"{pdb_file.stem}.pdbqt"
        
        # Write the molecule to PDBQT format
        obConversion.WriteFile(mol, str(output_file))
        
        logging.info(f"Converted {pdb_file.name} to PDBQT format.")

# Example usage
input_dir = "/blue/yanjun.li/vi.gade1/seabra-li/data/pdb_files/"  # Update this path to your input directory
output_dir = "/blue/yanjun.li/vi.gade1/seabra-li/data/pdb_files/"  # Update this path to your output directory
convert_pdb_to_pdbqt(input_dir, output_dir)
