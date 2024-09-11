import os
from openbabel import pybel
from pathlib import Path
import sys
import logging

# Set up logging
log_file = Path(__file__).parent / "sdf_pdbqt_monitor.log"
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(message)s')

def convert_sdf_to_pdbqt(sdf_file, output_dir):
    """
    Convert SDF file to PDBQT format.

    Parameters:
    - sdf_file (str): Path to the input SDF file.
    - output_dir (str): Directory to save the converted PDBQT files.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Open the SDF file using Pybel
    mols = pybel.readfile("sdf", sdf_file)
    
    for mol in mols:
        # Add hydrogen atoms if needed
        mol.addh()

        # Set the output file path
        output_file = Path(output_dir) / f"{mol.title}.pdbqt"
        
        # Convert to PDBQT format
        mol.write("pdbqt", str(output_file), overwrite=True)
        logging.info(f"Converted {mol.title} to PDBQT format.")
        print(f"Converted {mol.title} to PDBQT format.")

if __name__ == "__main__":
    # Check if the script is called with the necessary arguments
    if len(sys.argv) < 3:
        logging.error("Not enough arguments provided. Usage: python sdf_pdbqt.py <sdf_file> <output_dir>")
        print("Usage: python sdf_pdbqt.py <sdf_file> <output_dir>")
        sys.exit(1)

    # Read arguments
    sdf_file = sys.argv[1]
    output_dir = sys.argv[2]

    # Convert the SDF to PDBQT
    convert_sdf_to_pdbqt(sdf_file, output_dir)
