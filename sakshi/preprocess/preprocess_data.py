#from data import raw
import logging
import sys
from pathlib import Path
import subprocess
import yaml

# Set up logging
log_file = Path(__file__).parent / "preprocess_monitor.log"
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(message)s')

def read_config(config_path):
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)

def run_script(script_path, args=[]):
    """
    Runs a Python script with the provided arguments.
    
    Parameters:
    - script_path (str): Path to the Python script.
    - args (list): List of arguments to pass to the script.
    """
    result = subprocess.run(['python3', script_path] + args)
    if result.returncode != 0:
        logging.error(f"Error running {script_path} with args {args}")
        print(f"Error running {script_path} with args {args}")
        sys.exit(result.returncode)
    logging.info(f"Script {script_path} ran successfully with args {args}")

def check_and_convert_dataset(directory_path, converter_script,pdbqt_output_path):
    """
    Check each file in the directory for its format and convert if necessary.

    Parameters:
    - directory_path (str): Path to the directory containing dataset files.
    - converter_script (str): Path to the conversion script (e.g., 'sdf_pdbqt.py').
    """
    
    # Loop through each file in the directory
    for file_path in Path(directory_path).glob('*'):
        # Check if it's a file (not a directory)
        if file_path.is_file():
            # Get the file extension
            file_extension = file_path.suffix.lower()

            if file_extension == '.sdf':
                logging.info(f"Detected SDF file format for {file_path}. Converting to PDBQT format.")
                try:
                    # Run the conversion script
                    run_script(converter_script, [str(file_path), str(pdbqt_output_path)])
                except subprocess.CalledProcessError as e:
                    logging.error(f"Error during conversion of {file_path}: {e}")

if __name__ == "__main__":
    module_dir = Path(__file__).parent
    config_path = module_dir / "../configure.yaml"  # Path to your configuration file
    config = read_config(config_path)
    
    # Get dataset path and converter script from the config file
    dataset_path = config['task']['preprocess']['dataset_sdf_input']
    pdbqt_output_path = config['task']['preprocess']['dataset_pdbqt']
    converter_script = module_dir / '../docking/operations/sdf_pdbqt.py'
    f=config['task']['preprocess']['sdf_flag']
    # Check and convert dataset files
    if f == 1:
        check_and_convert_dataset(dataset_path, converter_script,pdbqt_output_path)
    elif f == 0:
        print("sdf_flag is set to 0, change it to 1 if you have sdf files and want them covert to pdbqt, ignore if you already have pdbqt")
    else:
        print("it can only be 1 or 0: if you have sdf then 1, else 0")


'''
import argparse

def preprocess_data(receptor_file, ligand_file, output_file):
    # Your preprocessing logic here
    print(f"Preprocessing receptor file: {receptor_file}")
    print(f"Preprocessing ligand file: {ligand_file}")
    print(f"Output will be saved to: {output_file}")
    # Add your preprocessing code here

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess data")
    parser.add_argument('--input', nargs=2, required=True, metavar=('RECEPTOR', 'LIGAND'), help="Specify the receptor and ligand input files")
    parser.add_argument('--output', required=True, metavar='OUTPUT', help="Specify the output file")
    args = parser.parse_args()
    
    preprocess_data(args.input[0], args.input[1], args.output)
'''