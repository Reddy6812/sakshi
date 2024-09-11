import os
import argparse
import shutil
import yaml
import sys
# Define the directory structure
directory_structure = {
    "ULVS": {
        "set-selection": {
            "isim": ["isim.py"]
        },
        "docking": {
            "vina": ["vina.py"],
            "autodock": ["AutoDock.py"],
            "gnina": ["gnina.py"],
            "add": ["Add_file.py"]
        },
        "rescoring": ["rescoring.py"],
        "ml_models": {
            "random_forest": ["random_forest.py"],
            "chembert": ["chembert.py"],
            "pinns": ["pinns.py"]
        },
        "properties": {
            "fingerprints": ["Fingerprints.py"],
            "logp": ["LogP.py"]
        },
        "db_scan": ["DB_Scan.py"]
    },
    "data": {
        "raw": [],
        "processed": ["feature_selected_data.csv", "combined_data.csv"],
        "docking": {
            "vina_scan": {
                "input": {
                    "ligands": [],
                    "protein.pdbqt": []
                },
                "output": {
                    "vina_results": [],
                    "pdbqt_files": []
                }
            },
            "deepatom": {
                "input": {
                    "chembert_input": ["chembert_input.smi"],
                    "unique_pdbqt_files": []
                },
                "output": {
                    "deepatom_results": ["deepatom_results.csv"],
                    "deepatom_rescored_pdb_files": []
                }
            },
            "chembert": {
                "input": ["chembert_input_da.smi"],
                "output": ["chembert_scores.csv"]
            }
        },
        "top_results": {
            "vina": ["top_results.csv", "top_pdbqt_files"],
            "deepatom": ["top_results.csv", "top_pdb_files"]
        }
    },
    "scripts": {
        "src": ["get_vina_results.py","prepare_chembert_input_from_vina.py","prepare_chembert_input_from_DA.py","select_top_results_Vina.py","select_top_results_DA.py","get_deepatom_results.py"],
        "analysis": ["preprocess_data.py", "feature_extraction.py", "docking_analysis.py", "drug_repurposing_analysis.py", "pinn_analysis.py", "model_training.py", "model_testing.py", "store_docking_results.py", "store_model_results.py", "get_top_scoring_results.py", "store_top_results.py"]
    },
    "results": ["model_accuracy.txt", "top_candidates.csv", "analysis_reports", "random_forest_model.joblib", "chembert_finetuned"],
    "configs": ["vina_config.yml","gnina_config.yml","AutoDock.yml"],
}

# Function to create the directory structure
def create_directory_structure(base_path):
    for key, value in directory_structure.items():
        create_subdirectories(os.path.join(base_path, key), value)

def create_subdirectories(base_path, structure):
    if isinstance(structure, dict):
        os.makedirs(base_path, exist_ok=True)
        for key, value in structure.items():
            create_subdirectories(os.path.join(base_path, key), value)
    elif isinstance(structure, list):
        os.makedirs(base_path, exist_ok=True)
        for item in structure:
            if isinstance(item, str):
                open(os.path.join(base_path, item), 'a').close()
            elif isinstance(item, list):
                os.makedirs(os.path.join(base_path, item[0]), exist_ok=True)

# Functions to add and delete files and folders
def add_file(directory, filename):
    file_path = os.path.join(directory, filename)
    print(f"Creating file at: {file_path}")

    if os.path.isdir(file_path):
        raise IsADirectoryError(f"The specified path {file_path} is a directory, not a file.")

    os.makedirs(directory, exist_ok=True)
    with open(file_path, 'a') as f:
        pass
    print(f"File {filename} created in {directory}")

def add_folder(directory, foldername):
    os.makedirs(os.path.join(directory, foldername), exist_ok=True)

def delete_file(filepath):
    if os.path.isfile(filepath):
        os.remove(filepath)

def delete_folder(folderpath):
    if os.path.isdir(folderpath):
        shutil.rmtree(folderpath)

# Function to load and check parameters
def load_and_check_parameters(config_file):
    with open(config_file, 'r') as file:
        params = yaml.safe_load(file)

    # Example parameter checks
    if 'email' not in params or not params['email']:
        raise ValueError("Email parameter is missing or empty")

    if 'partition' not in params or not params['partition']:
        raise ValueError("Partition parameter is missing or empty")

    if 'qos' not in params or not params['qos']:
        raise ValueError("QoS parameter is missing or empty")

    if 'gpus' not in params or params['gpus'] <= 0:
        raise ValueError("GPUs parameter is missing or invalid")

    if 'cpus_per_task' not in params or params['cpus_per_task'] <= 0:
        raise ValueError("CPUs per task parameter is missing or invalid")

    if 'time_limit' not in params or not params['time_limit']:
        raise ValueError("Time limit parameter is missing or empty")

    return params

# Function to handle input files
def handle_input(input_file, destination_dir):
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"The specified input file {input_file} does not exist.")
    os.makedirs(destination_dir, exist_ok=True)
    shutil.copy(input_file, destination_dir)
    print(f"Copied input file {input_file} to {destination_dir}")

# Function to handle output files
def handle_output(output_file, destination_dir):
    os.makedirs(destination_dir, exist_ok=True)
    file_path = os.path.join(destination_dir, output_file)
    open(file_path, 'a').close()
    print(f"Created output file {file_path}")

if __name__ == "__main__":
    print("hi from utils")
    def cpu_type(x):
        x = int(x)
        if x <= 0:
            raise argparse.ArgumentTypeError("Number of CPUs must be a positive integer")
        return x

    def filepath_type(x):
        if not os.path.exists(x):
            raise argparse.ArgumentTypeError(f"File {x} does not exist")
        return x

    parser = argparse.ArgumentParser(description="Manage project structure and parameters")
    parser.add_argument('--interactive', action='store_true', help="Run in interactive mode")
    parser.add_argument('-f', '--add-file', nargs=2, metavar=('DIRECTORY', 'FILENAME'), help="Add a file to a directory")
    #parser.add_argument('-afd', '--add-folder', nargs=2, metavar=('DIRECTORY', 'FOLDERNAME'), help="Add a folder to a directory")
    parser.add_argument('-d', '--delete-file', metavar='FILEPATH', help="Delete a file")
    #parser.add_argument('-dfd', '--delete-folder', metavar='FOLDERPATH', help="Delete a folder")
    #parser.add_argument('-cs', '--create-structure', metavar='BASEPATH', help="Create the directory structure")
    parser.add_argument('--config', metavar='CONFIGFILE', help="Path to the configuration file")
    parser.add_argument('-i', '--input', nargs=2, metavar=('INPUTFILE', 'DESTINATION'), help="Specify an input file and its destination directory")
    parser.add_argument('-o', '--output', nargs=2, metavar=('OUTPUTFILE', 'DESTINATION'), help="Specify an output file and its destination directory")
    parser.add_argument('-c', '--ncpu', type=cpu_type, help='Number of CPUs to use')
    parser.add_argument('-p', '--program', choices=['vina', 'autodock', 'gnina'], help="Select docking program")
    parser.add_argument('-m', '--ml-model', choices=['random_forest', 'chembert', 'pinns'], help="Select machine learning model")
    # Allow any string for the stage argument, we'll validate it manually
    parser.add_argument('-s', '--stage', type=str, help="Select stages (e.g., 0,1,2,preprocess,docking, operations, ml-train,ml-inference)")
    #parser.add_argument('-t', '--test', choices=['0', '1'], help="Select if you want to test before continue: 0=No, 1=Yes")
    
    args = parser.parse_args()
    
    # Split the stages by comma if provided and validate them
    if args.stage:
        stages = args.stage.split(',')
        for stage in stages:
            if stage not in ['0', '1', '2', '3', '4', 'preprocess','docking','operations', 'ml-train','ml-inference']:
                print(f"Invalid stage: {stage}")
                sys.exit(1)
    else:
        stages = []

    # Example processing based on selected stages
    if stages:
        for stage in stages:
            print(f"Selected stage: {stage}")
            
            
    #if args.create_structure:
    #    create_directory_structure(args.create_structure)

    if args.add_file:
        add_file(args.add_file[0], args.add_file[1])
    
    #if args.add_folder:
    #    add_folder(args.add_folder[0], args.add_folder[1])
    
    if args.delete_file:
        delete_file(args.delete_file)

    #if args.delete_folder:
        delete_folder(args.delete_folder)

    #if args.config:
        params = load_and_check_parameters(args.config)
        print("Loaded parameters:", params)
        
    if args.input:
        handle_input(args.input[0], args.input[1])

    if args.output:
        handle_output(args.output[0], args.output[1])

    if args.program:
        print(f"Selected docking program: {args.program}")

    if args.ml_model:
        print(f"Selected machine learning model: {args.ml_model}")
        
    if args.ncpu:
        print(f"Number of CPUs to use: {args.ncpu}")

    if args.stage:
        print(f"Selected stage: {args.stage}")

    if args.interactive:
        while True:
            action = input("Choose an action: [add-file, add-folder, delete-file, delete-folder, create-structure, input, output, program, ml-model, stage, exit]: ")
            if action == 'exit':
                break
            elif action == 'add-file':
                directory = input("Enter directory path: ")
                filename = input("Enter filename: ")
                add_file(directory, filename)
            elif action == 'add-folder':
                directory = input("Enter directory path: ")
                foldername = input("Enter foldernamet: ")
                add_folder(directory, foldername)
            elif action == 'delete-file':
                filepath = input("Enter file path: ")
                delete_file(filepath)
            elif action == 'delete-folder':
                folderpath = input("Enter folder path: ")
                delete_folder(folderpath)
            elif action == 'create-structure':
                basepath = input("Enter base path: ")
                create_directory_structure(basepath)
            elif action == 'input':
                input_file = input("Enter input file path: ")
                destination = input("Enter destination directory: ")
                handle_input(input_file, destination)
            elif action == 'output':
                output_file = input("Enter output file name: ")
                destination = input("Enter destination directory: ")
                handle_output(output_file, destination)
            elif action == 'program':
                program = input("Select docking program [vina, autodock, gnina]: ")
                print(f"Selected docking program: {program}")
            elif action == 'ml-model':
                ml_model = input("Select machine learning model [random_forest, chembert, pinns]: ")
                print(f"Selected machine learning model: {ml_model}")
            elif action == 'stage':
                stage = input("Select starting stage [0, 1, 2, 3, 4, preprocess,docking,operations,ml-train,ml-inference: ")
                print(f"Selected stage: {stage}")
