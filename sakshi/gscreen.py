import subprocess
import sys
from pathlib import Path

# Get the current module directory
module_dir = Path(__file__).parent
import logging
# Set up logging
log_file = Path(__file__).parent /"monitor.log"
logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s - %(message)s')


def call_utils(args):
    # Check if help flag is present in the arguments
    if '-h' in args or '--help' in args:
        # Call the utils.py script with the help flag
        result = subprocess.run(['python3', module_dir / 'utils.py', '--help'], text=True)
        sys.exit(0)  # Exit after showing help
    else:
        # Call the utils.py script with the provided arguments
        print(module_dir)
        result = subprocess.run(['python3', module_dir / 'utils.py'] + args, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"Error in utils.py: {result.stderr}")
            sys.exit(result.returncode)
        return result.stdout

def process_parameters(output):
    # Process the output from utils.py to get the selected parameters
    params = {}
    lines = output.split('\n')
    for line in lines:
        if line.startswith("Selected docking program:"):
            params['program'] = line.split(":")[1].strip()
        elif line.startswith("Selected machine learning model:"):
            params['ml_model'] = line.split(":")[1].strip()
        elif line.startswith("Number of CPUs to use:"):
            params['ncpu'] = int(line.split(":")[1].strip())
        elif line.startswith("Selected stage:"):
            # Split the stage into a list
            params['stage'] = line.split(":")[1].strip().split(',')
        
        logging.info(f"Parameters processed: {params}")
        
    return params

def run_script(script_path, args=[]):
    result = subprocess.run(['python3', script_path] + args)
    if result.returncode != 0:
        logging.error(f"Error running {script_path}")
        print(f"Error running {script_path}")
        sys.exit(result.returncode)
    logging.info(f"Script {script_path} ran successfully")

def stage_0(params):
    print("*******Preprocessing data started")
    logging.info("Stage 0: Preprocessing data started")
    # Step 2: Preprocess/Clean Data
    run_script(module_dir / 'preprocess/preprocess_data.py')
    print("*****Preprocessing data completed")

def stage_1(params):
    # Start from Docking stage
    # Step 3: Start Docking
    print("\n *****Starting docking:   ")
    #if 'program' in params and 'ncpu' in params:
    #    run_script(module_dir / f'ULVS/docking/{params["program"]}/{params["program"]}.py', [f'--ncpu={params["ncpu"]}'])
    if 'program' in params:
        run_script(module_dir / f'docking/{params["program"]}/{params["program"]}_scan.py')
    # Step 4: Read Docking Data
    #run_script(module_dir / 'scripts/src/docking_analysis.py')
    print("\n ***** docking completed**** ")
    # Step 5: Store Docking Results
    #run_script(module_dir / 'scripts/src/store_docking_results.py')
    
def stage_2(params):
    
    
    print("\n *****Running Operations: ")
    run_script(module_dir / f'docking/operations/operations.py')

    # Step 10: Store Top Results in different file(create new file in src)
    #run_script(module_dir / 'scripts/src/store_top_results.py')
    print("\n *****Operations completed******* ")

    
    
def stage_3(params):
    print("\n *****Starting ml-models:   ")
    # Start from ML Models stage
    if 'ml_model' in params:
        run_script(module_dir / f'ml/training/{params["ml_model"]}/{params["ml_model"]}.py')
    else:
        print("Skipping ML Model stage: 'ml_model' parameter not provided")

    print("\n *****ml-models completed**** ")
    
def stage_4(params):
    print("\n *****Starting ml-Inferance:   ")
    # Start from ML Models stage
    if 'ml_model' in params:
        run_script(module_dir / f'ml/inference/{params["ml_model"]}/{params["ml_model"]}.py')
    else:
        print("Skipping ML Model stage: 'ml_model' parameter not provided")

    print("\n *****ml-Inferance completed and extracting top results**** ")
    
def main():
    # Get all the command line arguments except the script name
    utils_args = sys.argv[1:]

    # Call utils.py and get the parameters
    output = call_utils(utils_args)
    params = process_parameters(output)

    # Debug print
    print(f"Parameters: {params}")

    # Determine the stages to execute
    stages = params.get('stage')
    
    if stages:
        print(f"Stages to execute: {stages}")
        
        #'preprocess','docking','operations', 'ml-train','ml-inference'
        for stage in stages:
            if stage == '0':
                stage_0(params)
            elif stage == 'preprocess':
                stage_0(params)
            elif stage == '1':
                stage_1(params)
            elif stage == 'docking':
                stage_1(params)
            elif stage == '2':
                stage_2(params)
            elif stage == 'operations':
                stage_2(params)
            elif stage == '3':
                stage_3(params)
            elif stage == 'ml-train':
                stage_3(params)
            elif stage == '4':
                stage_4(params)
            elif stage == 'ml-inference':
                stage_4(params)
            else:
                print(f"Unknown stage: {stage}")
                sys.exit(1)
    else:
        print("No stages specified.")
        sys.exit(1)

if __name__ == "__main__":
    main()
