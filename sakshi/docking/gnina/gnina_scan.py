#!/usr/bin/env python

import os
import subprocess
import yaml
from pathlib import Path

# Load configuration
config_file = Path(__file__).parent / "../../configure.yaml"  # Adjust this path if necessary

def load_config(config_file):
    """Loads the YAML configuration file."""
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

def run_gnina_docking(receptor, ligand_dir, output_dir, grid_center, grid_size, exhaustiveness, n_poses, gpu, gpu_id):
    """Runs GNINA docking for each ligand."""
    ligands = list(Path(ligand_dir).glob("*.pdbqt"))  # Find all .pdbqt ligands
    receptor = Path(receptor)  # Receptor file

    if not ligands:
        print("No ligands found in the specified directory!")
        return

    for ligand in ligands:
        ligand_name = ligand.stem
        output_file = Path(output_dir) / f"{ligand_name}_out.pdbqt"
        
        # Specify the path to GNINA executable
        gnina_executable = "/blue/yanjun.li/vi.gade1/seabra-li/docking/gnina/gnina"
        
        # Build the GNINA command
        gnina_command = [
            gnina_executable,  # Full path to the GNINA executable
            "--receptor", str(receptor),
            "--ligand", str(ligand),
            "--center_x", str(grid_center[0]),
            "--center_y", str(grid_center[1]),
            "--center_z", str(grid_center[2]),
            "--size_x", str(grid_size[0]),
            "--size_y", str(grid_size[1]),
            "--size_z", str(grid_size[2]),
            "--exhaustiveness", str(exhaustiveness),
            "--num_modes", str(n_poses),
            "--out", str(output_file)
        ]
        
        # Add GPU flag if enabled
        if gpu:
            gnina_command.extend(["--device", str(gpu_id)])

        # Optional: Reduce CNN memory usage
        gnina_command.extend(["--cnn_resolution", "1.0"])  # Increase CNN resolution to reduce memory usage
        gnina_command.extend(["--cnn_scoring", "rescore"])  # Use only CNN rescoring instead of full refinement
        
        # Run the docking
        print(f"Running GNINA docking for {ligand_name}...")
        result = subprocess.run(gnina_command, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error docking {ligand_name}: {result.stderr}")
        else:
            print(f"Docking completed for {ligand_name}. Results saved to {output_file}")

if __name__ == "__main__":
    # Load configuration
    config = load_config(config_file)

    # Control variables from the config file
    exhaustiveness = config['task']['docking']['gnina']['control']['exhaustiveness']
    n_poses = config['task']['docking']['gnina']['control']['n_poses']
    gpu = config['task']['docking']['gnina']['control']['gpu']
    gpu_id = config['task']['docking']['gnina']['control']['gpu_id']

    # Grid definition
    center_x = config['task']['docking']['gnina']['grid']['center_x']
    center_y = config['task']['docking']['gnina']['grid']['center_y']
    center_z = config['task']['docking']['gnina']['grid']['center_z']

    size_x = config['task']['docking']['gnina']['grid']['size_x']
    size_y = config['task']['docking']['gnina']['grid']['size_y']
    size_z = config['task']['docking']['gnina']['grid']['size_z']

    # Paths from the config file
    receptor_file = config['task']['receptor']['protein_file']
    ligand_dir = config['task']['preprocess']['dataset_pdbqt']
    output_dir = config['task']['docking']['gnina']['docking_output_gnina']

    # Ensure the output directory exists
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Run GNINA docking
    run_gnina_docking(
        receptor_file,
        ligand_dir,
        output_dir,
        [center_x, center_y, center_z],
        [size_x, size_y, size_z],
        exhaustiveness,
        n_poses,
        gpu,
        gpu_id
    )
