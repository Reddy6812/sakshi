#!/usr/bin/env python

import multiprocessing as mp
from ctypes import c_int
from functools import partial
from pathlib import Path
import shutil
import time
import yaml  # Import yaml to read the configuration file
from vina import Vina
from tqdm import tqdm
import torch  # Import torch to check for GPU availability
import sys  # For error handling

# Get the current module directory
module_dir = Path(__file__).parent

def load_config(config_file):
    """Loads the configuration YAML file."""
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

def dock_ligands(target, docked_dir, output_dir, center, size, ligand, n_cpus, exh, n_poses, e_range):
    """ Creates decoys for a given complex using Vina. """
    ligand_filename = Path(ligand).name
    print(f"Processing ligand: {ligand_filename}")

    # Check if this molecule has already been docked
    existing_files = {f.name for f in docked_dir.glob('*.pdbqt')}  # Use a set for better performance
    docked_filename = docked_dir / ligand_filename

    if ligand_filename in existing_files:
        with counter_lock:
            n_old_mols.value += 1
        result = output_dir / ligand_filename
        result.symlink_to(docked_filename)
        print(f"Ligand {ligand_filename} already docked, symlink created.")
    else:
        with counter_lock:
            n_new_mols.value += 1

        target_file = Path(target)
        ligand_file = Path(ligand)
        ligand_name = ligand_file.stem

        try:
            # Dock this molecule using Vina
            docker = Vina(sf_name='Vina', cpu=n_cpus, verbosity=0)
            docker.set_receptor(str(target_file))
            docker.set_ligand_from_file(str(ligand_file))
            docker.compute_vina_maps(center, size)
            docker.dock(exhaustiveness=exh, n_poses=n_poses)

            # Save the docking results
            docker.write_poses(f"{str(output_dir)}/{ligand_name}.pdbqt", n_poses, energy_range=e_range, overwrite=True)

            # Move the file to the collection dir and avoid unnecessary symlink if moved successfully
            result = output_dir / ligand_filename
            shutil.move(result, docked_filename)
            print(f"Ligand {ligand_filename} docked and results saved.")
        except Exception as e:
            print(f"Error docking ligand {ligand_filename}: {e}", file=sys.stderr)

    return

if __name__ == "__main__":

    start = time.time()  # Initialize the start time

    print(module_dir)
    # Load configuration
    config_path = module_dir / "../../configure.yaml"  # Adjust this path if necessary
    config = load_config(config_path)

    # Control variables
    EXHAUSTIVENESS = config['task']['docking']['vina']['control']['exhaustiveness']
    N_POSES = config['task']['docking']['vina']['control']['n_poses']
    PARALLEL_PROCESSES = config['task']['docking']['vina']['control']['parallel_processes']
    VINA_CPUS = config['task']['docking']['vina']['control']['vina_cpus']
    ENERGY_RANGE = config['task']['docking']['vina']['control']['energy_range']

    # Grid Definition
    center_x = config['task']['docking']['vina']['grid']['center_x']
    center_y = config['task']['docking']['vina']['grid']['center_y']
    center_z = config['task']['docking']['vina']['grid']['center_z']

    size_x = config['task']['docking']['vina']['grid']['size_x']
    size_y = config['task']['docking']['vina']['grid']['size_y']
    size_z = config['task']['docking']['vina']['grid']['size_z']

    print("Docking library with Vina:")
    print(f"  Exhaustiveness={EXHAUSTIVENESS}, N poses={N_POSES}")
    print(f"  Parallel processes={PARALLEL_PROCESSES}, Vina CPUs={VINA_CPUS}, Energy range={ENERGY_RANGE}")

    # Define paths
    ligand_dir = Path(config['task']['preprocess']['dataset_pdbqt'])
    receptor_file = Path(config['task']['receptor']['protein_file'])

    output_dir = Path(config['task']['docking']['vina']['docking_output_vina'])
    output_dir.mkdir(exist_ok=True)
    docked_dir = Path(module_dir, "dockedenamine")
    docked_dir.mkdir(exist_ok=True)

    ligands = list(ligand_dir.glob("*.pdbqt"))

    print(f"  Found: Target: {receptor_file} \n"
          f"         {len(ligands)} ligands in {ligand_dir}.\n")

    c = 1
    # Keeps track of how many molecules were already docked previously
    n_new_mols = mp.Value(c_int, 0)
    n_old_mols = mp.Value(c_int, 0)
    counter_lock = mp.Lock()

    # Generate decoys for each complex. Use multiprocessing to speed up the process.
    dock_ligands_partial = partial(dock_ligands, 
                                    receptor_file,
                                    docked_dir,
                                    output_dir,
                                    [center_x, center_y, center_z],
                                    [size_x, size_y, size_z],
                                    n_cpus=VINA_CPUS,
                                    exh=EXHAUSTIVENESS,
                                    n_poses=N_POSES,
                                    e_range=ENERGY_RANGE)
    
    # Start multiprocessing pool
    with mp.Pool(processes=PARALLEL_PROCESSES) as pool:
        with tqdm(total=len(ligands)) as pbar:
            for _ in pool.imap_unordered(dock_ligands_partial, ligands):
                pbar.update()

    # Calculate elapsed time
    finish = time.time()
    elapsed = finish - start

    print(f"From a total of {len(ligands)} molecules in the library:")
    print(f"    - {n_old_mols.value} have already been docked previously, and")
    print(f"    - {n_new_mols.value} have been docked now.\n")
    print(f"Docking process completed in {elapsed:.2f} seconds.")
