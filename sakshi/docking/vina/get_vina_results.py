#!/usr/bin/env python
"""
This script reads the PDBQT files with results from Vina/QVina, and outputs
a '.csv' file with the columns:

    'File':       File name
    'MoleculeID': Molecule ID number
    'Scores':     Vina Scores

ATTENTION: this only works well with pdbqt files with a SINGLE MODEL. If there are
more models in the file, there is no guarantee which model is being returned!
"""



# -- General
import numpy as np
import pandas as pd
from pathlib import Path
import sys
from tqdm import tqdm
import shutil

#-----------------------------------
def parse_pdbqt(pdbqt_file):
    """Gets information from the FIRST MODEL in the  PDBQT file """
    mol = {'name':'', 'vina_affinity': 0.0}
    with open(pdbqt_file, 'r') as infile:
        for line in infile:
            if line[:12] == 'REMARK  Name':
                mol['name'] = line[15:].strip()
            elif line[:18] == 'REMARK VINA RESULT':
                mol['vina_affinity'] = float(line[19:].split()[0])
            elif 'ENDMDL' in line:
                break
    return mol

# -------------------
if __name__ == '__main__':


    lig_files = []
    lig_names = []
    scores    = []
    troubled_mols = []

    pdbqt_dirs = [d.resolve() for d in list(Path('./').glob('**')) if 'out' in str(d)]
    print("Reading docking results from PDBQT files in:")
    for d in pdbqt_dirs:
        print(f"\t{d}")

    with tqdm(pdbqt_dirs) as pbar_dirs:
        for pdbqt_dir in pbar_dirs:
            pbar_dirs.set_description(f"Processing dir {str(pdbqt_dir)}")

            ligand_files = list(Path(pdbqt_dir).glob('*.pdbqt'))

            with tqdm(ligand_files) as pbar_files:

                for lig_file in pbar_files:
                    pbar_files.set_description(f"Processing file {str(lig_file.stem)}")
                    this_lig = parse_pdbqt(lig_file)

                    if this_lig is not None:
                        lig_files.append(lig_file)
                        lig_names.append(this_lig['name'])
                        scores.append(this_lig['vina_affinity'])
    
                    else:
                        troubled_mols.append(f"{lig_file} \t {failure} \n")
    
    results_df = pd.DataFrame({'File':lig_files,'Variant':lig_names, 'Scores':scores})
    results_df['MoleculeID'] = results_df.Variant.str.split('-').str[0]
    results_df = results_df.sort_values(by=['Scores','MoleculeID'])
    results_df.to_csv('docking_results_all.csv', index=False)
    
    # Troubled Molecules
    # Those are the molecules rdkit could not read for some reason.
    if len(troubled_mols) > 0:
        troubled_file = "troubled_mols.txt"
        print(f"Found {len(troubled_mols)} problematic molecules. Writing to file {troubled_file}")
        with open(troubled_file,'w') as of:
            for mol in troubled_mols:
                of.write(str(mol))

    # For every MoleculeID, select only the best result
    results_df.drop_duplicates(subset='MoleculeID', keep='first', inplace=True, ignore_index=True)
    results_df.to_csv('docking_results_unique_molid.csv', index=False)

    # Get n_best molecules for next step
    # n_best=1_000
    # results_df[:n_best].to_csv('top_1k.csv', index=False)
    # Path('./top_1k').mkdir(exist_ok=True)
    # for docked_file in tqdm(results_df.File[:n_best]):
    #    origin = Path(docked_file)
    #    destination = Path('top_1k',origin.name)
    #    shutil.copy(origin,destination)


