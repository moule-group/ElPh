#!/usr/bin/env python3
import ase.io
import glob
import os
import numpy as np
import networkx as nx
import string
import sys
from itertools import combinations
from scipy.spatial.distance import pdist, squareform
from ase.calculators.gaussian import Gaussian
from ase.neighborlist import natural_cutoffs, neighbor_list, NeighborList
from scipy import sparse   
from pathlib import Path
import time

def number_to_letters(n):
    """Convert 0 -> A, 1 -> B, ..., 25 -> Z, 26 -> AA, ..."""
    result = ""
    n += 1
    while n > 0:
        n, rem = divmod(n - 1, 26)
        result = chr(65 + rem) + result
    return result

def getGeometry(path):
    """ Using glob function in python to find the structure
    file in the current path
    The type of the structure files: "POSCAR", "CONTCAR", "*.xyz"
    Args:
        path: The current directory (use os.getcwd())
    ----------------------------------------------
    Return:
        file: The structure file in the path
    """
    file = glob.glob(path + "/POSCAR") + glob.glob(path + "/CONTCAR") + glob.glob(path + "/*.xyz") 
    if len(file) == 0:
        raise FileNotFoundError
    return file[0]

def mol_in_cell(atoms):
    """ Get the number of atoms in the unit cell
    Args:
        atoms (ASE Atoms object): molecule or dimer 
    -----------------------------------------------
    Return:
        nmol_in_cell: number of molecules in the unit cell
    """
    cutoff = natural_cutoffs(atoms)
    nl = NeighborList(cutoff, self_interaction=False, bothways=True)
    nl.update(atoms)
    matrix = nl.get_connectivity_matrix()
    nmol_in_cell, _ = sparse.csgraph.connected_components(matrix) # nmol_in_cell is how many molecules in unitcell
    return nmol_in_cell

def neighbor(atoms_unitcell, supercell_array, nmols):
    """ Use ase and networkx to find the neighbors in the crystal structure and return the molecules
    Args: 
        atoms_unitcell: Atoms from the unitcell
        supercell_array: Supercell size
        nmols: num of molecules that are extracted
    -----------------------------------------------
    Return:
        atoms : ASE atoms object with supercell
        full_mols: list of full molecules with atomic indices 
        nearest_idx: list of indices of the nearest neighbors
        coms_array: array of center of mass for each full molecule
    """
    natoms_in_cell = len(atoms_unitcell) # number of atoms in the unit cell
    nmol_in_cell = mol_in_cell(atoms_unitcell) # number of molecules in the unit cell
    atoms = atoms_unitcell * supercell_array # atoms in the supercell
    cutoff = natural_cutoffs(atoms) # The cutoff between atoms
    full_mols = [] # The empty list to save a complete molecules (sometimes molecules would be cut)

    print("Start to find the molecules in the supercell!!!\n"
          f"The expected atom count for each molecule: {natoms_in_cell / nmol_in_cell}.")
    attempt = 0
    max_attempts = 3
    while attempt < max_attempts:
        i,j,S = neighbor_list(quantities='ijS', a=atoms, cutoff=cutoff) # i: atom 1 index, j: atom 2 index, S: shift vector
        atom_index = []
        neighbor_index = []

        for a, b, s in zip(i,j,S):
            if sum(s) == 0: # This means 2 atoms are in the same cell, and it must be neighbors
                atom_index.append(int(a))
                neighbor_index.append(int(b))

        G = nx.Graph() # Initiate networkx
        for i, j in zip(atom_index, neighbor_index):
            G.add_edge(int(i), int(j)) 

        molecules = list(nx.connected_components(G)) # These are molecules
        for mol in molecules:
            if len(mol) == (natoms_in_cell / nmol_in_cell):
                full_mols.append(mol)

        if len(full_mols) >= nmols: # If full molecules are found, break the loop
            break
        
        attempt += 1
        cutoff_h = [np.float64(x+0.03*attempt) if x == 0.31 else x for x in natural_cutoffs(atoms)] # increase Hydrogen cutoff
        cutoff_c = [np.float64(x+0.03*attempt) if x == 0.76 else x for x in cutoff_h] # increase Carbon cutoff
        cutoff_n = [np.float64(x+0.03*attempt) if x == 0.71 else x for x in cutoff_c] # increase Nitrogen cutoff
        cutoff_s = [np.float64(x+0.03*attempt) if x == 1.05 else x for x in cutoff_n] # increase Sulfur cutoff
        cutoff_si = [np.float64(x+0.03*attempt) if x == 1.11 else x for x in cutoff_s] # increase silicon cutoff
        cutoff = [np.float64(x+0.03*attempt) if x == 0.57 else x for x in cutoff_si] # increase flourine cutoff

    if not full_mols:
        print("Failed to find any full molecules after all retries. Exiting.")
        sys.exit(1)

    coms = []  # center of mass
    for mol_indices in full_mols: # calculate the center of mass for each full molecule
        mol = atoms[list(mol_indices)]  # extract the molecule
        com = mol.get_center_of_mass()
        coms.append(com)

    coms_array = np.array(coms)  # shape (n_mols, 3)
    distance_matrix = squareform(pdist(coms_array)) # Calulate distance matrix 

    distances_matrix_0 = distance_matrix[0] # reference molecule (select first row)
    distances_matrix_0[0] = np.inf # set the diagonal to infinity to ignore self-distance
    min_dist = np.min(distances_matrix_0) # The minimum distance between reference molecule and the other molecules
    num_min_dist = len(np.where(distances_matrix_0 == min_dist)[0]) # Number of the molecules that are closest to the reference molecule

    sorted_idx = np.argsort(distances_matrix_0) # sort the distance from short to long by index
    start = num_min_dist - 1 # sometimes it happens the nearest neighbors have 2, so instead of taking twice, we pick only 1
    stop = start + (nmols - 1) # grab nmols-1 neighbors
    nearest_idx = sorted_idx[start:stop] # get the indices of the nearest neighbors
    nearest_idx = np.insert(nearest_idx, 0, 0)  # insert the first molecule (itself) at the beginning of the list

    np.savez_compressed('center_of_mass', coms_array[nearest_idx,:]) 
    return atoms, full_mols, nearest_idx, coms_array

def map_to_middle(pos, cell):
    """ Map the coordinates to the middle of the unit cell
    Args:
        pos (np.ndarray): Coordinates of the atoms
        cell (np.ndarray): Cell vectors
    -----------------------------------------------
    Return:
        pos_cell (np.ndarray): Coordinates of the atoms wrapped into [0, 1)
    """
    inv_cell = np.linalg.inv(cell)
    pos_cell = pos @ inv_cell
    pos_cell = pos_cell - np.floor(pos_cell) # wrap cell coordinates into [0, 1)
    return pos_cell

def mapping_atom(new_atoms_pos, cell, old_atoms_pos, tol):
    """ Mapping atoms in .xyz file with POSCAR unitcell for the following phonopy modulation
    Args:
        new_atoms_pos (np.ndarray): atomic coordinates for the atoms to be mapped
        cell (np.ndarray): Cell vectors of the unitcell
        old_atoms_pos (np.ndarray): Coordinates (scaled) of the atoms in the unit cell
        tol (float): Tolerance for the mapping (Defaults to 1e-4)
    -----------------------------------------------
    Return:
        mapping (list): List of indices of the atoms in the unit cell that resemble to the atoms in the supercell
    """
    mapping = []
    new_atoms_scaled_pos = map_to_middle(new_atoms_pos, cell)
    for pos in new_atoms_scaled_pos:
        diffs = np.abs(old_atoms_pos - pos) # Compute distance to all unit cell atoms
        dists = np.linalg.norm(diffs, axis=1) # The distance matrix
        idx = np.argmin(dists)
        if dists[idx] < tol:
            mapping.append(int(idx))
        else:
            mapping.append(None)  # or raise error / warning
    return mapping

def unwrap_molecule_dimer(structure_path, supercell_array, nmols):
    """ Get single molecule and molecular pairs (dimer) files.
    Args:
        structure_path (str): structure file path
        supercell_martix (tuple): supercell size
        nmols (int): number of molecules that are extracted
    -----------------------------------------------
    Return:
        molecule_{x}.xyz file, where x is the numbering 
        dimer_{A}.xyz, where A is the labeling 
    """
    label_list = []
    atoms_unitcell = ase.io.read(structure_path) # Load structure
    cell = atoms_unitcell.get_cell() # Get cell vectors of the unit cell
    atoms_unitcell_pos = atoms_unitcell.get_scaled_positions() # Get scaled positions of the atoms in the unit cell
    atoms, full_mols, nearest_idx, _ = neighbor(atoms_unitcell, supercell_array, nmols) # full_mols is a list save the atomic index for each molecule

    allmols_index = np.concatenate([
    list(full_mols[idx]) for idx in nearest_idx[:nmols]
    ]) # This finds nearest neighbors molecules, each atomic index is written in the list

    newmol = atoms[allmols_index]
    ase.io.write('allmols.xyz', newmol) # Check the geometry of the molecules
   
    for i in range(nmols):
        os.mkdir(f'{i+1}')
        label_list.append(f'{i+1}')
        name_mol = os.path.join(str(i+1), f"monomer_{i+1}.xyz")
        atoms_id = list(full_mols[nearest_idx[i]])
        mol = atoms[atoms_id]
        mol.set_pbc((False, False, False))
        ase.io.write(name_mol, mol)

    pairs = list(combinations(nearest_idx, 2)) 
    for j, pair in enumerate(pairs):
        letter = number_to_letters(j)
        os.mkdir(letter)
        label_list.append(letter)
        name_dimer = os.path.join(letter, f"dimer_{letter}.xyz")
        atoms_id1 = list(full_mols[pair[0]])
        atoms_id2 = list(full_mols[pair[1]])
        dimer = atoms[atoms_id1] + atoms[atoms_id2]
        dimer.set_pbc((False, False, False))
        ase.io.write(name_dimer, dimer) 

        mapping = mapping_atom(dimer.get_positions(), cell, atoms_unitcell_pos, tol=1e-4)
        np.savez_compressed(os.path.join('mapping', f'map_{letter}.npz'), mapping=mapping) # Save the mapping of atoms

    return label_list

def mol_orbital(bset, functional, disper_corr):
    """ Run Gaussian to compute the molecular orbitals coefficient and energy for the system (Single point calculation)
    Args:
        bset (str): Basis set for Gaussian calculation 
        functional (str): Functional for Gaussian calculation 
        disper_corr (int): 1 is D3, 2 is D3BJ
    -----------------------------------------------
    Return:
        .pun file which contains molecular orbital for the cacluation later for J 
        .log file output file from Gaussian
    """
    path = os.getcwd()
    geometry = getGeometry(path)
    atoms = ase.io.read(geometry)

    if disper_corr == 0:
        dispersion = ''
    elif disper_corr == 1:
        dispersion = 'EmpiricalDispersion=GD3'
    else:
        dispersion = 'EmpiricalDispersion=GD3BJ'
    
    atoms.calc = Gaussian(mem='24GB',
                 nprocshared=12,
                 label='mo',
                 method=functional,
                 basis=bset,  
                 scf='tight',
                 pop='full',
                 extra=f'nosymm punch=mo iop(3/33=1) {dispersion}') # iop(3/33=1) output one-electron integrals to log file.
    atoms.calc.write_input(atoms)

def write_run_script():
    """ Write the run script in each subdirectory for Gaussian simulation
    """
    string = """#!/bin/bash
folder=$(basename "$PWD")
g16 < mo.com > mo.log
mv fort.7 "${folder}.pun"
    """
    run_file = Path("run.sh")
    run_file.write_text(string)
    os.chmod(run_file, 0o755)

def run_init(basis_set, functional, dispersion_correction, supercell_array, nmols):
    """ Main function to generate the initial files
    Args:
        basis_set (str): DFT Basis sets
        functional (str): DFT functional
        dispersion_correction (int): Dispersion correction (0-> no dispersion, 1-> D3, 2-> D3-BJ)
        supercell_array (tuple): The supercell matrix (need large enough to find the neighbors)
        nmols (int): Number of molecules to construct nearest neighbors
    ------------------------------------------------------------
    Return:
       Initials files for Gaussian calculations (.com) and .xyz files for monomers and dimers.
    """    
    path = os.getcwd() # Main directory which contain all subfolders
    j0_file = glob.glob(os.path.join(path, 'j', 'j0_eff.json'))
    xyz_file = glob.glob(os.path.join(path, '1', 'monomer_1.xyz'))
    if not j0_file: # If j_0.json is not exists, run the following simulation
        try:
            geometry = getGeometry(path) # Get the geometry file
            os.makedirs(os.path.join(path, 'j'), exist_ok=True) # Create a directory for J_ij

        except FileNotFoundError:
            print("Crystal structure (CONTCAR; POSCAR ...) file not found in the current directory. Exiting.") 
            sys.exit(1)

        if not xyz_file:
            os.makedirs(os.path.join(path, 'mapping'), exist_ok=True) # Create a directory for J_ij
            label_list = unwrap_molecule_dimer(geometry, supercell_array, nmols) # Unwrap the crystal to get single molecule and dimers
    
        for folder in label_list:
            os.chdir(folder)
            mol_orbital(basis_set, functional, dispersion_correction) 
            write_run_script()
            os.chdir(path)

    print(" --- Initial Files Generations Finished! --- ")

def ask_with_default(prompt, default, type=str):
    """ Ask user for input with a default value. If user presses Enter, return default.
    Args:
        prompt (str): The prompt to display to the user
        default: The default value to return if user presses Enter
        type: The type to cast the user input to (default is str)
    """
    user_input = input(f"{prompt} [{default}]: ").strip()
    if user_input == "":
        return default
    return type(user_input)

if __name__ == "__main__":

    print(" ----------------------------------------------------------------- ")
    print(" Starting on Electron Phonon Coupling (ElPh) for Organic Semiconductors !!! ")
    print(" ElPh code made by Moule Group at UC Davis. Author: Ray Chiang ")
    print(" #########   $           $ ########   $        $  ")
    print(" $           $           $        $   $        $  ")
    print(" $########   $           $ #######    $########$  ")
    print(" $           $           $            $        $  ")
    print(" #########   ##########  $            $        $  ")
    print(" ----------------------------------------------------------------- ")

    print("This script allows to generate the initial files for J_eff calculation\n" 
          "Gaussian input files for molecular orbital calculation will be generated in each folder\n")

    nmols = ask_with_default(
    "Number of molecules to construct nearest neighbors",
    4,int)

    supercell_input = ask_with_default(
    "Supercell array (comma separated, e.g. 4,4,4)",
    "4,4,4",str)

    supercell_array = [int(x) for x in supercell_input.split(",")]

    basis_set = ask_with_default(
    "Basis set. Example: Def2SVP, TZVP, Def2TZVP .... Default is TZVP",
    "TZVP")

    functional = ask_with_default(
    "Functional. Example: B3LYP, PBE0, M062X .... Default is M062X",
    "M062X")

    dispersion_correction = ask_with_default(
    "Dispersion correction (0-> no dispersion, 1-> D3, 2-> D3BJ)",
    1,int)

    print("\n"
          "Input parameters:\n"
          " ----------------------- \n"
          f"Number of molecules to construct nearest neighbors: {nmols}\n"
          f"Basis sets: {basis_set}\n"
          f"Functional: {functional}\n"
          f"Dispersion correction (0-> no dispersion, 1-> D3, 2-> D3-BJ): {dispersion_correction}\n"
          " ----------------------- \n"
    )

    run_init(basis_set, functional, dispersion_correction, supercell_array, nmols) 

    print("")
    print(" --------------------------------------------------------- ")
    print(
    """                 
    #########  ##    #  #########
    $          # #   #  $       $
    ########   #  #  #  $      $
    $          #   # #  $     $ 
    #########  #    ##  ######
    """)
    print(" --------------------------------------------------------- ")
    print(time.ctime())
