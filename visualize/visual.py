############## Visualization for electron-phonon coupling simulation ################
import ase.io
import glob
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import os
import sys
from ase.neighborlist import natural_cutoffs, neighbor_list
from scipy import sparse
from collections import OrderedDict, defaultdict
from itertools import product

def getGeometry(path):
    """ Using glob function in python to find the structure
    file in the current path
    The type of the structure files: ".cif"
    Args:
    path: The current directory (use os.getcwd())
    #########################################
    Return:
    file: The structure file in the path
    """
    file = glob.glob(path + "/*.cif") + glob.glob(path + "/*.vasp") + glob.glob(path + "/CONTCAR") + glob.glob(path + "/*.xyz")
    if len(file) == 0:
        raise FileNotFoundError
    
    return file[0]

def extract_Z_from_cif(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if '_cell_formula_units_Z' in line:
                return int(line.split()[-1])
    return 1  # Fallback

def scale_incrementer():
    base = (1, 1, 1)
    
    while True:
        if all(b >= 20 for b in base):
            sys.exit("Reached maximum scale value: (20, 20, 20)")
        for delta in product([0, 1], repeat=3):
            yield tuple(b + d for b, d in zip(base, delta))
        base = tuple(b + 1 for b in base)

def visualize():
    """ Unwrap the crystal structure and visualize the molecules in crystal file
    Return:
    A plot for visualization
    """
    main_path = os.getcwd() # Main directory 
    atomic_structure = ase.io.read(getGeometry(main_path)) # Load structure
    formula_units_Z = extract_Z_from_cif(getGeometry(main_path))
    atoms_per_molecule = len(atomic_structure) // int(formula_units_Z)
    print(f'Atoms in molecule: {atoms_per_molecule}')

    scale_gen = scale_incrementer()
    scale = next(scale_gen)
    natoms = 0
    mol_count = 0
    while natoms != atoms_per_molecule or mol_count < 12:
        print(f'\nAttempting atom scale of {scale}')
        atoms = atomic_structure * scale
        i,j,S = neighbor_list(quantities='ijS', a=atoms, cutoff=natural_cutoffs(atoms)) # i: atom index, j: neighbor index, S: pbc

        atom_index = []
        neighbor_index = []

        for a, b, s in zip(i,j,S):
            if sum(s) == 0: # Only consider the nearest neighbors (remove periodic boundary condition)
                atom_index.append(int(a))
                neighbor_index.append(int(b))

        G = nx.Graph() # Initiate networkx
        for i, j in zip(atom_index, neighbor_index):
            G.add_edge(int(i), int(j))

        molecules = list(nx.connected_components(G))
        molecule_positions = {} # empty dic for saving molecule POS
    
        for mol_idx, atom_idx in enumerate(molecules): # Loop through every molecules and store positions of the molecule
            molecule = atoms[list(atom_idx)]
            molecule_positions[mol_idx] = molecule.get_positions() 

        # Find the correct number of atoms for 1 molecule
        for v in molecule_positions.values():
            if len(v) > natoms: # only full molecules will be added 
                natoms = len(v)
        print(f'Atoms found: {natoms}')

        mol_count = 0
        for positions in molecule_positions.values():
            if len(positions) == natoms:
                mol_count += 1
        print(f'Molecule count: {mol_count}')

        scale = next(scale_gen)

    print(f'\nSuccess!')
    # Plot each molecule with a label
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    markers = ['o','*','v','p']
    for mol_idx, positions in molecule_positions.items():
        if len(positions) == natoms: # Only full molecules will be plotted
            positions = np.array(positions)
            ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], 
                    label=f"Molecule {mol_idx + 1}",
                    marker=markers[mol_idx%len(markers)])

    # Labels and visualization options
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")
    ax.set_title("Molecule Numbering Visualization")
    ax.legend()
    plt.show()

if __name__ == "__main__":
    visualize()