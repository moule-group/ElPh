############## Visualization for electron-phonon coupling simulation ################
import ase.io
import glob
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import os
from ase.neighborlist import natural_cutoffs, neighbor_list
from scipy import sparse
from collections import OrderedDict, defaultdict

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
    file = glob.glob(path + "/*.cif") + glob.glob(path + "/*.xyz")
    if len(file) == 0:
        raise FileNotFoundError
    
    return file[0]

def visualize():
    """ Unwrap the crystal structure and visualize the molecules in crystal file
    Return:
    A plot for visualization
    """
    main_path = os.getcwd() # Main directory 
    atoms = ase.io.read(getGeometry(main_path)) # Load structure
    atoms *= (2,2,2)

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
    natoms = 0
    for v in molecule_positions.values():
        if len(v) > natoms: # only full molecules will be added 
            natoms = len(v)

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