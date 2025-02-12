############## Visualization for electron-phonon coupling simulation ################

import ase.io
import matplotlib.pyplot as plt
import numpy as np
import elph as ep
import os
from collections import OrderedDict, defaultdict


def visualize(structure_path):
    """ Unwrap the crystal structure and visualize the molecules in crystal file
    Args:
    structure_path (str): structure file (make sure crystal is large enough to capture 
    the 3 nearest neightbors)
    #############################################
    Return:
    A plot for visualization
    """
    atoms = ase.io.read(structure_path) # Load structure
    atoms *= (2,2,2)
    n_components, component_list = ep.neighbor_list(atoms)
    
    # Group atoms by their molecule index
    molecules = defaultdict(list)
    for atom_idx, mol_idx in enumerate(component_list):
        molecules[mol_idx].append(atom_idx)

    molecule_positions = {} # empty dic for saving molecule POS
    # Save each molecule to a separate file
    for mol_idx, atom_idx in molecules.items():
        molecule = atoms[atom_idx]
        # Store positions of the molecule
        molecule_positions[mol_idx] = molecule.get_positions()

    # Prepare figure
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot each molecule with a label
    markers = ['o','*','v','p']
    for mol_idx, positions in molecule_positions.items():
        positions = np.array(positions)
        ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], 
                   label=f"Molecule {mol_idx + 1}",
                   marker=markers[mol_idx%len(markers)])
        # Annotate molecule positions with numbers
        #center = positions.mean(axis=0)
        #ax.text(center[0], center[1], center[2], f"{mol_idx + 1}", fontsize=12, color='red', fontweight='bold')

    # Labels and visualization options
    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")
    ax.set_title("Molecule Numbering Visualization")
    ax.legend()
    plt.show()

if __name__ == '__main__':
    path = os.getcwd()
    file = ep.getGeometry(path)
    visualize(file)