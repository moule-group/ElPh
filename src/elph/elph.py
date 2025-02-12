################### Electron Phonon Coupling molecular semiconductors (Pentacene, Rubrene) #####################
# Author: Ray
# Begin date: 2024/11/14
# 1st modified: Mapping for atoms is not correct, supercell size in phonopy modulation should be 1 1 1;
# Cannot use .xyz file (since .xyz file gives you many atoms outside unitcell)

import ase.io
import json
import argparse
import glob
import os
import numpy as np
import phonopy
import subprocess
import matplotlib.pyplot as plt
import elph.utils as ut
from time import time
from pathlib import Path
from phonopy.cui.create_force_sets import create_FORCE_SETS
from ase.calculators.gaussian import Gaussian
from ase.visualize import view
from ase.neighborlist import natural_cutoffs, NeighborList
from ase import Atoms
from ase.build import sort
from scipy import sparse 
from collections import OrderedDict, defaultdict
    
def phonon(mesh):
    """ Obtain FORCE_CONSTANTS file for specific calculator from Phonopy and return normal mode frequencies;
        Run Phonopy modulation to create atomic displacement and return displacement_list.
    Args:
    mesh (list): Need define a mesh grid. (Defaults to [8,8,8])
    #########################################
    Output:
    phonopy_params.yaml file
    mesh.yaml file (frequency)
    displacement_list
    """
    if os.path.exists("FORCE_SETS"):
        print(" FORCE_SETS already exists, continue calculation! ")
        force_sets = "./FORCE_SETS"
    
    else:
        ut.print_error(" FORCE_SETS file is missing ")
        sys.exit(1)
    
    phonon = phonopy.load('phonopy_disp.yaml')
    force_constant = phonon.produce_force_constants()
    print(" Creating phonopy_params.yaml file which saves force constants. ")
    phonon.save(settings={'force_constants': True})
    print(" Done creating phonopy_params.yaml. ")
    
    phonon.run_mesh(mesh,with_eigenvectors=True)
    print(" Creating mesh.yaml file which saves normal mode frequencies. ")
    phonon.write_yaml_mesh() 
    print(" Done creating mesh.yaml. ")
    
    mesh_data = phonon.get_mesh_dict()
    
    phonon_modes = [] # Empty list to save phonon modes
    displacement_list = [] # Empty list to save displacement
    
    amplitude = 0.01  # Displacement amplitude (Defaults to 0.01 according to phonopy)
    phase = 0.0 # Phase factor (Defaults to 0 according to phonopy)
    
    print(" Starting Modulation using Phonopy ")
    for i, q in enumerate(mesh_data['qpoints']):
        frequencies = mesh_data['frequencies'][i]
        eigenvectors = mesh_data['eigenvectors'][i]
    
        for band_index in range(len(frequencies)):
            mode = [q, band_index, amplitude, phase] 
            mode_ = [mode]
            phonon_modes.append(mode_)
    
    for pm in phonon_modes:
        phonon.run_modulations(dimension=[1,1,1],phonon_modes=pm) # Change here, dimension should be [1x1x1] to fit the number of phonon modes.
        modulation, supercell = phonon.get_modulations_and_supercell()
        displacement = modulation
        displacement_list.append(np.real(displacement))
        
    print(" Finish Modulation using Phonopy ")
    
    return displacement_list

def neighbor_list(atoms):
    """ Find the neighbors for molecules in molecular crystals
    Args:
    atoms: ase Atoms object
    #########################################
    Return:
    n_components: number of molecules
    component_list: contains molecule number of each atom in the system
    """
    nl = NeighborList(natural_cutoffs(atoms), self_interaction=False, bothways=True) # Determine the indices of each molecule using NeighborList
    nl.update(atoms)

    conn_matrix = nl.get_connectivity_matrix(False) # Connectivity matrix
    
    n_components, component_list = sparse.csgraph.connected_components(conn_matrix) 

    return n_components, component_list

def unwrap_molecule_dimer(structure_path, mol1, mol2, mol3):
    """ Get single molecule and molecular pairs (dimer) files.
    Args:
    structure_path (str): structure file path
    mol1: The numbering of first neighbor molecule
    mol2: The numbering of second neighbor molecule
    mol3: The numbering of third neighbor molecule (translation of center molecule)
    Return:
    molecule_{x}.xyz file, where x is the numbering (3 files)
    dimer_{A}.xyz, where A is the labeling (3 files)
    """
    atoms = ase.io.read(structure_path) # Load structure
    atoms *= (2,2,2) # Supercell
    n_components, component_list = neighbor_list(atoms)
    
    # Group atoms by their molecule index
    molecules = defaultdict(list)
    for atom_idx, mol_idx in enumerate(component_list):
        molecules[mol_idx].append(atom_idx)

    os.mkdir('1')
    os.mkdir('2')
    os.mkdir('3')
    name_mol1 = "1/monomer_1.xyz"
    name_mol2 = "2/monomer_2.xyz"
    name_mol3 = "3/monomer_3.xyz"
    
    mono1 = atoms[molecules[mol1-1]]
    mono2 = atoms[molecules[mol2-1]]
    mono3 = atoms[molecules[mol3-1]]

    mono1.set_pbc((False, False, False))
    mono2.set_pbc((False, False, False))
    mono3.set_pbc((False, False, False))
    
    ase.io.write(name_mol1, mono1)
    ase.io.write(name_mol2, mono2)
    ase.io.write(name_mol3, mono3)
    
    pairs = {'1':molecules[mol1-1],
             '2':molecules[mol2-1],
             '3':molecules[mol3-1]}
    
    with open('monomer.json', 'w', encoding='utf-8') as f: # Create atom mapping for electron-phonon coupling term (match phonons)
        json.dump(pairs, f, ensure_ascii=False, indent=4)

    os.mkdir('A')
    os.mkdir('B')
    os.mkdir('C')
    name_dimerA = "A/dimer_A.xyz"
    name_dimerB = "B/dimer_B.xyz"
    name_dimerC = "C/dimer_C.xyz"

    di1 = atoms[molecules[mol1-1]] + atoms[molecules[mol2-1]]
    di2 = atoms[molecules[mol1-1]] + atoms[molecules[mol3-1]]
    di3 = atoms[molecules[mol2-1]] + atoms[molecules[mol3-1]]

    di1.set_pbc((False, False, False))
    di2.set_pbc((False, False, False))
    di3.set_pbc((False, False, False))
    
    ase.io.write(name_dimerA, di1) 
    ase.io.write(name_dimerB, di2)
    ase.io.write(name_dimerC, di3)
    
    ########### Create atom mapping for electron-phonon coupling term (match phonons) ###########
    #counter = 0
    #atom_mapping = {}
    #idx_1 = molecules[mol1-1]
    #idx_2 = molecules[mol2-1]
    #idx_3 = molecules[mol3-1]
    
    #idx_total = idx_1 + idx_2 + idx_3
    
    #for idx in idx_total:
       # atom_mapping[idx] = counter
       # counter += 1

    #with open('atom_mapping.json', 'w') as f:
       # f.write(json.dumps(OrderedDict(sorted(atom_mapping.items(), key=lambda t: t[1])), indent=2))
        
def get_displacement(atoms):
    """ Get numbering of displaced atom, displacement direction (x,y,z) and sign (+,-) 
    Args:
    atoms (ASE Atom object): molecule or dimer 
    #########################################
    Return:
    na: numbering of displaced atom
    vec: displacement direction (x,y,z) in (0,1,2)
    sign: (+,-) in (+1, -1)
    """
    num_atoms = len(atoms) # number of atoms
    for na in range(num_atoms):
        for vec in range(3):
            for sign in [-1, 1]:
                yield (na,vec,sign) 
        
def create_displacement(delta=0.01):
    """ Create atomic displacement and return each displaced structures
    Args: 
    delta (float): Magnitude of displacement in Angstrom (Defaults to 0.01A)
    """
    main_path = os.getcwd()
    
    folder_list = ['1', '2', '3', 'A', 'B', 'C']

    # Create the folders under the "displacements" directory
    for folder in folder_list:
        os.makedirs(os.path.join(folder, 'displacements'),exist_ok=True)
    print(" Finish creating displacements folder! ")
    
    mol_list = ['1/molecule_1.xyz', '2/molecule_2.xyz', '3/molecule_3.xyz', 'A/dimer_A.xyz', 'B/dimer_B.xyz', 'C/dimer_C.xyz']
    
    for folder in folder_list:
	
        os.chdir(folder)
        path = os.getcwd() # current path (1/; 2/; 3/; A/; B/; C/)
        print(f" Currently simulated in folder {path} ")
        
        geometry = getGeometry(path)
        atoms = ase.io.read(geometry)
        print(f" Now using this molecule file {geometry} ... ")
        
        for na,vec,sign in get_displacement(atoms):
            disp_atoms = atoms.copy()
            pos = disp_atoms.get_positions()
            pos[na, vec] += delta * sign
            disp_atoms.set_positions(pos)
	
            disp_name = f"disp_{na}_{vec}_{sign}"
            os.makedirs(os.path.join(f'./displacements', disp_name),exist_ok=True)
            os.chdir(f'./displacements/{disp_name}')
            ase.io.write(disp_name+'.xyz', disp_atoms)
            os.chdir('../..')
	   
        os.chdir(main_path) 

def mol_orbital(atoms=None):
    """ Run Gaussian to compute the molecular orbitals for the system
    Args:
    atoms (ASE atoms object): optional, if not specified, it will run getGeometry 
    ########################################
    Return:
    .pun file which contains molecular orbital for the cacluation later for J 
    """
    if not atoms:
        path = os.getcwd()
        print(f" Now working in this directory {path} ... ")
        geometry = getGeometry(path)
        atoms = ase.io.read(geometry)
    
    pun_files = glob.glob('*.pun')

    if not pun_files: # If there is no Gaussian output, it will run Gaussian
        atoms.calc = Gaussian(mem='16GB',
                              nprocshared=12,
                              label='mo',
                              save=None,
                              method='b3lyp',
                              basis='3-21G*',
                              scf='tight',
                              pop='full',
                              extra='nosymm punch=mo iop(3/33=1)') # iop(3/33=1) output one-electron integrals to log file.

        atoms.get_potential_energy()
        os.rename('fort.7', os.path.basename(path) + '.pun')
        
    print(f' Gaussian simulation for {path} molecular orbitals is done! ')

def run_catnip(path1, path2, path3, path4, path5, path6):
    """ Run Catnip to calculate the transfer integral
    Args:
    path1 (str): Path to the first Gaussian .pun file
    path2 (str): Path to the second Gaussian .pun file
    path3 (str): Path to the third Gaussian .pun file (pair)
    path4 (str): Path to the first Gaussian .log file
    path5 (str): Path to the second Gaussian .log file
    path6 (str): Path to the third Gaussian .log file (pair)
    ###################################################
    Returns:
    output.decode('ascii').split()[-2]: Transfer integral J_eff (Effective Transfer Integral) for the system, units is eV
    output.decode('ascii').split()[-13]: Transfer integral J for the system, units is eV
    """
    
    cmd = f"calc_J -p_1 {path1} -p_2 {path2} -p_P {path3} -l_1 {path4} -l_2 {path5} -l_P {path6}"
    output = subprocess.check_output(cmd, shell=True)
    
    return output.decode('ascii').split()[-2], output.decode('ascii').split()[-13]  # return J_eff and J

def get_deri_Jmatrix(j_list, delta=0.01):
    """ Calculate derivative of transfer integral J and return as electron-phonon coupling matrix 
    Args:
    j_list (list): list of transfer integrals for each dimer 
    delta (float): size of displacement (Defaults to 0.01 Angstrom)
    ########################################################
    Returns:
    matrix: matrix containing gradient of J_ij
    """
    num_atoms = len(j_list) // 6 # 6 because each atom have 6 displaced treatment, 3 directions (x,y,z) and 2 sign (+ & -)
    
    # array with j- (negative delta)
    j_n = np.empty([num_atoms, 3]) 
    j_n[:, 0] = j_list[0::6] # atom x -
    j_n[:, 1] = j_list[2::6] # atom y -
    j_n[:, 2] = j_list[4::6] # atom z -
    
    # array with j+ (positive delta)
    j_p = np.empty([num_atoms, 3])
    j_p[:, 0] = j_list[1::6] # atom x +
    j_p[:, 1] = j_list[3::6] # atom y +
    j_p[:, 2] = j_list[5::6] # atom z +

    matrix = (np.abs(j_p) - np.abs(j_n)) / (2 * delta)

    return matrix
        

    
 