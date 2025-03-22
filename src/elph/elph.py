################### Electron Phonon Coupling molecular semiconductors (Pentacene, Rubrene) #####################
# Author: Ray
# Begin date: 2024/11/14
# 1st modified: Mapping for atoms is not correct, supercell size in phonopy modulation should match the supercell size (here is 2x2x2)
# Cannot use .xyz file (since .xyz file gives you many atoms outside unitcell)

import ase.io
import cclib
import json
import glob
import os
import numpy as np
import networkx as nx
import phonopy
import subprocess
import sys
import matplotlib.pyplot as plt
import elph.utils as ut
from time import time
from pathlib import Path
from phonopy.cui.create_force_sets import create_FORCE_SETS
from ase.calculators.gaussian import Gaussian
from ase.visualize import view
from ase.neighborlist import natural_cutoffs, neighbor_list
from ase import Atoms
from ase.build import sort
from scipy import sparse 
from scipy.constants import hbar, k
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
    file = glob.glob(path + "/*.cif") + glob.glob(path + "/*.vasp") + glob.glob(path + "/CONTCAR") + glob.glob(path + "/*.xyz")
    if len(file) == 0:
        raise FileNotFoundError
    
    return file[0]
    
def phonon(natoms,mesh,sc):
    """ Obtain FORCE_CONSTANTS file for specific calculator from Phonopy and return normal mode frequencies;
        Run Phonopy modulation to create atomic displacement and return modulation and frequency.
    Args:
    natoms (int): Number of atoms in the system (Defaults to None)
    mesh (list): Need define a mesh grid. (Defaults to [8,8,8])
    sc (list): Need define a supercell size. (Defaults to [2,2,2])
    #########################################
    Output:
    phonopy_params.yaml file
    mesh.yaml file (frequency)
    mod, freq: displacement_list and normal mode frequencies
    """
    if os.path.exists("FORCE_SETS"):
        print(" FORCE_SETS already exists, continue calculation! ")
        force_sets = "./FORCE_SETS"
    
    else:
        ut.print_error(" FORCE_SETS file is missing ")
        sys.exit(1)
    
    phonon = phonopy.load('phonopy_disp.yaml')
    phonon.run_mesh(mesh,with_eigenvectors=True)
    mesh_data = phonon.get_mesh_dict()
    freq = mesh_data['frequencies'].flatten()
    qpts = mesh_data['qpoints']
    eigenv = mesh_data['eigenvectors']
    a=sc[0]
    b=sc[1]
    c=sc[2]
    mod = np.zeros((int(freq.shape[0]*qpts.shape[0]),natoms*a*b*c,3)) # modulations

    mode = [[q, band_index, 1, 0.0] for q in qpts for band_index in range(natoms*3)]
 
    phonon.run_modulations(dimension=(a,b,c),phonon_modes=mode) 
    modulation, supercell = phonon.get_modulations_and_supercell()
    mod = np.real(modulation)

    index = np.where(freq>0)
    freq = freq[index]
    eigenv = eigenv[index]
    mod = mod[index]

    np.savetxt("frequencies.txt", freq, header="Phonon frequencies (THz)") # save frequencies to txt file
    np.savetxt("eigvectors.txt", eigenv, header="Eigenvectors") # save eigenvectors
    print(" Finish Modulation using Phonopy ")
    
    return mod, freq

def neighbor(atoms):
    """ Use ase and networkx to find the neighbors in the crystal structure and return the molecules
    Args: 
    atoms: ASE atoms object
    """
    i,j,S = neighbor_list(quantities='ijS', a=atoms, cutoff=natural_cutoffs(atoms)) # i: atom index, j: neighbor index, S: pbc

    atom_index = []
    neighbor_index = []

    for a, b, s in zip(i,j,S):
        if sum(s) == 0:
            atom_index.append(int(a))
            neighbor_index.append(int(b))

    G = nx.Graph() # Initiate networkx
    for i, j in zip(atom_index, neighbor_index):
        G.add_edge(int(i), int(j))

    molecules = list(nx.connected_components(G))

    return molecules

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
    
    # Group atoms by their molecule index
    neighbors = neighbor(atoms) 
    molecules = defaultdict(list) # empty dic for saving molecule POS
    for mol_idx, atom_idx in enumerate(neighbors): # Loop through every molecules and 
        molecules[mol_idx].extend(list(atom_idx))

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
            os.makedirs(os.path.join(f'./displacements', disp_name), exist_ok=True)
            os.chdir(f'./displacements/{disp_name}')
            ase.io.write(disp_name+'.xyz', disp_atoms)
            os.chdir('../..')
	   
        os.chdir(main_path) 

def mol_orbital(bset, atoms=None):
    """ Run Gaussian to compute the molecular orbitals for the system
    Args:
    opt (int): Run Gaussian to optimize the structure, Defaults to 0 (without optimization)
    bset (str): Basis set for Gaussian calculation (Defaults to 3-21G*)
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
                              basis=f'{bset}', # can use 6-31G* 
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
    output.decode('ascii').split()[-2]: Transfer integral J_eff,ab (Effective Transfer Integral) for the system, units is eV
    output.decode('ascii').split()[-13]: Transfer integral J_ab for the system, units is eV
    """
    cmd = f"calc_J -p_1 {path1} -p_2 {path2} -p_P {path3} -l_1 {path4} -l_2 {path5} -l_P {path6}"
    output = subprocess.check_output(cmd, shell=True)
    
    return output.decode('ascii').split()[-2], output.decode('ascii').split()[-13]

def onsiteE(homo):
    """ Use cclib to parse the Gaussian log file and get the onsite energy for the system.
    Args:
    homo (bootlean): Defaults to True, most OSCs are p-type (hole carrier transport in HOMO)
    """
    filename = 'mo.log'
    data = cclib.io.ccread(filename)
    moenergy = data.moenergies[0]
    if homo:
        homo_index = data.homos 
        onsite_eng = moenergy[homo_index]
    else:
        lumo_index = data.homos + 1
        onsite_eng = moenergy[lumo_index]

    return onsite_eng

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

def variance(freqs, g, temp):
    """ Calculate the variance of the transfer integral J
    Apply the formula: Var(J) = <J^2> - <J>^2
    Args:
    freqs: The numpy array of phonon frequencies
    g: electron-phonon coupling matrix
    temp: Temperture in Kelvin
    """
    var = (g**2/2)/np.tanh((hbar*freqs*1e12)/(2*k*temp)) # freqs in Phonopy is THz, so need to convert to Hz

    return var

        

    
 