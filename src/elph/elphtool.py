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
import string
import sys
import elph.utils as ut
from itertools import combinations
from scipy.spatial.distance import pdist, squareform
from ase.calculators.gaussian import Gaussian, GaussianOptimizer
from ase.neighborlist import natural_cutoffs, neighbor_list, NeighborList
from scipy import sparse 
from scipy.constants import h, k
from collections import OrderedDict, defaultdict

def getGeometry(path):
    """ Using glob function in python to find the structure
    file in the current path
    The type of the structure files: ".cif"
    Args:
    path: The current directory (use os.getcwd())
    ----------------------------------------------
    Return:
    file: The structure file in the path
    """
    file = glob.glob(path + "/POSCAR") +  glob.glob(path + "/*.xyz") 
    if len(file) == 0:
        raise FileNotFoundError
    
    return file[0]
    
def phonon(natoms, mesh):
    """ Obtain FORCE_CONSTANTS file for specific calculator from Phonopy and return normal mode frequencies;
        Run Phonopy modulation to create atomic displacement and return modulation and frequency.
    Args:
    natoms (int): Number of atoms in the system (Defaults to None)
    mesh (list): Need define a mesh grid. (Defaults to [8,8,8])
    ----------------------------------------------
    Output:
    phonopy_params.yaml file
    mod: displacement_list in unitcell 
    freq: normal mode frequencies in THz
    nqpts: number of q points
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
    freq = mesh_data['frequencies'].flatten() # Shape is (nqpts*natoms*3,)
    qpts = mesh_data['qpoints'] 
    nqpts = len(qpts) # number of q points
    mod = np.zeros((int(freq.shape[0]),natoms,3)) # modulations

    mode = [[q, band_index, 1, 0.0] for q in qpts for band_index in range(natoms*3)]
 
    phonon.run_modulations(dimension=(1,1,1), 
                           phonon_modes=mode) 
    modulation, _ = phonon.get_modulations_and_supercell()
    mod = np.real(modulation)

    index = np.where(freq>0)
    freq = freq[index]
    mod = mod[index]

    np.savetxt("frequencies.txt", freq, header="Phonon frequencies (THz)") # save frequencies to txt file
    print(" Finish Modulation using Phonopy ")
    
    return mod, freq, nqpts

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

def neighbor(atoms_unitcell, supercell_array, nmols=3):
    """ Use ase and networkx to find the neighbors in the crystal structure and return the molecules
    Args: 
    atoms: ASE atoms object
    supercell_array: Supercell size
    nmols: num of molecules that are extracted
    -----------------------------------------------
    Return:
    atoms : ASE atoms object with supercell
    full_mols: list of full molecules with atomic indices 
    nearest_idx: list of indices of the nearest neighbors
    """
    natoms_in_cell = len(atoms_unitcell) # number of atoms in the unit cell
    nmol_in_cell = mol_in_cell(atoms_unitcell) # number of molecules in the unit cell
    atoms = atoms_unitcell * supercell_array
    cutoff = natural_cutoffs(atoms)

    full_mols = []

    attempt = 0
    max_attempts = 3
    while attempt < max_attempts:

        i,j,S = neighbor_list(quantities='ijS', a=atoms, cutoff=cutoff) # i: atom index, j: neighbor index, S: pbc

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

        for mol in molecules:
            if len(mol) == (natoms_in_cell / nmol_in_cell):
                full_mols.append(mol)

        if len(full_mols) >= 3: # If full molecules are found, break the loop
            break
        
        attempt += 1
        ut.print_error("No molecules matched the expected atom count of {natoms_in_cell / nmol_in_cell}. Increase the cutoff distance and try again.")
        cutoff_h = [np.float64(x+0.02*attempt) if x == 0.31 else x for x in natural_cutoffs(atoms)] # increase Hydrogen cutoff
        cutoff_c = [np.float64(x+0.01*attempt) if x == 0.76 else x for x in cutoff_h] # increase Carbon cutoff
        cutoff_n = [np.float64(x+0.01*attempt) if x == 0.71 else x for x in cutoff_c] # increase Nitrogen cutoff
        cutoff = [np.float64(x+0.01*attempt) if x == 1.05 else x for x in cutoff_n] # increase Sulfur cutoff

    if not full_mols:
        ut.print_error("Failed to find any full molecules after all retries. Exiting.")
        sys.exit(1)

    coms = [] 
    for mol_indices in full_mols: # calculate the center of mass for each full molecule
        mol = atoms[list(mol_indices)]  # extract the molecule
        com = mol.get_center_of_mass()
        coms.append(com)

    coms_array = np.array(coms)  # shape (n_mols, 3)
    distance_matrix = squareform(pdist(coms_array)) # Calulate distance matrix 

    distances_matrix_0 = distance_matrix[0] # reference molecule (select first row)
    distances_matrix_0[0] = np.inf # set the diagonal to infinity to ignore self-distance
    nearest_idx = np.argsort(distances_matrix_0)[:nmols-1] # get the indices of the nearest neighbors
    nearest_idx = np.insert(nearest_idx, 0, 0) # insert the first molecule (itself) at the beginning of the list

    np.savez_compressed('center_of_mass', coms_array[nearest_idx,:]) 

    return atoms, full_mols, nearest_idx

def map_to_middle(coordinates, cell):
    """ Map the coordinates to the middle of the unit cell
    Args:
    coordinates (np.ndarray): Coordinates of the atoms in the supercell
    cell (np.ndarray): Cell vectors of the supercell
    -----------------------------------------------
    Return:
    coordinates_cell (np.ndarray): Coordinates of the atoms wrapped into [0, 1)
    """
    inv_cell = np.linalg.inv(cell)
    coordinates_cell = coordinates @ inv_cell
    
    coordinates_cell = coordinates_cell - np.floor(coordinates_cell) # wrap cell coordinates into [0, 1)

    return coordinates_cell

def mapping_atom(coordinates, cell, unitcell, tol=1e-4):
    """ Mapping atoms in .xyz file with POSCAR unitcell for the following phonopy modulation
    Args:
    cell (np.ndarray): Cell vectors of the unitcell
    coordinates (np.ndarray): Coordinates of (This ) the atoms in the supercell
    unitcell (np.ndarray): Coordinates (scaled) of the atoms in the unit cell
    tol (float): Tolerance for the mapping (Defaults to 1e-4)
    -----------------------------------------------
    Return:
    mapping (list): List of indices of the atoms in the unit cell that resemble to the atoms in the supercell
    """
    mapping = []
    scaled_coordinates = map_to_middle(coordinates, cell)
    for pos in scaled_coordinates:
        # Compute distance to all unit cell atoms
        diffs = np.abs(unitcell - pos)
        dists = np.linalg.norm(diffs, axis=1) # The distance matrix
        idx = np.argmin(dists)
        if dists[idx] < tol:
            mapping.append(int(idx))
        else:
            mapping.append(None)  # or raise error / warning
    
    return mapping

def unwrap_molecule_dimer(structure_path, supercell_array, nmols=3):
    """ Get single molecule and molecular pairs (dimer) files.
    Args:
    structure_path (str): structure file path
    supercell_martix (tuple): supercell size
    nmols (int): number of molecules that are extracted (Defaults to 3)
    -----------------------------------------------
    Return:
    molecule_{x}.xyz file, where x is the numbering (3 files)
    dimer_{A}.xyz, where A is the labeling (3 files)
    """
    atoms_unitcell = ase.io.read(structure_path) # Load structure
    cell = atoms_unitcell.get_cell() # Get cell vectors of the unit cell
    unitcell = atoms_unitcell.get_scaled_positions() # Get scaled positions of the atoms in the unit cell
    atoms, full_mols, nearest_idx = neighbor(atoms_unitcell, supercell_array) 

    allmols_index = np.concatenate((list(full_mols[nearest_idx[0]]),
                                    list(full_mols[nearest_idx[1]]),list(full_mols[nearest_idx[2]])))
    newmol = atoms[allmols_index]
    ase.io.write('allmols.xyz', newmol) # Check the geometry of the molecules
   
    for i in range(nmols):
        os.mkdir(f'{i+1}')
        name_mol = os.path.join(str(i+1), f"monomer_{i+1}.xyz")
        atoms_id = list(full_mols[nearest_idx[i]])
        mol = atoms[atoms_id]
        mol.set_pbc((False, False, False))
        ase.io.write(name_mol, mol)

    pairs = list(combinations(nearest_idx, 2)) 
    for j, letter in enumerate(string.ascii_uppercase[:len(pairs)]):
        os.mkdir(letter)
        name_dimer = os.path.join(letter, f"dimer_{letter}.xyz")
        atoms_id1 = list(full_mols[pairs[j][0]])
        atoms_id2 = list(full_mols[pairs[j][1]])
        dimer = atoms[atoms_id1] + atoms[atoms_id2]
        dimer.set_pbc((False, False, False))
        ase.io.write(name_dimer, dimer) 

        mapping = mapping_atom(dimer.get_positions(), cell, unitcell, tol=1e-4)
        np.savez_compressed(os.path.join('mapping', f'map_{letter}.npz'), letter=mapping) # Save the mapping of atoms
        
def get_displacement(atoms):
    """ Get numbering of displaced atom, displacement direction (x,y,z) and sign (+,-) 
    Args:
    atoms (ASE Atom object): molecule or dimer 
    -----------------------------------------------
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
        
def create_displacement(delta=0.01, nmols=3):
    """ Create atomic displacement and return each displaced structures
    Args: 
    delta (float): Magnitude of displacement in Angstrom (Defaults to 0.01A)
    nmols (int): Number of molecules that are extracted (Defaults to 3)
    """
    main_path = os.getcwd()
    
    # List of folders to create
    if nmols == 3:
        folder_list = ['1', '2', '3', 'A', 'B', 'C']

    # Create the folders under the "displacements" directory
    for folder in folder_list:
        os.makedirs(os.path.join(folder, 'displacements'), exist_ok=True)
    print(" Finish creating displacements folder! ")
    
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
            os.makedirs(os.path.join('./displacements', disp_name), exist_ok=True)
            os.chdir(f'./displacements/{disp_name}')
            ase.io.write(disp_name+'.xyz', disp_atoms)
            os.chdir('../..')
	   
        os.chdir(main_path) 

def gaussian_opt(atoms, bset, label, functional, ncharge=0):
    """ Run Gaussian simulation to get the onsite energy for the system
    Args:
    atoms (ASE atoms object)
    bset (list): Basis set for Gaussian calculation 
    label (str): Label for the Gaussian calculation
    functional (list): Functional for Gaussian calculation
    ncharged (int): Charge of the system (0: neutral, 1: hole transport, -1: electron transport)
    -----------------------------------------------
    """
    opt_cmd = Gaussian(mem='16GB',
                       chk=f'{label}.chk',
                       nprocshared=12,
                       label=label,
                       charge=ncharge,
                       mult=2*0.5*ncharge+1, # 2S+1
                       save=None,
                       method=functional[0],
                       basis=bset[0], 
                       scf='tight',
                       pop='full',
                       extra='nosymm EmpiricalDispersion=GD3BJ freq') 

    opt_calc = GaussianOptimizer(atoms, opt_cmd)
    opt_calc.run(fmax='tight',steps=60)
    ase.io.write(f'{label}.xyz', atoms)

def hr_factor(bset, functional):
    """ Create Gaussian input files for Huang-Rhys factor calculation and run Gaussian
    Args:
    bset (list): Basis set for Gaussian calculation
    functional (list): Functional for Gaussian calculation
    ------------------------------------------------
    """
    with open('hr_cation.com', 'w') as f1: # This is from neutral to cation
        cmds = ["%mem=16GB\n",
                "%oldchk=cation.chk\n",
                "%chk=fc.chk\n",
                f"#P {functional[0]}/{bset[0]} EmpiricalDispersion=GD3BJ Geom=AllCheck Freq=(READFC,FC,ReadFCHT)\n",
                "\n",
                "Initial=Source=Chk Final=Source=Chk\n",
                "print=(huangrhys,matrix=JK)\n",
                "\n",
                "neutral.chk\n",
                "cation.chk\n",
                "\n"
               ]
          
        f1.writelines(cmds)

    with open('hr_neutral.com', 'w') as f2: # This is from cation to neutral 
        cmds = ["%mem=16GB\n",
                "%oldchk=neutral.chk\n",
                "%chk=fc.chk\n",
                f"#P {functional[0]}/{bset[0]} EmpiricalDispersion=GD3BJ Geom=AllCheck Freq=(READFC,FC,ReadFCHT)\n",
                "\n",
                "Initial=Source=Chk Final=Source=Chk\n",
                "print=(huangrhys,matrix=JK)\n",
                "\n",
                "cation.chk\n",
                "neutral.chk\n",
                "\n"
               ]
          
        f2.writelines(cmds)

    subprocess.run('g16 < hr_cation.com > hr_cation.log', shell=True)
    subprocess.run('g16 < hr_neutral.com > hr_neutral.log', shell=True)

def parse_log(logfile1, logfile2):
    """ Parse the Gaussian log file to get (1) onsite energy (2) frequencies (3) Huang-Rhys factors
    (4) reorganization energy (5) local EPC
    Args:
    logfile1 (str): Path to the Gaussian log file (for cclib)
    logfile2 (str): Path to the Gaussian log file (for bottom part)
    -------------------------------------------
    Returns:
    eng (float): Onsite energy for the system (eV)
    freq (list): Frequencies of the system (cm^-1)
    huangrhys (list): Huang-Rhys factors for the system (unitless)
    reorg_eng (list): Reorganization energy for the system in eV
    gii (list): Square of local electron-phonon coupling for the system in (eV^2)
    """
    data = cclib.io.ccread(logfile1)
    moenergy = data.moenergies[0]
    homo_index = data.homos 
    eng = moenergy[homo_index]
    freqs = data.vibfreqs # units: cm^-1
    vibdisp_cart = data.vibdisps
    vibdisp_cart_squared = vibdisp_cart**2
    vibdisp_squared = np.sum(vibdisp_cart_squared, axis=1)

    jtoev = 6.241509074460763e+18 
    cm_1tohz = 3e10   
    huangrhys = []  
    gii_squared = []
    
    with open(logfile2) as log:
    
        for l in log.readlines():
            L = l.split()
    
            if len(L)==6 and L[0]=='Mode' and L[1]=='num.' and L[3]=='-' and L[4]=='Factor:':
                hr=L[5]
                str1=hr[0:8]
                str2=hr[-3:]
                hr=str1+'E'+str2
                huangrhys.append(float(hr)) # unitless         

        # lambda_i = hbar * w_i * S (S is Huang Rhys factor)
        reorg_eng = freqs * cm_1tohz * h * jtoev * huangrhys # Unit here should be eV

    # lambda_i = g_ii^2 / (hbar*w_i)
    gii_squared_cart = np.zeros((len(freqs), 3))
    for i in range(len(freqs)):
        gii_2 = reorg_eng[i] * h * jtoev * freqs[i] * cm_1tohz
        gii_squared.append(gii_2)
        denominator = np.sum(vibdisp_squared[i])
        gii_squared_x = gii_2 * vibdisp_squared[i][0] / denominator
        gii_squared_y = gii_2 * vibdisp_squared[i][1] / denominator
        gii_squared_z = gii_2 * vibdisp_squared[i][2] / denominator
        gii_squared_cart[i,:] = [gii_squared_x, gii_squared_y, gii_squared_z]

    return eng, freqs, huangrhys, reorg_eng, gii_squared, gii_squared_cart

def mol_orbital(bset, functional, atoms=None):
    """ Run Gaussian to compute the molecular orbitals coefficient and energy for the system (Single point calculation)
    Args:
    bset (list): Basis set for Gaussian calculation 
    functional (list): Functional for Gaussian calculation 
    atoms (ASE atoms object): optional, if not specified, it will run getGeometry 
    -----------------------------------------------
    Return:
    .pun file which contains molecular orbital for the cacluation later for J 
    .log file output file from Gaussian
    """
    if not atoms:
        path = os.getcwd()
        geometry = getGeometry(path)
        atoms = ase.io.read(geometry)
    
    pun_files = glob.glob('*.pun')
    if not pun_files: # If there is no Gaussian output, it will run Gaussian
        atoms.calc = Gaussian(mem='16GB',
                              nprocshared=12,
                              label='mo',
                              save=None,
                              method=functional[1],
                              basis=bset[1], # can use 6-31G* 
                              scf='tight',
                              pop='full',
                              extra='nosymm punch=mo EmpiricalDispersion=GD3BJ iop(3/33=1)') # iop(3/33=1) output one-electron integrals to log file.

        atoms.get_potential_energy()
        os.rename('fort.7', os.path.basename(path) + '.pun')

def run_catnip(path1, path2, path3, path4, path5, path6):
    """ Run Catnip to calculate the transfer integral
    Args:
    path1 (str): Path to the first Gaussian .pun file
    path2 (str): Path to the second Gaussian .pun file
    path3 (str): Path to the third Gaussian .pun file (pair)
    path4 (str): Path to the first Gaussian .log file
    path5 (str): Path to the second Gaussian .log file
    path6 (str): Path to the third Gaussian .log file (pair)
    ------------------------------------------------
    Returns:
    output.decode('ascii').split()[-2]: Transfer integral J_eff,ab (Effective Transfer Integral) for the system, units is eV
    output.decode('ascii').split()[-13]: Transfer integral J_ab for the system, units is eV
    """
    cmd = f"calc_J -p_1 {path1} -p_2 {path2} -p_P {path3} -l_1 {path4} -l_2 {path5} -l_P {path6}"
    output = subprocess.check_output(cmd, shell=True)
    
    return output.decode('ascii').split()[-2], output.decode('ascii').split()[-13]

def get_deri_Jmatrix(j_list, delta=0.01):
    """ Calculate derivative of transfer integral J and return as electron-phonon coupling matrix 
    Args:
    j_list (list): list of transfer integrals for each dimer 
    delta (float): size of displacement (Defaults to 0.01 Angstrom)
    ------------------------------------------------
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

def variance(freqs, g2, nqpts, temp, unit='THz'):
    """ Calculate the variance of the transfer integral J
    Apply the formula: Var(J) = <J^2> - <J>^2
    Args:
    freqs: The numpy array of phonon frequencies
    g2: squared electron-phonon coupling matrix
    nqpts (int): total qpts used in phonon calculation
    temp (float): Temperture in Kelvin
    unit (str): unit of the frequencies (THz or cm^-1, defaults to THz)
    ------------------------------------------------
    Returns:
    var (array): variance of the transfer integral J
    sigma (float): standard deviation of the transfer integral J (eV)
    """ 
    cm_1tohz = 3e10  
    if unit == 'THz':
        boseein = 1 / np.tanh((h*freqs*1e12)/(2*k*temp)) # Bose-Einstein distribution
        var = (g2/2) * boseein # freqs in Phonopy is THz, so need to convert to Hz
        sigma = (np.sum(var)/nqpts)**0.5 # Square root of variance, have to do normalization over the number of q points 

    if unit == 'cm-1':
        boseein = 1 / np.tanh((h*freqs*cm_1tohz)/(2*k*temp))
        var = (g2/2) * boseein # freqs in Phonopy is THz, so need to convert to Hz
        sigma = (np.sum(var)/nqpts)**0.5 # Square root of variance, have to do normalization over the number of q points 

    return var, sigma

        

    
 