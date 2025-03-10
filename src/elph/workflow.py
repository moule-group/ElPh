import ase
import glob
import json
import numpy as np
import os
import sys
import yaml
import elph.utils as ut
import elph.elph as ep
from elph.mobility import Mobility

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
    file = glob.glob(path + "/*.cif") +  glob.glob(path + "/*.vasp") + glob.glob(path + "/CONTCAR") + glob.glob(path + "/*.xyz")
    if len(file) == 0:
        raise FileNotFoundError
    
    return file[0]

def run_j0(mol_list, basis):
    """ Main function for running Gaussian and Catnip to get transfer integral J_0
    Args:
    mol_list (list): After visualization, you need to specify 3 molecules. The order is: 1st; 2nd; 3rd
    basis (str): Gaussian basis sets
    ########################################################################
    Here is the example 2D figure below, 1 and 2 are molecules. You have to reference the numbering of Ref; and 3 molecules with * sign.
    ----------------------------
          *      *     *
          
      #      2#      3#
                   
          *      1*    *
                        
      #      #       #
    ---------------------------- 
  
    ########################################################################
    Return:
    j_A, j_B, j_C as j0.json and j0_eff.json file
    """    
    main_path = os.getcwd()
    json_file = glob.glob(os.path.join(main_path, 'j0_eff.json'))
    xyz_file = glob.glob(os.path.join(main_path,'1', 'monomer_1.xyz'))
    main_path = os.getcwd() # Main directory which contain all subfolders
    if not json_file: # If j_0.json is not exists, run the following simulation
        try:
            geometry = getGeometry(main_path) # Get the geometry file

        except FileNotFoundError:
                ut.print_error("CIF file not found in the current directory. Exiting.") 
                sys.exit(1)  # Exit the script with an error

        if not xyz_file:
            ep.unwrap_molecule_dimer(geometry, mol_list[0], mol_list[1], mol_list[2]) # Unwrap the crystal to get single molecule and dimers
        
        path_list = ['./1','./2','./3','./A','./B','./C']
        for path in path_list:
            os.chdir(path)
            ep.mol_orbital(bset=basis) # Run Gaussian to get molecular orbitals
            os.chdir(main_path)

        # Calculate J 
        jA_eff, jA = ep.run_catnip('./1/1.pun', './2/2.pun', './A/A.pun', './1/mo.log', './2/mo.log', './A/mo.log')
        jB_eff, jB = ep.run_catnip('./1/1.pun', './3/3.pun', './B/B.pun', './1/mo.log', './3/mo.log', './B/mo.log')
        jC_eff, jC = ep.run_catnip('./2/2.pun', './3/3.pun', './C/C.pun', './2/mo.log', './3/mo.log', './C/mo.log')

        print(f' Done calculation on J_A = {jA_eff} eV ')
        print(f' Done calculation on J_B = {jB_eff} eV ')
        print(f' Done calculation on J_C = {jC_eff} eV ')
    
        j0_eff = {'A':f'{jA_eff}', 
                  'B':f'{jB_eff}',
                  'C':f'{jC_eff}'
                 }
        
        j0 = {'A':f'{jA}', 
              'B':f'{jB}',
              'C':f'{jC}'
             }
    
        with open('j0_eff.json', 'w', encoding='utf-8') as f1:
            json.dump(j0_eff, f1, ensure_ascii=False, indent=4)
            
        with open('j0.json', 'w', encoding='utf-8') as f2:
            json.dump(j0, f2, ensure_ascii=False, indent=4)

def run_disp_j(basis):
    """ Main function for running Gaussian and Catnip to get transfer integral J for displaced molecules and dimers
    Return:
    j_list (list): The transfer integral list for displaced molecules and dimers!
    """
    main_path = os.getcwd()
    if not os.path.exists(f"{main_path}/C/displacements"):
        print(' Creating displaced molecules and dimers ... ')
        ep.create_displacement()
        
    print(" Displacement folders are finished! ")
    
    path_list = ['./1/displacements','./2/displacements','./3/displacements',
                 './A/displacements','./B/displacements','./C/displacements']
    
    for path in path_list:
        root, dirs, files = next(os.walk(path)) # Get the directory under each displacements/ folder
        os.chdir(path)  
        for d in dirs:
            os.chdir(d)
            ep.mol_orbital(opt=0, bset=basis) # Run Gaussian to get molecular orbitals
            os.chdir(os.pardir)
        os.chdir(main_path)
        
    #### Calculate J ####
    
    if not os.path.exists(f"{main_path}/C_disp_J.npz"): # Check whether it is necessary to run Catnip
        
        dimer_dict = {'A':[1,2],
                      'B':[1,3],
                      'C':[2,3]} # The dict for the dimer A, B, C
      
        for key, values in dimer_dict.items():
            j_list = [] # Reset j_list for each dimer
            mol_1 = values[0] # molecule 1
            mol_2 = values[1] # molecule 2
        
            molecules = ase.io.read(f'{main_path}/{mol_1}/monomer_{mol_1}.xyz')
            dimers = ase.io.read(f'{main_path}/{key}/dimer_{key}.xyz')
            offset = len(molecules) # This is because the for loop below only loop through molecule, and number of atoms in dimer is 2 times of a molecule
        
            for na, vec, sign in ep.get_displacement(molecules):  
                j_1, _ = ep.run_catnip(f'./{mol_1}/displacements/disp_{na}_{vec}_{sign}/disp_{na}_{vec}_{sign}.pun', 
                                    f'./{mol_2}/{mol_2}.pun', 
                                    f'./{key}/displacements/disp_{na}_{vec}_{sign}/disp_{na}_{vec}_{sign}.pun', 
                                    f'./{mol_1}/displacements/disp_{na}_{vec}_{sign}/mo.log', 
                                    f'./{mol_2}/mo.log', 
                                    f'./{key}/displacements/disp_{na}_{vec}_{sign}/mo.log')
                j_list.append(j_1)
        
                j_2, _ = ep.run_catnip(f'./{mol_1}/{mol_1}.pun', 
                                    f'./{mol_2}/displacements/disp_{na}_{vec}_{sign}/disp_{na}_{vec}_{sign}.pun', 
                                    f'./{key}/displacements/disp_{na+offset}_{vec}_{sign}/disp_{na+offset}_{vec}_{sign}.pun', 
                                    f'./{mol_1}/mo.log', 
                                    f'./{mol_2}/displacements/disp_{na}_{vec}_{sign}/mo.log', 
                                    f'./{key}/displacements/disp_{na+offset}_{vec}_{sign}/mo.log')
                j_list.append(j_2)
        
            data = {'J_ij': j_list} 
            np.savez_compressed(key + '_disp_J.npz', **data)
            print(f" Successfully create {key}_disp_J.npz file which saves J_ij!!! ")
        
def run_matrix(mesh,sc):
    """ Calculate electron phonon coupling matrix and then connect with each phonon mode (from Phonopy)
    Dependency: 
    phonon(mesh), this function is to get phonon modes
    ########################################################
    Args:
    mesh (list): Need define a mesh grid. (Defaults to [8,8,8])
    sc (list): The supercell matrix. (Defaults to [2,2,2])
    """
    ####### Calculate J matrix #########
    jlist_A = np.load('A_disp_J.npz')['J_ij'] # - +
    jlist_B = np.load('B_disp_J.npz')['J_ij'] # - +
    jlist_C = np.load('C_disp_J.npz')['J_ij'] # - + 
 
    matrix_A = ep.get_deri_Jmatrix(jlist_A)
    matrix_B = ep.get_deri_Jmatrix(jlist_B)
    matrix_C = ep.get_deri_Jmatrix(jlist_C)
    
    ####### Connection with Phonon modes ########
    with open('monomer.json', 'r') as j:
        mapping = list(json.load(j).values()) # The numbering of atoms for dimers
        
    mapping_A = mapping[0] + mapping[1] # molecule 1 + molecule 2
    mapping_B = mapping[0] + mapping[2] # molecule 1 + molecule 3
    mapping_C = mapping[1] + mapping[2] # molecule 2 + molecule 3
    
    main_path = os.getcwd()
    atoms = ase.io.read(getGeometry(main_path)) # Read the structure file
    natoms = len(atoms) # Number of atoms in the supercell
    displacement, freqs = ep.phonon(natoms,mesh,sc) # Run phonopy modulation to create eigendisplacements list 
    # the shape of displacement_list is [ phonon modes(number of q point * number of atoms in unitcell * 3), number of atoms in supercell, 3 (x,y,z) ]
    
    displacement_A = np.zeros((displacement.shape[0],len(mapping_A),3))
    displacement_B = np.zeros((displacement.shape[0],len(mapping_B),3)) 
    displacement_C = np.zeros((displacement.shape[0],len(mapping_C),3))
    
    for i, a in enumerate(mapping_A): # i is the index of the dimer, a is the atom index of the supercell
        displacement_A[:, i, :] = displacement[:,a,:]  # Assign the corresponding displacement

    for i, b in enumerate(mapping_B):
        displacement_B[:, i, :] = displacement[:,b,:]  # Assign the corresponding displacement

    for i, c in enumerate(mapping_C):
         displacement_C[:, i, :] = displacement[:,c,:]  # Assign the corresponding displacement
    
    ep_couplingA = np.einsum('ij,kij->k',matrix_A, displacement_A)  # i is number of atoms, j is 3 (x,y,z); k is the index of phonon modes
    ep_couplingB = np.einsum('ij,kij->k',matrix_B, displacement_B)  # We get coefficient g_ij here           
    ep_couplingC = np.einsum('ij,kij->k',matrix_C, displacement_C)               
    
    ep_coupling = {'A':ep_couplingA,
                   'B':ep_couplingB,
                   'C':ep_couplingC}
    
    # Save the electron-phonon coupling matrix asp a numpy .npz file.
    np.savez_compressed('ep_coupling' + '.npz', **ep_coupling)

    variA = ep.variance(freqs, ep_couplingA, 298) # Variance for dimer A
    variB = ep.variance(freqs, ep_couplingB, 298) # Variance for dimer B
    variC = ep.variance(freqs, ep_couplingC, 298) # Variance for dimer C

    variance = {'A':variA,
                'B':variB,
                'C':variC}

    np.savez_compressed('variance' + '.npz', **variance)

def run_tlt_mobility(filename="mobility.json", output="tlt_mobility"):
    """
    Run TLT mobility simulation
    """
    print(" Running TLT mobility simulation using parameters in mobility.json ... ")

    mobility = Mobility(mob_file=filename)
    mobilityx, mobilityy, mobility_average = mobility.tlt_mobility()

    mobility = {'Mobility on X direction (cm^/Vs)': round(mobilityx, 3), 
                'Mobility on Y direction (cm^/Vs)': round(mobilityy, 3),
                'Average mobility ': round(mobility_average, 3)
    }
    
    with open(f'{output}.json', 'w', encoding='utf-8') as f:
        json.dump(mobility, f, ensure_ascii=False, indent=4)

   
    