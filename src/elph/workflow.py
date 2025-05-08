import ase
import glob
import json
import numpy as np
import os
import sys
import yaml
from scipy.constants import hbar
import elph.utils as ut
import elph.elph as ep
import elph.svdprojection as svd
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

def run_j0(mol_list, supercell_matrix, basis):
    """ Main function for running Gaussian and Catnip to get transfer integral J_0
    Args:
    mol_list (list): After visualization, you need to specify 3 molecules. The order is: 1st; 2nd; 3rd
    supercell_matrix (tuple): The supercell matrix
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
    main_path = os.getcwd() # Main directory which contain all subfolders
    j0_file = glob.glob(os.path.join(main_path, 'j0_eff.json'))
    xyz_file = glob.glob(os.path.join(main_path, '1', 'monomer_1.xyz'))
    if not j0_file: # If j_0.json is not exists, run the following simulation
        try:
            geometry = getGeometry(main_path) # Get the geometry file

        except FileNotFoundError:
                ut.print_error("CIF file not found in the current directory. Exiting.") 
                sys.exit(1)  # Exit the script with an error

        if not xyz_file:
            ep.unwrap_molecule_dimer(geometry, supercell_matrix, mol_list[0], mol_list[1], mol_list[2]) # Unwrap the crystal to get single molecule and dimers
        
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

def run_lambda(basis):
    """ Run onsite energy calculation using normal mode analysis to get reorganization energy &
    local EP coupling
    Args:
    basis (str): Gaussian basis sets
    -----------------------------------------------
    Return:
    local.json file which saves onsite energy & reorganization energy (units: eV)
    """
    main_path = os.getcwd()
    orginal_atoms = ase.io.read(f'{main_path}/1/monomer_1.xyz') # Read the structure file

    if not os.path.exists(f'{main_path}/local/local_epc.json'):
        print(' Running local EPC calculation ... ')
        os.mkdir(os.path.join(main_path, 'local'))
        os.chdir(os.path.join(main_path, 'local'))
    
        ep.gaussian_opt(atoms=orginal_atoms, bset=basis, label='neutral', ncharge=0)
        ep.gaussian_opt(atoms=orginal_atoms, bset=basis, label='cation', ncharge=1)

        ep.hr_factor()
        eng_n, freq_n, huangrhys_n, reorg_n, gii_n, gii_n_cart = ep.parse_log('neutral.log', 'hr_neutral.log')
        eng_c, freq_c, huangrhys_c, reorg_c, gii_c, gii_c_cart = ep.parse_log('cation.log', 'hr_cation.log')
        # 4-point method
        # Calculate onsite energy for relaxed neutral molecule (Rn)
        #os.mkdir(os.path.join(main_path, 'reorgE', 'neutral_opt'))
        #os.chdir(os.path.join(main_path, 'reorgE', 'neutral_opt'))
        #ep.onsite_eng(atoms=orginal_atoms, bset=basis, opt=True)
        #neutral_rn_eng = ep.parse_log('mo.log')
        #os.chdir(os.pardir)
        # Calculate onsite energy for relax charged molecule (Rc)
        #os.mkdir(os.path.join(main_path, 'reorgE', 'charged_opt'))
        #os.chdir(os.path.join(main_path, 'reorgE', 'charged_opt'))
        #ep.onsite_eng(atoms=orginal_atoms, bset=basis, ncharge=1, opt=True)
        #charged_rc_eng = ep.parse_log('mo_opt.log')
        #os.chdir(os.pardir)

        # Calculate onsite energy for charged molecule at Rn
        #os.mkdir(os.path.join(main_path, 'reorgE', 'charged_at_Rn'))
        #os.chdir(os.path.join(main_path, 'reorgE', 'charged_at_Rn'))
        #cmol_rn = ase.io.read(os.path.join(main_path, 'onsiteE', 'charged_opt', 'opt.xyz')) # Read the structure file
        #ep.onsite_eng(atoms=cmol_rn, bset=basis, ncharge=1)
        #charged_rn_eng = ep.parse_log('mo_opt.log')
        #os.chdir(os.pardir)
        # Calcualte onsite energy for neutral molecule at Rc
        #os.mkdir(os.path.join(main_path, 'reorgE', 'neutral_at_Rc'))
        #os.chdir(os.path.join(main_path, 'reorgE', 'neutral_at_Rc'))
        #mol_rc = ase.io.read(os.path.join(main_path, 'reorgE', 'neutral_opt', 'opt.xyz')) # Read the structure file
        #ep.onsite_eng(atoms=mol_rc, bset=basis)
        #neutral_rc_eng = ep.parse_log('mo_opt.log')
        #os.chdir(os.pardir) 

        #reorg_eng = charged_rn_eng + neutral_rc_eng - neutral_rn_eng - charged_rc_eng

        #onsiteE = {'neutral_rc': neutral_rc_eng,
        #           'neutral_rn': neutral_rn_eng,
        #           'charged_rn': charged_rn_eng,
        #           'charged_rc': charged_rc_eng,
        #           'reorg': reorg_eng
        #        }

        local = {'eng_n': eng_n,
                 'eng_c': eng_c,
                 'reorg_eng_n': sum(reorg_n),
                 'reorg_eng_c': sum(reorg_c),
                }
        
        with open('local_epc.json', 'w', encoding='utf-8') as f:
            json.dump(local, f, ensure_ascii=False, indent=4)

        data = {'gii_n': gii_n,
                'gii_c': gii_c,
                'gii_n_cart': gii_n_cart,
                'gii_c_cart': gii_c_cart,
                'freq_n': freq_n,
                'freq_c': freq_c,
                'huangrhys_n': huangrhys_n,
                'huangrhys_c': huangrhys_c,
                'reorg_eng_n': reorg_n,
                'reorg_eng_c': reorg_c
               }
        np.savez_compressed('local_epc.npz', **data) 
        os.chdir(main_path)

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
            ep.mol_orbital(bset=basis) # Run Gaussian to get molecular orbitals
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
            #dimers = ase.io.read(f'{main_path}/{key}/dimer_{key}.xyz')
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
        
def run_matrix(mesh, sc):
    """ Calculate electron phonon coupling matrix and then connect with each phonon mode (from Phonopy)
    Dependency: 
    phonon(mesh), this function is to get phonon modes
    ------------------------------------------------------------
    Args:
    mesh (list): Need define a mesh grid. (Defaults to [8,8,8])
    sc (list): The supercell matrix. (Defaults to [2,2,2])
    """
    ####### Calculate J_ij matrix #########
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
    mod, freqs, nqpts = ep.phonon(natoms,mesh) # Run phonopy modulation to create eigendisplacements list 
    displacement = np.tile(mod, (1, sc[0]*sc[1]*sc[2], 1)) # The shape of displacement is [ phonon modes(number of q points * number of atoms in unitcell * 3), number of atoms in supercell, 3 (x,y,z) ]
    
    # the shape of displacement_list is [ phonon modes(number of q points * number of atoms in unitcell * 3), number of atoms in supercell, 3 (x,y,z) ]
    displacement_A = displacement[:,mapping_A,:]  # Assign the corresponding displacement
    displacement_B = displacement[:,mapping_B,:]  # Assign the corresponding displacement
    displacement_C = displacement[:,mapping_C,:]  # Assign the corresponding displacement
         
    #epcA = np.einsum('ij,kij->k', matrix_A, displacement_A)  # i is number of atoms, j is 3 (x,y,z); k is the index of phonon modes
    #epcB = np.einsum('ij,kij->k', matrix_B, displacement_B)  # We get coefficient g_IJ here           
    #epcC = np.einsum('ij,kij->k', matrix_C, displacement_C)  

    gii_n_squared = np.load('local/local_epc.npz')['gii_n'] # Load the gii from local_epc.json file
    gii_c_squared = np.load('local/local_epc.npz')['gii_c'] 
    gii_n_cart_squared = np.load('local/local_epc.npz')['gii_n_cart'] # cart means the cartesian coordinates
    gii_c_cart_squared = np.load('local/local_epc.npz')['gii_c_cart']
    epcL_squared = np.concatenate((gii_n_squared, gii_c_squared), axis=0)
    epcL_cart_squared = np.concatenate((gii_n_cart_squared, gii_c_cart_squared), axis=0)

    epcA_cart = np.einsum('ij,kij->kj', matrix_A, displacement_A)  # i is number of atoms, j is 3 (x,y,z; cartesian); k is the index of phonon modes
    epcB_cart = np.einsum('ij,kij->kj', matrix_B, displacement_B)  # The shape here is k,3           
    epcC_cart = np.einsum('ij,kij->kj', matrix_C, displacement_C) 

    epcA_cart_squared = epcA_cart**2
    epcB_cart_squared = epcB_cart**2
    epcC_cart_squared = epcC_cart**2

    epcA_squared = np.sum(epcA_cart_squared, axis=1)
    epcB_squared = np.sum(epcB_cart_squared, axis=1)
    epcC_squared = np.sum(epcC_cart_squared, axis=1)
                                      
    epc_squared = {'L':epcL_squared,
                   'A':epcA_squared,
                   'B':epcB_squared,
                   'C':epcC_squared}
    
    epc_cart_squared = {'L':epcL_cart_squared,
                        'A':epcA_cart_squared, # This is the epc matrix for running SVD projection
                        'B':epcB_cart_squared,
                        'C':epcC_cart_squared}
    
    # Save the electron-phonon coupling matrix asp a numpy .npz file.
    np.savez_compressed('epc' + '.npz', **epc_squared)
    np.savez_compressed('epc_cart' + '.npz', **epc_cart_squared) # Save the svd electron-phonon coupling matrix as a numpy .npz file.

    # Calculate the variance of the electron-phonon coupling matrix
    variL, sigmaL = ep.variance(freqs, epcL_squared, nqpts, 298) # Variance for local el-ph coupling matrix
    variA, sigmaA = ep.variance(freqs, epcA_squared, nqpts, 298) # Variance for dimer A
    variB, sigmaB = ep.variance(freqs, epcB_squared, nqpts, 298) # Variance for dimer B
    variC, sigmaC = ep.variance(freqs, epcC_squared, nqpts, 298) # Variance for dimer C

    variance = {'vL':variL,
                'sL':sigmaL,
                'vA':variA,
                'sA':sigmaA,
                'vB':variB,
                'sB':sigmaB,
                'vC':variC,
                'sC':sigmaC} 

    np.savez_compressed('variance' + '.npz', **variance)

def run_tlt_mobility(filename="mobility.json", output="tlt_mobility"):
    """
    Run TLT mobility simulation
    Args:
    filename (str): The mobility.json file to run TLT mobility simulation
    output (str): The output name for TLT mobility simulation (Defaults to tlt_mobility)
    """
    print(" Running TLT mobility simulation using parameters in mobility.json ... ")

    mobility = Mobility(mob_file=filename)
    avglx2, avgly2, mobilityx, mobilityy, mobility_average = mobility.tlt_mobility()

    mobility = {'Localization length in X': round(avglx2, 3),
                'Localization length in Y': round(avgly2, 3),
                'Mobility on X direction (cm^/Vs)': round(mobilityx, 3), 
                'Mobility on Y direction (cm^/Vs)': round(mobilityy, 3),
                'Average mobility': round(mobility_average, 3)
    }
    
    with open(f'{output}.json', 'w', encoding='utf-8') as f:
        json.dump(mobility, f, ensure_ascii=False, indent=4)

def run_svd_projection(matrix, nqpts):
    """
    Run SVD projection using numpy
    Args:
    matrix (str): The matrix to run SVD projection (epc or var)
    nqpts (int): The number of q points to run SVD projection
    """
    print(" Running phonon modes projections using SVD ... ")

    main_path = os.getcwd()
    atoms = ase.io.read(getGeometry(main_path)) # Read the structure file
    natoms = len(atoms)
    nmodes = 3 * natoms
    if matrix == 'epc':
        svd_epc, f_sys, f_bath, coeff_sys, coeff_bath = svd.svd_projection(num_modes=nmodes, nqpts=nqpts, matrix='epc')

        result = {'epc': svd_epc,
                   'freq_sys':f_sys,
                   'freq_bath':f_bath,
                   'coeff_sys':coeff_sys,
                   'coeff_bath':coeff_bath}
    
        np.savez_compressed('svd_epc_result' + '.npz', **result) 

    if matrix == 'epcbe':
        svd_epc, f_sys, f_bath, coeff_sys, coeff_bath = svd.svd_projection(num_modes=nmodes, nqpts=nqpts, matrix='epcbe')

        result = {'epc': svd_epc,
                   'freq_sys':f_sys,
                   'freq_bath':f_bath,
                   'coeff_sys':coeff_sys,
                   'coeff_bath':coeff_bath}
    
        np.savez_compressed('svd_epcbe_result' + '.npz', **result)
    
   
    