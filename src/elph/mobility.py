# Charge carrier mobility using Transient Localization Theory (TLT)
import numpy as np
from scipy.constants import e, k


def hamiltonian(atoms, j_ij, sigma, nx, ny):
    """ Define the tight-binding Hamiltonian matrix for the charge carrier.
    H = H_el + H_ph + H_elph,l + H_elph,nl
      = (H_ii + H_ij) + H_ph + H_elph,l + H_elph,nl
    In TLT limit, because we are considering short timing, we neglect phonon term. We only need to include electronic 
    & e-p coupling term. (Onsite energy can be set to zero or not, depending on the accuracy level you want to achieve)
    Args:
    atoms: ase atoms object
    j_ij: Inter-molecular transfer integral (J_a, J_b, J_c)
    sigma: dynamic disorder (sigma_a, sigma_b, sigma_c)
    nx (int): supercell size in x-direction
    ny (int): supercell size in y-direction
    """
    nmol = len(atoms)
    H = np.zeros((nx*nmol,ny*nmol), dtype=complex)
    
    H = 

    return H

def local_length():
    """
    Calculate the localization length of the charge carrier.
    Diagonalization of the Hamiltonian matrix with static disorder
    """
    path_list = ['./1/displacements','./2/displacements','./3/displacements',
                 './A/displacements','./B/displacements','./C/displacements']
    
    for path in path_list:
        root, dirs, files = next(os.walk(path)) # Get the directory under each displacements/ folder
        os.chdir(path)  
        for d in dirs:
            os.chdir(d)
            os.chdir(os.pardir)
        os.chdir(main_path)
    
    
    
    partition = np.exp(1/k*temp) 

    H = # Hamiltonian matrix
    eigenvalues, eigenvectors = np.linalg.eig(H)


def relax_time(energy=5e-3):
    """
    Calculate the relaxation time of the charge carrier.
    energy (float): Energy of intermolecular phonon modes (Defaults to 5 meV)
    """
    tau = energy * hbar

    return tau

def mobility(temp):
    """
    TLT mobility calculation.
    """

    l_length = local_length()
    tau = relax_time() # relaxation time \tau




    mobility = (e/(k*temp)) * (l_length**2/(2*tau) )