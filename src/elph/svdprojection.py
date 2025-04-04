import numpy as np
from scipy.constants import hbar, k
import yaml
#from qutip import * # quantum toolbox in python


def svd_projection(num_modes, nqpts, matrix, threshold=1e-9, temp=298):
    """
    Singular value decomposition (SVD) of the electron-phonon coupling matrix
    and projection of the phonon modes into system and bath modes.

    Parameters
    ----------
    num_modes : int
        Number of phonon modes.
    nqpts : int
        Number of q-points.
    threshold : float
        Threshold for singular values to be considered zero.
    matrix : str
        Select the matrix to be used for SVD, can be 'epc' or 'epcbe'.
    temp: float
        Temperature in Kelvin. (Default is 298 K)

    Returns
    -------
    svd_epc: ndarray
        SVD of electron-phonon coupling matrix
    f_sys : ndarray
        Frequencies of system phonon modes.
    f_bath : ndarray
        Frequencies of bath phonon modes.
    coeff_sys : ndarray
        Coefficients of system phonon modes.
    coeff_bath : ndarray
        Coefficients of bath phonon modes.
    """
    freq = np.loadtxt('frequencies.txt')[0:num_modes*nqpts]
    print(f"Frequency shape is {freq.shape}")

    # Numpy SVD 
    cp = np.load('epc_for_svd.npz') # Load epc
    epcA = cp['A'][0:num_modes*nqpts]
    epcB = cp['B'][0:num_modes*nqpts]
    epcC = cp['C'][0:num_modes*nqpts]
    
    freqs = np.tile(freq[:, np.newaxis], (1, 3))
    b_e = 1 / np.tanh((hbar*freqs*1e12)/(2*k*temp))
    b_e = b_e[0:num_modes*nqpts]
    
    if matrix == 'epc':
        epc = epcA + epcB + epcC
        print(f"EPC shape is {epc.shape}")
        U, S, Vh = np.linalg.svd(epc, full_matrices=True)  # Reduced SVD (if full_matrices is False): U is rotational orthogonal matrix; 
        # S is singular vectors (return singular value in 1D array); Vh is rotational orthogonal matrix
    
    if matrix == 'epcbe':
        epc = epcA*b_e + epcB*b_e + epcC*b_e
        print(f"EPC shape is {epc.shape}")
        U, S, Vh = np.linalg.svd(epc, full_matrices=True)
    
    print(f'Singular values are {S}') 
    print(f"Shape of left orthogonal matrix {U.shape}")
    print(f"Shape of singular vectors {S.shape}")
    print(f"Shape of right orthogonal matrix {Vh.shape}")

    Snonzero = np.where(S > threshold)[0] 
    print('Indices of non-zero singular values=', Snonzero)

    P = U[:, Snonzero] @ U[:, Snonzero].T # Projection operator
    I = np.eye(num_modes * nqpts)
    Q = I - P # Complement operator
    print(f"Shape of projection operator {P.shape}")
    print(f"Projection operator test: ||P^2 - P|| = {np.linalg.norm(P @ P - P)}")
    print(f"Complement operator test: ||Q^2 - Q|| = {np.linalg.norm(Q @ Q - Q)}")

    hessian = np.diag(freq**2) # Define Hessian matrix
    print(f"Shape of Hessian matrix {hessian.shape}")
    hessian_sys = P @ hessian @ P # system mode
    hessian_bath = Q @ hessian @ Q # bath mode
    eigval_sys, eigvecs_sys = np.linalg.eigh(hessian_sys)
    eigval_bath, eigvecs_bath = np.linalg.eigh(hessian_bath)

    Pnonzero = np.where(eigval_sys > threshold)[0]
    Qnonzero = np.where(eigval_bath > threshold)[0]

    coeff_sys = eigvecs_sys[Pnonzero] # Eigenvectors for each phonon modes to construct the 3 system phonon modes (coefficient)
    coeff_bath = eigvecs_bath[Qnonzero] # Eigenvectors for each phonon modes to construct the remaining phonon modes
    f_sys = np.sqrt(eigval_sys[Pnonzero]) # Frequency of system phonon modes
    f_bath = np.sqrt(eigval_bath[Qnonzero]) # Frequency of bath phonon modes
    thztocm_1 = 33.356
    f_sys *= thztocm_1
    f_bath *= thztocm_1
    print(f"Number of total phonon modes are {num_modes*nqpts}")
    print(f"Number of system phonon modes are {len(Pnonzero)}")
    print(f"Number of bath phonon modes are {len(Qnonzero)}")
    print(f"System phonon modes (LCAO of orginal phonon modes) frequencies are {f_sys} cm-1")
    print(f"Shape of system phonon modes coefficient {coeff_sys.shape}")
    print(f"Shape of bath phonon modes coefficient {coeff_bath.shape}")

    # Explanation: svd_epcA1 equals to epcA in x axis * coeff_sysmode1 + epcA in y axis * coeff_sysmode1 + epcA in z axis * coeff_sysmode1
    
    vecs = np.vstack((coeff_sys, coeff_bath))
    svd_epcx = epc[:,0] @ vecs
    svd_epcy = epc[:,1] @ vecs
    svd_epcz = epc[:,2] @ vecs
    svd_epc = np.array([svd_epcx, svd_epcy, svd_epcz])
      
    return svd_epc, f_sys, f_bath, coeff_sys, coeff_bath
    
    
    
    