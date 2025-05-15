import numpy as np
import os
import elph.elphtool as ep
#from qutip import * # quantum toolbox in python

def svd_projection(freq_array, svd_array, threshold=1e-8):
    """
    Singular value decomposition (SVD) of the electron-phonon coupling matrix
    and projection of the phonon modes into system and bath modes.

    Args:
    freq_array : ndarray
        Frequencies of the phonon modes.
    svd_array : ndarray
        Electron-phonon coupling matrix.
    threshold : float
        Threshold for singular values to be considered zero.
    ----------------------------------------------------------
    Returns:
    S: ndarray
        Singular values of the SVD.
    f_sys : ndarray
        Frequencies of system phonon modes.
    f_bath : ndarray
        Frequencies of bath phonon modes.
    coeff_sys : ndarray
        Coefficients of system phonon modes.
    coeff_bath : ndarray
        Coefficients of bath phonon modes.
    """
    print(" Running SVD projection ")
    U, S, Vh = np.linalg.svd(svd_array, full_matrices=True)
    print(f'Singular values are {S}') 
    print(f"Shape of left orthogonal matrix {U.shape}")
    print(f"Shape of singular vectors {S.shape}")
    print(f"Shape of right orthogonal matrix {Vh.shape}")

    Snonzero = np.where(S > threshold)[0] 
    print('Indices of non-zero singular values=', Snonzero)

    P = U[:, Snonzero] @ U[:, Snonzero].T # Projection operator
    I = np.eye(svd_array.shape[0]) # Identity operator
    Q = I - P # Complement operator
    print(f"Shape of projection operator {P.shape}")
    print(f"Projection operator test: ||P^2 - P|| = {np.linalg.norm(P @ P - P)}")
    print(f"Complement operator test: ||Q^2 - Q|| = {np.linalg.norm(Q @ Q - Q)}")

    hessian = np.diag(freq_array**2) # Define Hessian matrix
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
    print(f"Number of total phonon modes are {len(Pnonzero) + len(Qnonzero)}")
    print(f"Number of system phonon modes are {len(Pnonzero)}")
    print(f"Number of bath phonon modes are {len(Qnonzero)}")
    print(f"System phonon modes (LCAO of orginal phonon modes) frequencies are {f_sys} cm-1")
    print(f"Shape of system phonon modes coefficient {coeff_sys.shape}")
    print(f"Shape of bath phonon modes coefficient {coeff_bath.shape}")
      
    return S, f_sys, f_bath, coeff_sys, coeff_bath
    
    
    
    