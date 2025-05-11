import numpy as np
import os
import elph.elphtool as ep
#from qutip import * # quantum toolbox in python
cm_1tothz = 3e-2  

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
    main_path = main_path = os.getcwd()
    freq_l_n = np.load(os.path.join(main_path, 'local', 'local_epc.npz'))['freq_n']
    freq_l_n *= cm_1tothz # Convert to Hz
    freq_l_c = np.load(os.path.join(main_path, 'local', 'local_epc.npz'))['freq_c']
    freq_l_c *= cm_1tothz # Convert to Hz
    freq_l = np.concatenate((freq_l_n, freq_l_c), axis=0)
    freq_nl = np.loadtxt('frequencies.txt')[0:num_modes*nqpts]
    freq_tot = np.concatenate((freq_l, freq_nl), axis=0)

    # Numpy SVD 
    epc_cart = np.load('epc_cart.npz') # Load epc
    epcL_cart = epc_cart['L']
    epcA_cart = epc_cart['A'][0:num_modes*nqpts]
    epcB_cart = epc_cart['B'][0:num_modes*nqpts]
    epcC_cart = epc_cart['C'][0:num_modes*nqpts]
    
    freqs_l = np.tile(freq_l[:, np.newaxis], (1, 3))
    freqs_nl = np.tile(freq_nl[:, np.newaxis], (1, 3))

    print(" Running SVD projection ")

    if matrix == 'epc':
        epcnl_cart = epcA_cart + epcB_cart + epcC_cart
        epc = np.concatenate((epcL_cart, epcnl_cart), axis=0)
        print(f"EPC shape is {epc.shape}")
        U, S, Vh = np.linalg.svd(epc, full_matrices=True)  # Reduced SVD (if full_matrices is False): U is rotational orthogonal matrix; 
        # S is singular vectors (return singular value in 1D array); Vh is rotational orthogonal matrix
    
    if matrix == 'epcbe':
        boseein_l, _, _ = ep.variance(freqs_l, epcL_cart, 1, temp, unit='cm-1')
        boseein_nl, _, _ = ep.variance(freqs_nl, (epcA_cart + epcB_cart + epcC_cart), nqpts, temp, unit='THz')
        boseein_nl = boseein_nl[0:num_modes*nqpts]
        epcnl_cart = 0.5*(epcA_cart**2)*boseein_nl + 0.5*(epcB_cart**2)*boseein_nl + 0.5*(epcC_cart**2)*boseein_nl
        epc = np.concatenate((0.5*(epcL_cart**2)*boseein_l, epcnl_cart), axis=0)
        print(f"EPC shape is {epc.shape}")
        U, S, Vh = np.linalg.svd(epc, full_matrices=True)
    
    print(f'Singular values are {S}') 
    print(f"Shape of left orthogonal matrix {U.shape}")
    print(f"Shape of singular vectors {S.shape}")
    print(f"Shape of right orthogonal matrix {Vh.shape}")

    Snonzero = np.where(S > threshold)[0] 
    print('Indices of non-zero singular values=', Snonzero)

    P = U[:, Snonzero] @ U[:, Snonzero].T # Projection operator
    I = np.eye(epc.shape[0]) # Identity operator
    Q = I - P # Complement operator
    print(f"Shape of projection operator {P.shape}")
    print(f"Projection operator test: ||P^2 - P|| = {np.linalg.norm(P @ P - P)}")
    print(f"Complement operator test: ||Q^2 - Q|| = {np.linalg.norm(Q @ Q - Q)}")

    hessian = np.diag(freq_tot**2) # Define Hessian matrix
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

    # Explanation: svd_epcA1 equals to epcA in x axis * coeff_sysmode1 + epcA in y axis * coeff_sysmode1 + epcA in z axis * coeff_sysmode1
    
    vecs = np.vstack((coeff_sys, coeff_bath))
    svd_epcx = epc[:,0] @ vecs
    svd_epcy = epc[:,1] @ vecs
    svd_epcz = epc[:,2] @ vecs
    svd_epc = np.array([svd_epcx, svd_epcy, svd_epcz])
      
    return svd_epc, f_sys, f_bath, coeff_sys, coeff_bath
    
    
    
    