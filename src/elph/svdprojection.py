import numpy as np
import yaml
#from qutip import * # quantum toolbox in python


def svd_projection(num_modes, nqpts, threshold=1e-9, matrix='epc'):
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
        Select the matrix to be used for SVD, can be 'epc' or 'variance'.

    Returns
    -------
    svd_epcA: ndarray
        SVD of electron-phonon coupling matrix for dimer A.
    svd_epcB: ndarray
        SVD of electron-phonon coupling matrix for dimer B.
    svd_epcC: ndarray
        SVD of electron-phonon coupling matrix for dimer C. 
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
    if matrix == 'epc':
        cp = np.load('svd_ep_coupling.npz') # Load epc
        epcA = cp['A'][0:num_modes*nqpts]
        epcB = cp['B'][0:num_modes*nqpts]
        epcC = cp['C'][0:num_modes*nqpts]
        epc = epcA+epcB+epcC
        print(f"EPC shape is {epc.shape}")
        U, S, Vh = np.linalg.svd(epc, full_matrices=True)  # Reduced SVD (if full_matrices is False): U is rotational orthogonal matrix; 
        # S is singular vectors (return singular value in 1D array); Vh is rotational orthogonal matrix
    
    if matrix == 'variance':
        variance = np.load('variance_for_svd.npz') # Load variance
        varA = variance['A'][0:num_modes*nqpts]
        varB = variance['B'][0:num_modes*nqpts]
        varC = variance['C'][0:num_modes*nqpts]
        var = varA + varB + varC
        print(f"Variance shape is {var.shape}")
        U, S, Vh = np.linalg.svd(var, full_matrices=True)
    
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


    if matrix == 'epc':
        svd_epcA1 = epcA[:,0] @ coeff_sys[0,:] + epcA[:,1] @ coeff_sys[0,:] + epcA[:,2] @ coeff_sys[0,:] 
        svd_epcA2 = epcA[:,0] @ coeff_sys[1,:] + epcA[:,1] @ coeff_sys[1,:] + epcA[:,2] @ coeff_sys[1,:] 
        svd_epcA3 = epcA[:,0] @ coeff_sys[2,:] + epcA[:,1] @ coeff_sys[2,:] + epcA[:,2] @ coeff_sys[2,:] 
        svd_epcA = np.array([svd_epcA1, svd_epcA2, svd_epcA3])
        svd_epcB1 = epcB[:,0] @ coeff_sys[0,:] + epcB[:,1] @ coeff_sys[0,:] + epcB[:,2] @ coeff_sys[0,:]
        svd_epcB2 = epcB[:,0] @ coeff_sys[1,:] + epcB[:,1] @ coeff_sys[1,:] + epcB[:,2] @ coeff_sys[1,:]
        svd_epcB3 = epcB[:,0] @ coeff_sys[2,:] + epcB[:,1] @ coeff_sys[2,:] + epcB[:,2] @ coeff_sys[2,:]
        svd_epcB = np.array([svd_epcB1, svd_epcB2, svd_epcB3])
        svd_epcC1 = epcC[:,0] @ coeff_sys[0,:] + epcC[:,1] @ coeff_sys[0,:] + epcC[:,2] @ coeff_sys[0,:]
        svd_epcC2 = epcC[:,0] @ coeff_sys[1,:] + epcC[:,1] @ coeff_sys[1,:] + epcC[:,2] @ coeff_sys[1,:]
        svd_epcC3 = epcC[:,0] @ coeff_sys[2,:] + epcC[:,1] @ coeff_sys[2,:] + epcC[:,2] @ coeff_sys[2,:]
        svd_epcC = np.array([svd_epcC1, svd_epcC2, svd_epcC3])

        return svd_epcA, svd_epcB, svd_epcC, f_sys, f_bath, coeff_sys, coeff_bath
    
    if matrix == 'var':
        svd_varA1 = varA[:,0] @ coeff_sys[0,:] + varA[:,1] @ coeff_sys[0,:] + varA[:,2] @ coeff_sys[0,:]
        svd_varA2 = varA[:,0] @ coeff_sys[1,:] + varA[:,1] @ coeff_sys[1,:] + varA[:,2] @ coeff_sys[1,:]
        svd_varA3 = varA[:,0] @ coeff_sys[2,:] + varA[:,1] @ coeff_sys[2,:] + varA[:,2] @ coeff_sys[2,:]
        svd_varA = np.array([svd_varA1, svd_varA2, svd_varA3])
        svd_varB1 = varB[:,0] @ coeff_sys[0,:] + varB[:,1] @ coeff_sys[0,:] + varB[:,2] @ coeff_sys[0,:]
        svd_varB2 = varB[:,0] @ coeff_sys[1,:] + varB[:,1] @ coeff_sys[1,:] + varB[:,2] @ coeff_sys[1,:]
        svd_varB3 = varB[:,0] @ coeff_sys[2,:] + varB[:,1] @ coeff_sys[2,:] + varB[:,2] @ coeff_sys[2,:]
        svd_varB = np.array([svd_varB1, svd_varB2, svd_varB3])
        svd_varC1 = varC[:,0] @ coeff_sys[0,:] + varC[:,1] @ coeff_sys[0,:] + varC[:,2] @ coeff_sys[0,:]
        svd_varC2 = varC[:,0] @ coeff_sys[1,:] + varC[:,1] @ coeff_sys[1,:] + varC[:,2] @ coeff_sys[1,:]
        svd_varC3 = varC[:,0] @ coeff_sys[2,:] + varC[:,1] @ coeff_sys[2,:] + varC[:,2] @ coeff_sys[2,:]
        svd_varC = np.array([svd_varC1, svd_varC2, svd_varC3])

        return svd_varA, svd_varB, svd_varC, f_sys, f_bath, coeff_sys, coeff_bath
    
    