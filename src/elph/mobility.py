# Charge carrier mobility using Transient Localization Theory (TLT) and Marcus theory with KMC
import numpy as np
import sys
import json
from scipy.spatial import cKDTree

hbar = 6.582e-16 # eV s
kb = 8.6173e-5 # eV K-1  
e = 1 # Charge on electron in eV/V

class Mobility():
    """
    Args:
    site (np.array): containing the COM of the in the crystal unit cell as rows, expressed in fractional coordinate. 
    n (int): indicating the number of times the crystal unit cell is repeated along the coordinate axis. 
    r (float): the cutoff radius for finding nearest neighbors in fractional length (1 means 1 translation)
    lattice (np.array): Unit cell lattice vectors
    j_ij (list): Inter-molecular transfer integral (ex: J_a, J_b, J_c)
    sigma_ii (float): local electronic phonon coupling (ex: sigma_ii)
    sigma_ij (list): nonlocal electronic phonon coupling (dynamic disorder) (ex: sigma_a, sigma_b, sigma_c)
    temp (float): Temperature in Kelvin (Defaults to 298)
    inverse_htau (float): Inverse of the scattering time (hbar/tau) units in eV (Defaults to 5e-3)
    is_hole (bool): If True, hole transport, otherwise electron transport (Defaults to True)
    realizations (int): Number of realizations for average calculation (Defaults to 250)
    mob_file (str): The json file containing the mobility parameters (Defaults to "mobility.json")
    """
    def __init__(self, site=None, n=None, r=1, lattice=None, nearest_vecs=None, Lambda=0.0, Epsilon=0.0, j_ij=None, sigma_ii=0.0, sigma_ij=None, temp=298.0, inverse_htau=5e-3, is_hole=True, realizations=250, 
                 mob_file="mobility.json"):
        
        if mob_file:
            with open(mob_file, "r") as file:
                config = json.load(file)
            
            self.site = np.array(config.get("site", site))
            self.n = config.get("n", n)
            self.r = config.get("r", r)
            self.lattice = np.array(config.get("lattice", lattice))
            self.nearest_vecs = config.get("nearest_vecs", nearest_vecs)
            self.Lambda = config.get("Lambda", Lambda)
            self.Epsilon = config.get("Epsilon", Epsilon)
            self.j_ij = config.get("j_ij", j_ij)
            self.sigma_ii = config.get("sigma_ii", sigma_ii)
            self.sigma_ij = config.get("sigma_ij", sigma_ij)
            self.temp = config.get("temp", temp)
            self.inverse_htau = config.get("inverse_htau", inverse_htau)
            self.is_hole = config.get("is_hole", is_hole)
            self.realizations = config.get("realizations", realizations)

        else:
            self.site = np.array(site)
            self.n = n
            self.r = r
            self.lattice = np.array(lattice)
            self.nearest_vecs = np.array(nearest_vecs)
            self.Lambda = Lambda
            self.Epsilon = Epsilon
            self.j_ij = j_ij
            self.sigma_ii = sigma_ii
            self.sigma_ij = sigma_ij
            self.temp = temp
            self.inverse_htau = inverse_htau
            self.is_hole = is_hole
            self.realizations = realizations

    def generate_lattice(self):
        '''
        This function generates a lattice of atoms to populate a simulation cell.
        Returns:
        positions (np.array): containing the positions of the atoms in the
        simulation cell as rows, expressed in units of the lattice parameter.
        '''
        n_in_cell = self.site.shape[0]
        positions = np.zeros((n_in_cell * self.n * self.n, 2))
    
        count = 0
        for a in range(self.n):
            for b in range(self.n):
                positions[count:(count + n_in_cell),] = self.site + [a, b]
                count += n_in_cell
        return positions

    def interaction(self, sites):
        """ 
        Decide type of J based on geometry of sites
        args:
        sites (np.array): The positions of the center for each molecule in the simulation cell
        nearest_distance (list): Specific interaction distances to consider
        --------------------------------------------------------------------
        Returns:
        interaction_matrix (np.array): The interaction matrix
        dist_vecs (np.array): The distance vectors array
        """
        N = len(sites) # number of molecules in supercell
        interaction_matrix = np.zeros((N, N), dtype=int)

        tree = cKDTree(sites, boxsize=self.n)
        neighbor_indices = tree.query_ball_tree(tree, 1)
        for i, neighbors in enumerate(neighbor_indices):
            neighbor_indices[i] = [j for j in neighbors if j != i] # remove self index

        for i in range(N):
            for j in neighbor_indices[i]:
                dist = (sites[j] - sites[i])
                values = self.check_neighbors(dist)
                if values > 0:
                    interaction_matrix[i][j] = values
                else:
                    interaction_matrix[i][j] = 0

        return interaction_matrix

    def hamiltonian(self, sites):
        """ Define the tight-binding Hamiltonian matrix for the charge carrier.
        H = H_el + H_ph + H_elph
        in original TLT: H_ph = 0, H_ii = 0; but we can add H_ii and H_elph,l
        H = (H_ii + H_elph,l) + H_ij + H_elph,nl
        ---------------------------------------------
        Return:
        H: Hamiltonian matrix
        """
        interaction_matrix = self.interaction(sites)
        Hij_matrix = np.copy(interaction_matrix).astype(float) # Transfer integral matrix (J_ij)
        sigmaij_matrix = np.copy(interaction_matrix).astype(float) # Dynamic disorder matrix (in TLT, we treat this as static disorder)

        # Diagonal (H_ii)
        diag_eng = np.random.normal(loc=self.Epsilon, scale=self.sigma_ii, size=len(sites))

        # Inter-molecular transfer integral matrix (H_ij)
        j1 = self.j_ij[0]
        j2 = self.j_ij[1]
        j3 = self.j_ij[2]

        Hij_matrix[Hij_matrix==1] = j1
        Hij_matrix[Hij_matrix==2] = j2
        Hij_matrix[Hij_matrix==3] = j3

        s1 = self.sigma_ij[0]
        s2 = self.sigma_ij[1]
        s3 = self.sigma_ij[2]

        sigmaij_matrix[sigmaij_matrix==1] = s1
        sigmaij_matrix[sigmaij_matrix==2] = s2
        sigmaij_matrix[sigmaij_matrix==3] = s3

        #np.random.seed(42)  # Ensures same random values each time
        gaussian_matrix = np.random.normal(0, 1, size=interaction_matrix.shape)
        gaussian_matrix = np.tril(gaussian_matrix) + np.tril(gaussian_matrix, -1).T
    
        H = Hij_matrix + sigmaij_matrix * gaussian_matrix
        np.fill_diagonal(H, diag_eng)

        return H
    
    def ipr(self, eigenvecs):
        """ Calculate the Inverse Participation Ratio (IPR) of the charge carrier.
        The IPR is a measure of the localization of the wavefunction.
        Args:
        eigenvecs (np.array): The eigenvectors of the Hamiltonian matrix
        ----------------------------------------------
        Return:
        ipr (float): The Inverse Participation Ratio
        """
        ipr = 1.0 / np.sum(np.abs(eigenvecs[0])**4, axis=0)

        return ipr

    def localization(self, sites):
        """
        Calculate the localization length of the charge carrier.
        Args:
        dist_vecs (np.array): The distance vectors array from interactions()
        interaction_matrix (np.array): The interaction matrix from interactions()
        inverse_htau (float): Inverse of the scattering time (hbar/tau) units in eV
        h_ij (np.array): The Hamiltonian matrix from hamiltonian()
        -----------------------------------------------------------------
        Return:
        lx2 (float): The localization length in x direction
        ly2 (float): The localization length in y direction
        eigenvecs (np.array): The eigenvectors of the Hamiltonian matrix for IPR
        """
        factor = -1
        if not self.is_hole: # If hole transport, it will transport at the top edge of the valence band, Boltzmann factor will be positive
            factor = 1

        beta = 1 / (kb * self.temp) # Boltzmann factor 
        h_ij = self.hamiltonian(sites) # Create Hamiltonian matrix
        energies, eigenvecs = np.linalg.eigh(h_ij) # Solve eigenvalues & eigenvectors
        sites = sites @ self.lattice.T # Back to Cartesian 
        operx = np.diag(sites[:,0])
        opery = np.diag(sites[:,1])
        weights = np.exp(-factor * energies * beta)
        partition = np.sum(weights)
    
        mxX = (eigenvecs.conj().T @ operx @ eigenvecs) # <n|x|m>, where x is the position operator
        mxY = (eigenvecs.conj().T @ opery @ eigenvecs)

        eng_diff = energies[:, None] - energies[None, :]
        mxX *= eng_diff # (En-Em) * <n|x|m>
        mxY *= eng_diff

        lx2 = sum(sum(weights * (np.abs(mxX)**2) * 2 / (self.inverse_htau**2 + eng_diff**2)))
        ly2 = sum(sum(weights * (np.abs(mxY)**2) * 2 / (self.inverse_htau**2 + eng_diff**2)))

        lx2 /= partition
        ly2 /= partition

        # Deal with IPR, where equation is 
        
        #operatorx = np.matmul(eigenvecs.T, np.matmul(dist_vecs[:,:,self.plane[0]] * h_ij, eigenvecs))
        #operatorx -= np.matmul(eigenvecs.T, np.matmul(dist_vecs[:,:,self.plane[0]] * h_ij, eigenvecs)).T

        #operatory = np.matmul(eigenvecs.T, np.matmul( dist_vecs[:,:,self.plane[1]]* h_ij, eigenvecs))
        #operatory -= np.matmul(eigenvecs.T, np.matmul(dist_vecs[:,:,self.plane[1]] * h_ij, eigenvecs)).T

        return lx2, ly2, eigenvecs

    def avg_localization(self, sites):
        """ 
        Perform average of the localization length calculation.
        Args:
        positions (np.array): The positions of the atoms in the simulation cell
        lattice (np.array): Unit cell lattice vectors
        distances (list): Specific interaction distances to consider
        j_ij (list): Inter-molecular transfer integral (J_a, J_b, J_c)
        sigma (list): dynamic disorder (sigma_a, sigma_b, sigma_c)
        inverse_htau (float): Inverse of the scattering time (hbar/tau) units in eV (Defaults to 5e-3)
        temp (float): Temperature in Kelvin (Defaults to 300)
        -----------------------------------------------------------------
        Return:
        avg_lx2 (float): The average square localization length in x direction
        avg_ly2 (float): The average square localization length in y direction
        """
        avglx2_list = []
        avgly2_list = []
        ipr_list = []
        for n in range(self.realizations):
            lx2, ly2, eigenvecs = self.localization(sites) # Calculation lx^2 and ly^2 
            ipr = self.ipr(eigenvecs) # Calculate IPR
            ipr_list.append(ipr)
            avglx2_list.append(lx2)
            avgly2_list.append(ly2)

        avglx2 = sum(avglx2_list) / self.realizations
        avgly2 = sum(avgly2_list) / self.realizations
        avgIPR = sum(ipr_list) / self.realizations

        return avglx2, avgly2, avgIPR

    def tlt_mobility(self):
        """
        TLT mobility calculation.
        Args:
        avg_lx2 (float): The average localization length in x direction in Angstrom
        avg_ly2 (float): The average localization length in y direction
        inverse_htau (float): Inverse of the scattering time (hbar/tau) units in eV
        temp (float): Temperature in Kelvin
        ---------------------------------------------------------------
        Return:
        mobilityx
        mobilityy
        mobility_average
        """
        sites = self.generate_lattice()
        avglx2, avgly2, avgIPR = self.avg_localization(sites)
        tau = hbar / self.inverse_htau # unit: second
        mobilityx = 1e-16 * e * avglx2 / (2 * tau * kb * self.temp) # Unit is cm^2/Vs
        mobilityy = 1e-16 * e * avgly2 / (2 * tau * kb * self.temp)
        mobility_average = (mobilityx + mobilityy) / 2

        return avglx2, avgly2, avgIPR, mobilityx, mobilityy, mobility_average
    
    def check_neighbors(self, dist):
        """ Check the type for each neighbors, it is for assign trnasfer integral J value
        #       #       #
    
            *      2*     *        
                                   
        #      1#      3#         
                                   
            *       *    *     
        Return values:
        1 -> 2: values = 1
        1 -> 3: values = 2
        2 -> 3: values = 3
        """
        #if np.allclose(np.abs(dist),nearest_vecs[0], atol=1e-5):    
        if np.abs(dist[0]) > 1 or np.abs(dist[1]) > 1:
            dist -= np.round(dist) # sometimes PBC can make dist really large

        if np.allclose(np.abs(dist), self.nearest_vecs[0], atol=1e-5):
            if dist[0] * dist[1] > 0:
                return 1
            else:
                return 3

        if np.allclose(np.abs(dist), self.nearest_vecs[1], atol=1e-5):
            return 2

        else: 
            return 0
        
    def marcus(self, J, deltaE=0):
        """ Calulate hopping rate using Marcus theory 
        Args:
        J (float): Transfer integral list for nearest neighbors in eV
        Lambda (float): Reorgnization energy in eV
        T (float): temperature in K
        deltaE: Energy differnece between 2 sites (if same type, 0)
        ------------------
        Return:
        k_ij: Marcus theory rate 
        """
        k_ij = (J**2) * (2*np.pi/hbar)* np.sqrt(1/(4*np.pi*self.Lambda*kb*self.temp)) * np.exp(-(deltaE+self.Lambda)**2 / (4*kb*self.temp))
    
        return round(k_ij, 5)
    
    def runKMC(self):
        """ Run kinetic Monte Carlo simulate charge transport in OSCs
        Args:
        site: np.array with center of mass of molecules in the unitcell for 2D plane
        Lambda: reorgnization energy in eV
        J_list: Transfer integral list in eV
        n: 2D supercell size
        r: cutoff radius for nearest neighbors (fractional length)
        T (float): temperature in K
        MCS: Maximum number of Monte Carlo
        ---------------------------------------------------
        Return
        traj: trajectory
        time: in ps
        """  
        sites = self.generate_lattice()
        tree = cKDTree(sites, boxsize=self.n)
        neighbor_indices = tree.query_ball_tree(tree, self.r)
        for i, neighbors in enumerate(neighbor_indices):
            neighbor_indices[i] = [j for j in neighbors if j != i] # remove self index
    
        t = 0 # time starts at 0
        time = [t]
        total_sites = len(sites)

        center = np.array([self.n / 2.0] * sites.shape[1]) # Center of the supercell
        _, idx_center = tree.query(center)
        init_pos = sites[idx_center] # Start at the center of the supercell 
        traj = [tuple(init_pos)]
    
        idx = idx_center
        pos = init_pos.copy()
    
        mcs = 0 # Monte Carlo simulation step number
        while mcs < self.realizations:
            k_list = []
            eng_idx = np.random.normal(self.Epsilon, self.sigma_ii) # Site energy at current site
            for i in neighbor_indices[idx]:
                dist = (sites[i] - sites[idx])
                eng_i = np.random.normal(self.Epsilon, self.sigma_ii) # Site energy at future site i
                values = self.check_neighbors(dist)
                if values == 1:
                    J = self.j_ij[0]
                    J += np.random.normal(0, self.sigma_ij[0]) # add dynamic disorder
                elif values == 2:
                    J = self.j_ij[1]
                    J += np.random.normal(0, self.sigma_ij[1])
                elif values == 3:
                    J = self.j_ij[2]
                    J += np.random.normal(0, self.sigma_ij[2])
                else: 
                    J = 0
                deltaE = eng_i - eng_idx # Energy difference between two sites
                k_ij = self.marcus(J, deltaE) # calculate rate constant
                k_list.append(k_ij)
                
            sum_k = np.sum(k_list)
            probs = np.array(k_list) / sum_k # probability array
            cumsum_probs = np.cumsum(probs) # cumulative probability 
        
            r1 = np.random.uniform(0, 1)
        
            for j, prob in enumerate(cumsum_probs):
                if r1 < prob:
                    nidx = neighbor_indices[idx][j]
                    break
                
            pos = sites[nidx]
            traj.append(tuple(pos))
            r2 = np.random.uniform(0, 1)
            dt = -np.log(r2) / sum_k
            idx = nidx
            t += dt
            mcs += 1
            time.append(t*1e12) # ps

        return traj, time
    
    def msd(self, traj):
        """ Get mean square displacement <r(t)^2> = <(r(t) - r(0))^2>
        Args:
        lattice (array): lattice vectors in 2D
        traj: trajectory from kmc
        --------------------------------
        Return:
        msd_vals: Mean square displacement array 
        """
        N = len(traj)
        msd_vals = np.zeros(N)

        traj_cart = np.array(traj) @ np.array(self.lattice).T # convert to cartesian

        for tau in range(N):
            disp = traj_cart[tau] - traj_cart[0] # displacement from initial position
            squared_disp = np.sum(disp**2)
            msd_vals[tau] = squared_disp

        return msd_vals
    
    def einstein_mobility(self, D):
        """ Calculate mobility using Einstein relation
        mu = eD / kT
        Args:
        D (float): Diffusion constant in cm^2/s, obtained by linregression of MSD
        ---------------------------------------------------------------
        Return: 
        mobility (float): Mobility in cm^2/Vs
        """ 
        mobility = 1e-16 * e * D / (kb*self.temp) # cm^2/Vs
    
        return mobility