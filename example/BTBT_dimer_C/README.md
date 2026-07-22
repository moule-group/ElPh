# Example of BTBT dimer analysis

This folder contains several files that will be used to conduct electron-phonon coupling (EPC) and singular value decomposition (SVD) analysis.

## Necessary input files

1. Force constant file from phonopy: FORCE_SETS
2. Phonopy displacement file: phonopy_disp.yaml
3. Crystal geometry: CONTCAR (POSCAR)
4. Dimer geometry: dimer_C.xyz 

## Notebook to run SVD analysis

svd_phonons.ipynb is the notebook that contain all functions for EPC and SVD analysis

Input parameters:
1. temp = 300 # TEMPERATURE: Kelvin
2. mesh = [8,8,8] # phonopy mesh grid
3. dimer_label = 'c' # label for dimer: {A:(1,2)} {B:(1,3)} {C:(1,4)} {D:(2,3)} {E:(2,4)} {F:(3,4)}
4. monomer1 = 1 # label for monomer 1
5. monomer2 = 4 # label for monomer 2
6. molecule = 'BTBT' # Will be used on naming for the plots
7. crystal = 'CONTCAR' # The unitcell file path for phonopy simulation
8. xyz = 'dimer_C.xyz' # dimer xyz file path
9. n_atoms_xyz = 48 # The number of atoms in dimer file (first line in .xyz file)
10. dof = n_atoms_xyz * 3 # Degrees of freedom in dimer
11. xtb = 'n' # Transfer integral simulation is done by xtb or not (y/n)

The code will then request additional input to determine whether the two molecules in the dimer are related by lattice translation.

For example, if the lattice translation vector is

`Lattice translational vector: [0. 1. 0.]`

the user should enter the vector components as space-separated numbers:

`0 1 0`

It will also asked the dimer is formed by tranlsation molecules or not (y/n):

`y`

If the dimer consists of two translationally related molecules, the user must also provide the atomic indices corresponding to one of the molecules. For example:

`0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46`

If the two molecules are not related by lattice translation, the user should press the **Enter** key without entering any indices. The code will then select all atomic indices.

## Time-dependent Trajectories 

After running svd_phonons.ipynb, user can then run write_svd_traj.py to obtain time-dependent trajectory that can be convert into neutron spectra and density of states spectra. 

`python write_svd_traj.py --xyz dimer_C_map.xyz`

options: 
1. dt: Trajectory timestep in fs (default: 1)
2. nframes: Number of trajectory frames (default: 20001)
3. amp: Visualization scaling factor (default: 80)
4. output-prefix: Output xyz prefix (default: proj_traj)
5. start_index: The index for the first PM (default: 0)
