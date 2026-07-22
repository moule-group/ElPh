#!/usr/bin/env python3
import argparse
import ase.io
import numpy as np

def trajectory_proj_mode(amp, dt, freqs, nframes, output, p, population_iq, X_allq, xyz):
    """ This is the function to generate atomic trajectory for linear combination of orginal phonon modes for each projected modes based on the contribution coefficient
    Idea: The orginal phonon modes have distinct frequencies and the projected mode will be like Fourier Series.
    Args:
        amp (float): visualization scaling vector
        dt: timestep (fs)
        freqs: Frequencies from phonopy
        nframes (int): number of frames for the trajectory
        output: The output .xyz file name
        p (int): The projected mode index, user should select the index
        population_iq: original phonon mode population 
        X_allq: X matrix
        xyz (str): dimer xyz file name
    Return: 
        The .xyz file which saves the trajectory of the molecules
    """
    atoms_xyz = ase.io.read(xyz)
    n_atoms_xyz = len(atoms_xyz)
    atoms_element = atoms_xyz.get_chemical_symbols()
    atoms_coords = atoms_xyz.get_positions()

    # Generate frames and write
    with open(f'{output}.xyz', "w") as f:
        for t in range(1, nframes):
            time = t * dt * 1e-15 # in seconds
            omega_t = 2 * np.pi * freqs * 1e12 * time # freq * time
            factor = np.exp(1j * omega_t) 
            factor *= population_iq[p,:] # Multiply the coefficient (population of each original mode)

            disp_matrix = factor * X_allq 
            disp = np.sum(disp_matrix, axis=1) # Sum over all modes 
            disp_xyz = disp.reshape(n_atoms_xyz, 3)
            atoms_new_coords = atoms_coords + amp * disp_xyz.real

            f.write(f"{n_atoms_xyz}\n") # XYZ frame
            f.write(f"phonon frame={t} time={time} PM-index={p}\n")
            for ele, (x, y, z) in zip(atoms_element, atoms_new_coords):
                f.write(f"{ele} {x:.8f} {y:.8f} {z:.8f}\n")

    return output

def build_parser():
    parser = argparse.ArgumentParser(
        description="Write projected-mode trajectory xyz files from arrays saved by the notebook."
    )
    parser.add_argument(
        "--data",
        default="trajectory_inputs.npz",
        help="NPZ file containing X_allq, population_iq, and freqs_flatten",
    )
    parser.add_argument("--xyz", required=True, help="Mapped dimer xyz file, for example dimer_A_map.xyz")
    parser.add_argument("--dt", type=float, default=1, help="Trajectory timestep in fs")
    parser.add_argument("--nframes", type=int, default=20001, help="Number of trajectory frames")
    parser.add_argument("--amp", type=float, default=80, help="Visualization scaling factor")
    parser.add_argument("--output-prefix", default="proj_traj", help="Output xyz prefix")
    parser.add_argument("--start_index", type=int, default=0, help="The index for the first PM")  
    
    return parser

def main():
    args = build_parser().parse_args()
    data = np.load("trajectory_inputs.npz")
    X_allq = data["X"]
    population_iq = data["pop"]
    freqs_flatten = data["freq"]
    index = data["index"]
    index_user = index[args.start_index:]

    for i, i_value in enumerate(index_user):
        output = f"{args.output_prefix}_{i+args.start_index+1}"
        trajectory_proj_mode(
            args.amp,
            args.dt,
            freqs_flatten,
            args.nframes,
            output,
            i_value,
            population_iq,
            X_allq,
            args.xyz,
        )

    print("--- Finish generating trajectory file ! ---")

if __name__ == "__main__":
    main()