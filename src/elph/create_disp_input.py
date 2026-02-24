#!/usr/bin/env python3
from curses import echo

import ase.io
import glob
import os
from ase.calculators.gaussian import Gaussian
from pathlib import Path
import time

def getGeometry(path):
    """ Using glob function in python to find the structure
    file in the current path
    The type of the structure files: "POSCAR", "CONTCAR", "*.xyz"
    Args:
        path: The current directory (use os.getcwd())
    ----------------------------------------------
    Return:
        file: The structure file in the path
    """
    file = glob.glob(path + "/POSCAR") + glob.glob(path + "/CONTCAR") + glob.glob(path + "/*.xyz") 
    if len(file) == 0:
        raise FileNotFoundError
    return file[0]

def get_displacement(atoms):
    """ Get numbering of displaced atom, displacement direction (x,y,z) and sign (+,-) 
    Args:
        atoms (ASE Atom object): molecule or dimer 
    -----------------------------------------------
    Return:
        na: numbering of displaced atom
        vec: displacement direction (x,y,z) in (0,1,2)
        sign: (+,-) in (+1, -1)
    """
    num_atoms = len(atoms) # number of atoms
    for na in range(num_atoms):
        for vec in range(3):
            for sign in [-1, 1]:
                yield (na,vec,sign) 

def write_run_script():
    """ Write the run script in each subdirectory for G16 simulation
    """
    string = """#!/bin/bash
folder=$(basename "$PWD")
g16 < mo.com > mo.log
mv fort.7 "${folder}.pun"
    """
    run_file = Path("run.sh")
    run_file.write_text(string)
    os.chmod(run_file, 0o755)

def write_all_run_script(label):
    """ Write the run script in the directory to 
        run all the G16 simulations for displaced structures
    Args:
        label (str): Label for the folder to create displaced structures
    """
    base_dir = os.getcwd()
    cwd = os.path.join(base_dir, label, 'displacements')

    script_content = f"""#!/bin/bash
set -euo pipefail

# Base directory
CWD="{cwd}"
DISP_LIST="${{CWD}}/disp_list.txt"

cd "${{CWD}}"

echo "Starting serial displacement runs on $(hostname) at $(date)"

while read -r DIR; do
    echo "--------------------------------------------"
    echo "Running ${{DIR}} at $(date)"

    cd "${{CWD}}/${{DIR}}"
    bash run.sh

    # Record completion
    echo "${{DIR}}/" >> "${{CWD}}/DONE"

    cd "${{CWD}}"
done < "${{DISP_LIST}}"
echo "All displacements finished at $(date)"
"""
    script_path = Path("run_all.sh")
    script_path.write_text(script_content)
    script_path.chmod(0o755)

def mol_orbital(basis_set, functional, disper_corr):
    """ Run Gaussian to compute the molecular orbitals coefficient and energy for the system (Single point calculation)
    Args:
        basis_set (str): Basis set for DFT calculation 
        functional (str): Functional for DFT calculation 
        disper_corr (int): 1 is D3, 2 is D3BJ
    -----------------------------------------------
    Return:
        .pun file which contains molecular orbital for the cacluation later for J 
        .log file output file from Gaussian
    """
    path = os.getcwd()
    geometry = getGeometry(path)
    atoms = ase.io.read(geometry)

    if disper_corr == 1:
        dispersion = 'EmpiricalDispersion=GD3'
    else:
        dispersion = 'EmpiricalDispersion=GD3BJ'
    
    atoms.calc = Gaussian(mem='18GB',
                 nprocshared=12,
                 label='mo',
                 save=None,
                 method=functional,
                 basis=basis_set,  
                 scf='tight',
                 pop='full',
                 extra=f'nosymm punch=mo iop(3/33=1) {dispersion}') # iop(3/33=1) output one-electron integrals to log file.
    atoms.calc.write_input(atoms)

def create_displacement(label, basis_set, functional, dispersion_correction, delta=0.01):
    """ Create atomic displacement and return each displaced structure; generate G16 input files at each subdirectory
    Args: 
        label (str): Label for the folder to create displaced structures
        basis_sets (str): DFT basis sets
        functional (str): DFT functional
        dispersion_correction (int): Dispersion correction (1-> D3, 2-> D3-BJ)
        delta (float): Magnitude of displacement in Angstrom (Defaults to 0.01A)
    """
    os.chdir(label)
    path = os.getcwd()
    geometry = getGeometry(path)
    atoms = ase.io.read(geometry)
    os.makedirs('displacements', exist_ok=True)
    os.chdir('displacements')
             
    with open("disp_list.txt", "w") as f:
        for na, vec, sign in get_displacement(atoms):
            disp_atoms = atoms.copy()
            pos = disp_atoms.get_positions()
            pos[na, vec] += delta * sign
            disp_atoms.set_positions(pos)
            disp_name = f"disp_{na}_{vec}_{sign}"
            f.write(disp_name + "\n")
            os.makedirs(disp_name, exist_ok=True)
            os.chdir(disp_name)
            ase.io.write(disp_name+'.xyz', disp_atoms)
            mol_orbital(basis_set, functional, dispersion_correction)
            write_run_script()
            os.chdir('..')

def run_disp_mo(label, basis_sets, functional, dispersion_correction):
    """ Main function to generate Gaussian input file and run.sh for displaced molecules 
    and dimers, which can calculate molecular orbitals
    Args:
        label (str): Label for the folder to create displaced structures
        basis_sets (str): DFT basis sets
        functional (str): DFT functional
        dispersion_correction (int): Dispersion correction (1-> D3, 2-> D3-BJ)
    ------------------------------------------------------------
    Return:
    j_list (list): The transfer integral list for displaced molecules and dimers!
    """
    path = os.getcwd()
    write_all_run_script(label)
    if not os.path.exists(os.path.join(path, label, 'displacements')):
        print(' Creating displaced molecules and dimers ... ')
        create_displacement(label, basis_sets, functional, dispersion_correction)
        
    print(f" Displacement folder {label} is generated and run_all.sh script has been written! ")

def ask_with_default(prompt, default, type=str):
    """ Ask user for input with a default value. If user presses Enter, return default.
    Args:
        prompt (str): The prompt to display to the user
        default: The default value to return if user presses Enter
        type: The type to cast the user input to (default is str)
    """
    user_input = input(f"{prompt} [{default}]: ").strip()
    if user_input == "":
        return default
    return type(user_input)

if __name__ == "__main__":

    print(" ----------------------------------------------------------------- ")
    print(" Starting on Electron Phonon Coupling (ElPh) for Organic Semiconductors !!! ")
    print(" ElPh code made by Moule Group at UC Davis. Author: Ray Chiang ")
    print(" #########   $           $ ########   $        $  ")
    print(" $           $           $        $   $        $  ")
    print(" $########   $           $ #######    $########$  ")
    print(" $           $           $            $        $  ")
    print(" #########   ##########  $            $        $  ")
    print(" ----------------------------------------------------------------- ")

    print("This script allows to generate displaced structures and generate gaussian input files\n"
          "for molecular orbital calculation. The run script will be generated.\n")
    
    label = ask_with_default(
    "Label for the folder to create displaced structures",
    "1",str)

    basis_set = ask_with_default(
    "Basis set. Example: Def2SVP, TZVP, Def2TZVP .... Default is TZVP",
    "TZVP")

    functional = ask_with_default(
    "Functional. Example: B3LYP, PBE0, M062X .... Default is M062X",
    "M062X") 
    
    dispersion_correction = ask_with_default(
    "Dispersion correction (1-> D3, 2-> D3BJ)",
    1,int)

    print("")
    print("-------------------------------\n"
      "Input parameters:\n"
      f"Label for the folder to create displaced structures: {label}\n"
      f"basis sets: {basis_set}\n"
      f"functional: {functional}\n"
      f"dispersion correction (1: D3, 2: D3-BJ): {dispersion_correction}\n"
      "-------------------------------"
     )
    print("")

    run_disp_mo(label, basis_set, functional, dispersion_correction)

    print("")
    print(" ----------------------------------------------------------------- ")
    print(
    """                 
    #########  ##    #  #########
    $          # #   #  $       $
    ########   #  #  #  $      $
    $          #   # #  $     $ 
    #########  #    ##  ######
    """)
    print(" ----------------------------------------------------------------- ")
    print(time.ctime())