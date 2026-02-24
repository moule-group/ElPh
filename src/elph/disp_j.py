#!/usr/bin/env python3
import ase.io
import os
import time

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

def generate_J_task(dimer, m1, m2):
    """ The function to generate the task to run J for displaced molecules and dimers
    Args: 
        dimer (str): The label for the dimer
        m1 (int): The numbering of monomer 1
        m2 (int): The numbering of monomer 2
    """
    molecules = ase.io.read(os.path.join(f'{m1}', f'monomer_{m1}.xyz'))
    offset = len(molecules)
    disp_list = list(get_displacement(molecules))
    # Loop 1: displace mol_1 (dimer uses same na)
    for na, vec, sign in disp_list:
        yield {
            "which": "mol1",
            "na": na,
            "vec": vec,
            "sign": sign,
            "offset": offset,
            "p1": f"./{m1}/displacements/disp_{na}_{vec}_{sign}/disp_{na}_{vec}_{sign}.pun",
            "p2": f"./{m2}/{m2}.pun",
            "pP": f"./{dimer}/displacements/disp_{na}_{vec}_{sign}/disp_{na}_{vec}_{sign}.pun",
            "l1": f"./{m1}/displacements/disp_{na}_{vec}_{sign}/mo.log",
            "l2": f"./{m2}/mo.log",
            "lP": f"./{dimer}/displacements/disp_{na}_{vec}_{sign}/mo.log",
        }
    # Loop 2: displace mol_2 (dimer uses na+offset)
    for na, vec, sign in disp_list:
        naP = na + offset
        yield {
            "which": "mol2",
            "na": na,
            "vec": vec,
            "sign": sign,
            "offset": offset,
            "p1": f"./{m1}/{m1}.pun",
            "p2": f"./{m2}/displacements/disp_{na}_{vec}_{sign}/disp_{na}_{vec}_{sign}.pun",
            "pP": f"./{dimer}/displacements/disp_{naP}_{vec}_{sign}/disp_{naP}_{vec}_{sign}.pun",
            "l1": f"./{m1}/mo.log",
            "l2": f"./{m2}/displacements/disp_{na}_{vec}_{sign}/mo.log",
            "lP": f"./{dimer}/displacements/disp_{naP}_{vec}_{sign}/mo.log",
        }

def write_sh(dimer, m1, m2):
    """ Write the .sh file to run catnip
    Args:
        dimer (str): The label for the dimer
        m1 (int): The numbering of monomer 1
        m2 (int): The numbering of monomer 2
    """
    out_sh  = f"dj{dimer}_tasks.sh"
    tasks = list(generate_J_task(dimer, str(m1), str(m2)))   
    with open(out_sh, "w") as f:
        f.write("#!/usr/bin/env bash\n")
        for i, t in enumerate(tasks, start=1):
            out = f"{i:04d}.out"
            err = f"{i:04d}.err"
            cmd = (
                f'calc_J -p_1 {(t["p1"])} -p_2 {(t["p2"])} -p_P {(t["pP"])} '
                f'-l_1 {(t["l1"])} -l_2 {(t["l2"])} -l_P {(t["lP"])} '
                f'1> "{out}" 2> "{err}"\n'
            )
            f.write(f"# task_id={i} dimer={dimer} which={t['which']} na={t['na']} vec={t['vec']} sign={t['sign']}\n")
            f.write(cmd + "\n")
            f.write(f'JAB=$(awk \'/^J_ab/  {{print $2}}\' "{out}")\n')
            f.write(f'JEFF=$(awk \'/^J_eff/ {{print $2}}\' "{out}")\n')
            f.write(f"echo $JAB $JEFF >> disp_j/{dimer}.txt\n")
            f.write(f'rm -f "{out}" "{err}"\n\n')

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

    print(" This script allows to generate the task file to run J for displaced molecules and dimers\n"
            ", which can calculate the transfer integral change upon molecular displacement. ")

    dimer = ask_with_default(
    "Label for the dimer folder to create displaced structures",
    "A",str)
    monomer1 = ask_with_default(
    "Label for the monomer1 folder to create displaced structures",
    1,str)   
    monomer2 = ask_with_default(
    "Label for the monomer2 folder to create displaced structures",
    2,str)

    write_sh(dimer,monomer1,monomer2)
    print(" --------------------------------------------------------- ")
    print(f" The task file dj{dimer}_tasks.sh has been generated! ")
    print("")
    print(" --------------------------------------------------------- ")
    print(
    """                 
    #########  ##    #  #########
    $          # #   #  $       $
    ########   #  #  #  $      $
    $          #   # #  $     $ 
    #########  #    ##  ######
    """)
    print(" --------------------------------------------------------- ")
    print(time.ctime())
