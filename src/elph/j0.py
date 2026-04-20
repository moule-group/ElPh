#!/usr/bin/env python3
import subprocess
import time

def run_catnip(path1, path2, path3, path4, path5, path6):
    """ Run Catnip to calculate the transfer integral
    Args:
        path1 (str): Path to the first Gaussian .pun file
        path2 (str): Path to the second Gaussian .pun file
        path3 (str): Path to the third Gaussian .pun file (pair)
        path4 (str): Path to the first Gaussian .log file
        path5 (str): Path to the second Gaussian .log file
        path6 (str): Path to the third Gaussian .log file (pair)
    ----------------------------------------------------------------
    Returns:
        output.decode('ascii').split()[-2]: Transfer integral J_eff (Effective Transfer Integral) for the system, units is eV
        output.decode('ascii').split()[-13]: Transfer integral J for the system, units is eV
    """
    cmd = f"calc_J -p_1 {path1} -p_2 {path2} -p_P {path3} -l_1 {path4} -l_2 {path5} -l_P {path6}"
    output = subprocess.check_output(cmd, shell=True)
    
    return output.decode('ascii').split()[-2], output.decode('ascii').split()[-13]

def run_j0(monomer1, monomer2, dimer):
    """ Main function for running Catnip to get transfer integral J_0
    Args:
        monomer1 (int): The numbering of monomer 1
        monomer2 (int): The numbering of monomer 2
        dimer (str): The label for the dimer
    Return:
        {dimer}.out file for transfer integral
    """    
    print(" --- Start transfer integral J calculation ---")      # A, B, C, ...

    j_eff, j0 = run_catnip(
            f'./{monomer1}/{monomer1}.pun',
            f'./{monomer2}/{monomer2}.pun',
            f'./{dimer}/{dimer}.pun',
            f'./{monomer1}/mo.log',
            f'./{monomer2}/mo.log',
            f'./{dimer}/mo.log'
        )
    with open(f"{dimer}.out", "w") as f:
        f.write(f"J_eff: {j_eff} eV\n")
        f.write(f"J_0: {j0} eV\n")

    print(f" --- The result has been saved to the {dimer}.out file --- ")
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

    dimer = ask_with_default(
    "Label for the dimer folder",
    "A",str)
    monomer1 = ask_with_default(
    "Label for the monomer1 folder",
    1,str)   
    monomer2 = ask_with_default(
    "Label for the monomer2 folder",
    2,str)
    
    run_j0(monomer1, monomer2, dimer)