# Electron Phonon Coupling for Molecular Crystals

# Installation:

Requirement: Gaussian 09 or 16

We recommend using conda to create virtual environment.

```
    git clone https://github.com/moule-group/ElPh.git

```

```
    conda create -n elph
    cd Elph
    pip install .
```

Environment variables in .bashrc

```
export PATH="your_path/catnip/build/":$PATH
export GAUSS_EXEDIR="your_path/g16"
export GAUSS_SCRDIR="$your_path/GaussianScratch"
export PATH="$GAUSS_EXEDIR:$PATH"
export PATH="$GAUSS_SCRDIR:$PATH"
```

For installing Catnip (ChArge TraNsfer Integral Package), please refer to https://joshuasbrown.github.io/docs/CATNIP/catnip_downloads.html

Make scripts in src/elph executable.

```
chmod +x * 
```

We recommend users set the environment variable 

```
export ELPH=/pscratch/sd/r/.../ElPh/src/elph
export PATH=$ELPH:$PATH
```

Example code
```
$ELPH/create_jo_input.py 
```

# Usage:

## Transfer Integral

**First step:**  Prepare input files in the folder: **CONTCAR** or **POSCAR** (VASP structure format) ; **FORCE_SETS** and **phonopy_disp.yaml** from Phonopy simulation.

**Second step:** Run create_j0_input.py, which will generate the input files to calculate transfer integral J at reference state. Then loop through every subfolder to execute run.bash. This will conduct gaussian simulation.

**Third step:** Run j0.py. This will call Catnip to calculate transfer integral based on Gaussian output files. the result will be written into j0.json.

## Electron Phonon Coupling

**First step:** Run create_disp_input.py. This will generate displaced structures in "displacements" folder. It will also return "run_all.sh" script for users to conduct Gaussian simulations.

**Second step:** User should finish all Gaussian simulations for 1 dimer and 2 monomers. (Ex: dimer A; monomer 1 and 2. dimer B; monomer 1 and 3). Then run disp_j.py, which will generate a run script for user to run catnip calculations.

## Variance and Projection

Please check the jupyter notebook in example folder. There are 3 different materials as tutorials. 

# Theory:
This will divide into 3 parts. First part is transfer integral J, the second part is electron phonon coupling parameter g and the last part is transient localization theory.

## Transfer integral J, Electron-Phonon Coupling (EPC), Variance, and Singular Value Decomposition (SVD)

We refer the user to the research article we are preparing.

## Transient Localization Theory (TLT)

The mobility equation is shown below:
$\mu =  \frac{e}{kT} \frac{L^2_{x(y)}}{2\tau}$

where $\tau$ is the relaxation time, $L^2_{x(y)}$ is squared localization length, e is the charge, T is the temperature in K, k is the Boltzmann constant. The mobility unit is in $\frac{cm^2}{Vs}$

## Singular Value Decomposition (SVD) Analysis 

Please check the example folder, there is a tutorial of BTBT dimer.
