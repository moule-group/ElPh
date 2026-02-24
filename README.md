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

## Transfer integral J
The method we use is called dimer projection method (DIPRO) (Note: Some people call it Fragment orbital Method (FO)), it is proposed by D. Andrienko group in 2010 research paper. 
Transfer integral between 2 molecules i and j:
$J_{ij} = <i|H|j>$
The reason why we cannot simply the equation above is because overlap matrix S is not 0 in molecular crystals. Therefore, we have to apply effective transfer integral

Effective Transfer Integral:
$J_{eff}= \frac{J_{ij}-\frac{S_{ij}(E_i+E_j)}{2}}{1-S_{ij}^2}$

In order to use DIPRO to calculate transfer integral, we have to run 3 quantum-chemical simulations (2 monomers and 1 dimer), there are 9 quantum-chemical simulations in total.

## Electron Phonon Coupling Parameter g

It can be further written as 
$J_{ij} = J_{ij}^0 + \sum_{I} g_{ij}^IQ_{I}$,
where $J_{ij}$ is transfer integral at equilibrium geometry. 
$g_{ij}^I$ is the electron-phonon coupling parameter, 
$Q_{I}$ is normal coordinate at vibrational mode I.

In order to efficiently calculate the electron-phonon coupling parameter $g_{ij}^I$, 
we can do the conversion as below.

$g_{ij} = \frac{\partial J_{ij}}{\partial Q_{I}}$

Using chain rule, we get

$g_{ij} = \nabla J_{ij}\frac{\partial x_{k}}{\partial Q_{I}}$

where $x_{k}$ is Cartesian coordinate of the molecule. 
And $\frac{\partial x_{k}}{\partial Q_{I}}$ represent how each Cartesian coordinate changes with the normal mode coordinate, which is calculated by Phonopy modulation. The equation of modulation (creation of crystal structures with displacements) is 

![Equation](https://latex.codecogs.com/svg.image?\frac{A}{\sqrt{N_a&space;m_j}}\,\text{Re}\left[\exp(i\phi)e_j\exp(i\mathbf{q}\cdot\mathbf{r}_{jl})\right])

, where A is the amplitude (Defaults to 0.01 Angstrom), $N_{a}$ is the number of atoms in the supercell, $m_{j}$ is the mass of j-th atom. $r_{jl}$ is the position of the j-th atom in the l-th unit cell, $e_{j}$ is the j-th atom part of eigenvector, and $\phi$ is the phase.

To evaluate $\nabla J_{ij}$, a displacements -0.01 Å and 0.01 Å for each direction (x,y,z) of the gradient have been employed.

## Transient Localization Theory (TLT)

The mobility equation is shown below:
$\mu =  \frac{e}{kT} \frac{L^2_{x(y)}}{2\tau}$

where $\tau$ is the relaxation time, $L^2_{x(y)}$ is squared localization length, e is the charge, T is the temperature in K, k is the Boltzmann constant. The mobility unit is in $\frac{cm^2}{Vs}$
