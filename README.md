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

# Usage:

## Non Local Electron Phonon Coupling

First step:  Input the number of molecules (defaults to 3) needs to be extracted, the code will generate monomer and dimer structure files.

Second step: Prepare input files in the folder: CONTCAR from vasp output and/or cif file of materials; FORCE_SETS from Phonopy simulation; phonopy_disp.yaml from Phonopy simulation.

Third step: Run elph (provide the numbering of three monomers) to return electron phonon coupling parameter for further research.

Note: We consider 2D plane (high mobility plane of organic semiconductors) and only pick 3 nearest neighbors in this 2D plane. The 3 numbering monomers will be pair A (monomer 1 and 2); pair B (monomer 1 and 3); pair C (monomer 2 and 3), pair A and pair B will be transversed pairs and pair C will be parallel pairs (the shorter lattice parameter in 2D plane).

```
elph -w 1
```

## SVD Phonon Modes Projection

Prerequisite: Finish non-local epc simulation first.

```
elph -w 2 -svd 128
```

## Transient Localization Theory Charge Carrier Mobility

Prepare mobility.json file as the input, then run

```
elph -w 3
```

## Arguments

-w --workflow: Workflow selection (1: non local electron phonon coupling 2: SVD phonon mode projection 3: TLT mobility)

-q --mesh: Defining a mesh grid. (Defaults to [8,8,8])

-n --nmol: The number of molecules will be extracted (Defaults to 3)

-b --basis: Gaussian basis sets (Defaults to ['6-311G*','6-311G**']) for local and non-local simulation

-f --functional: Gaussian functional (Defaults to['b3lyp','b3lyp']) for local and non-local simulation

-s --supercell: The supercell matrix (Defaults to [2,2,2])

-homo --homo: P-type semiconductors: HOMO; N-type semiconductors: LUMO. (Defaults to True)

-o --output: Mobility calculation output name (Defaults to tlt_mobility.json)

-svd --svdqpts: Number of qpoints that SVD projection will apply (Defaults to 1)

-mx --matrix: The matrix (epc or variance) that be applied to SVD projection (epc or var, defaults to epc)

## mobility.json

In order to run mobility calculation, there are variables need to be specified. Please see the example here. [View Example File](example/mobility.json)


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
