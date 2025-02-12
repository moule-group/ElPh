# Electron Phonon Coupling for Molecular Crystals

# Installation:

We recommend using conda to create virtual environment.

```
    conda create -n elph
```
Then, install the necessary package.
```
    conda activate elph

    conda install phonopy numpy scipy matplotlib

    pip install --upgrade ase
```
For installing Catnip (ChArge TraNsfer Integral Package), please refer to https://joshuasbrown.github.io/docs/CATNIP/catnip_downloads.html

# Usage:

First step: Use visualize.py to identify the numbering of molecules, there are 3 molecules need to identify. The order of the numbering is shown in the figure.

Second step: Prepare input files in the folder: cif file of materials; FORCE_SETS from Phonopy simulation (8x8x8 grids); phonopy_disp.yaml from Phonopy simulation; scripts.

Third step: Run elph.py (provide the numbering of three monomers) to return electron phonon coupling parameter for further research.

# Theory:
This will divide into 2 parts. First part is transfer integral J and the second part is electron phonon coupling parameter g.

## Transfer integral J
The method we use is called dimer projection method (DIPRO) (Note: Some people call it Fragment orbital Method (FO)), it is first proposed by D. Andrienko group in 2010 research paper. 
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