##### This code is for visualizing the map of the electron-phonon coupling w.r.t phonon frequency #####

import numpy as np
import matplotlib.pyplot as plt
import yaml

# Define constants
#kb = 1.38 * 1e-23 # Boltzmann constant (J/K)
#temp = 298 # Temperature, we set as room temperature (K)
thztocm_1 = 33.356 # Conversion factor from THz to cm-1

# Load Phonon frequency from Phonopy mesh.yaml
with open("mesh.yaml", "r") as f1:
    data1 = yaml.safe_load(f1)

# Frequency of phonon modes

nqpts = data1["nqpoint"] # Number of q-points

freq_list = [] #  Need to check unit here
for q in range(nqpts):
    freqs = [band["frequency"] for band in data1["phonon"][q]["band"]]
    freq_list.append(freqs)

# Load electron-phonon coupling 
epc = np.load('ep_coupling.npz')



# the vibrational amplitude of the phonon mode

# We need to deal with reduced mass

# reduced_mass = (m1 * m2) / (m1 + m2) # Reduced mass of the two atoms

# amp_list = []

# for freq in freq_list:
#      amp = np.sqrt( (kb*temp)/(reduced_mass*(freq**2)) ) # Vibrational amplitude A_{k}
#      amp_list.append(amp)

# Fluctuation of the charge transfer integral sigma = g_{ij} * A_{k}

# sigma = 

fig, axes = plt.subplots(3, 1, sharex=True, figsize=(6, 8))

for q in range(len(freq_list)):
    for f in range(len(freq_list[q])):
        axes[0].bar(freq_list[q][f]*thztocm_1, epc['A'][len(freq_list[q])*q + f] * 10000, width=0.5, color='red')

axes[0].axhline(y=0, color='black', linestyle='-', linewidth=0.5) 
axes[0].set_ylabel(r'$g_{ij}$ x 1$e^4$ (meV/$\AA$)')
axes[0].set_title('Electron Phonon Coupling Parameter vs Frequency for Dimer A')
axes[0].set_ylim(-0.5, 0.5)

for q in range(len(freq_list)):
    for f in range(len(freq_list[q])):
        axes[1].bar(freq_list[q][f]*thztocm_1, epc['B'][len(freq_list[q])*q + f] * 10000, width=0.5, color='blue')

axes[1].axhline(y=0, color='black', linestyle='-', linewidth=0.5) 
axes[1].set_ylabel(r'$g_{ij}$ x 1$e^4$ (meV/$\AA$)')
axes[1].set_title('Electron Phonon Coupling Parameter vs Frequency for Dimer B')
axes[1].set_ylim(-0.5, 0.5)

for q in range(len(freq_list)):
    for f in range(len(freq_list[q])):
        axes[2].bar(freq_list[q][f]*thztocm_1, epc['C'][len(freq_list[q])*q + f] * 10000, width=0.5, color='green')

axes[2].axhline(y=0, color='black', linestyle='-', linewidth=0.5) 
axes[2].set_ylabel(r'$g_{ij}$ x 1$e^4$ (meV/$\AA$)')
axes[2].set_title('Electron Phonon Coupling Parameter vs Frequency for Dimer C')
axes[2].set_ylim(-0.5, 0.5)


axes[-1].set_xlim(-10, 3500) 
axes[-1].set_xlabel(r'$Energy$ (cm-1)')
plt.show()


