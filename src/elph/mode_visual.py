##### This code is for visualizing the map of the electron-phonon coupling w.r.t phonon frequency #####

import numpy as np
import matplotlib.pyplot as plt
import yaml

# Define constants
#kb = 1.38 * 1e-23 # Boltzmann constant (J/K)
#temp = 298 # Temperature, we set as room temperature (K)

# Load Phonon frequency from Phonopy mesh.yaml
with open("mesh.yaml", "r") as f1:
    data1 = yaml.safe_load(f1)

# Frequency of phonon modes

nqpts = data1["nqpoint"] # Number of q-points

freq_list = [] #  Need to check unit here
for q in range(nqpts):
    freqs = [band["frequency"] for band in data1["phonon"][q]["band"]]
    freq_list.append(freqs)

# Load electron-phonon coupling from elph.py
with open("ep_coupling.yaml", "r") as f2:
    data2 = yaml.safe_load(f2)



# the vibrational amplitude of the phonon mode

# We need to deal with reduced mass

# reduced_mass = (m1 * m2) / (m1 + m2) # Reduced mass of the two atoms

# amp_list = []

# for freq in freq_list:
#      amp = np.sqrt( (kb*temp)/(reduced_mass*(freq**2)) ) # Vibrational amplitude A_{k}
#      amp_list.append(amp)

# Fluctuation of the charge transfer integral sigma = g_{ij} * A_{k}

# sigma = 



# Plot
fig, ax = plt.subplots()

#for freq in freq_list:
 #   ax.bar(freq, ep * 30, width=2, color='black')

ax.set_xlabel(r'$Energy$ (THz)')
ax.set_ylabel(r'$g_{ij}$ (meV/$\AA$)')
ax.set_ylim(0, 35)
plt.show()
