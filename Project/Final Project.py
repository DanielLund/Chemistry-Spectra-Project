import numpy as np
from scipy.constants import h, k, Avogadro
from math import pi
import matplotlib.pyplot as plt


def reduced_mass(m1, m2):
    return (m1*m2)/(m1+m2)


def moment_of_inertia(mu, r):
    return mu*(r**2)


rotational_prefactor = h**2/(8*(pi**2))


def rotational_constant(moment_of_inertia):
    return rotational_prefactor/moment_of_inertia


def rotational_energy(j, rotational_constant):
    return rotational_constant*j*(j+1)


def transition_energy(j1, j2, rotational_constant):
    energy_1 = rotational_energy(j1, rotational_constant)
    energy_2 = rotational_energy(j2, rotational_constant)
    return energy_2 - energy_1


def transition_probability(j_initial, j_final):
    if j_final - j_initial == 1:
        return 1.0
    elif j_final - j_initial == -1:
        return 1.0
    else:
        return 0


def degeneracy(j):
    return 2*j + 1


def boltzmann(j, rotational_constant, temperature):
    return (2*j+1)*np.exp(-(rotational_constant*j*(j+1))/k*temperature)


def intensity(j, rotational_constant, temperature):
    return ((h*rotational_constant)/(k*temperature))*((2*j)+1)\
           * np.exp(-rotational_energy(j, rotational_constant)/(k*temperature))


m_H = 1.66e-27
m_Cl = 5.81e-26
r_HCl = 1.29e-10


frequencies_in_GHz = []
intensities = []
max_j = 11
temperature = 298.0

mu = reduced_mass(m_H, m_Cl)
I = moment_of_inertia(mu, r_HCl)
B = rotational_constant(I)

for j_initial in range(0, max_j):
    for j_final in range(j_initial+1, max_j):
        if transition_probability(j_initial, j_final) == 1:
            frequencies_in_GHz.append(transition_energy(j_initial, j_final, B)/(h*1e9))
            intensities.append(intensity(j_final, B, temperature))

plt.stem(frequencies_in_GHz, intensities, basefmt='-')
plt.show()


frequencies_in_GHz = []
intensities = []
max_j = 20
temperature = 100 # Kelvin

mu = reduced_mass(m_H, m_Cl)
I = moment_of_inertia(mu, r_HCl)
B = rotational_constant(I)

for j_initial in range(0, max_j):
    for j_final in range(j_initial+1, max_j):
        if transition_probability(j_initial, j_final) == 1:
            frequencies_in_GHz.append(transition_energy(j_initial, j_final, B)/(h*1e9))
            intensities.append(intensity(j_final, B, temperature))

plt.stem(frequencies_in_GHz, intensities, basefmt='-')
plt.plot(frequencies_in_GHz, intensities, '--', label=temperature)
plt.show()


frequencies_in_GHz = []
intensities = []
max_j = 20
temperature = 298

mu = reduced_mass(m_H, m_Cl)
I = moment_of_inertia(mu, r_HCl)
B = rotational_constant(I)

for j_initial in range(0, max_j):
    for j_final in range(j_initial+1, max_j):
        if transition_probability(j_initial, j_final) == 1:
            frequencies_in_GHz.append(transition_energy(j_initial, j_final, B)/(h*1e9))
            intensities.append(intensity(j_final, B, temperature))

plt.stem(frequencies_in_GHz, intensities, basefmt='-')
plt.plot(frequencies_in_GHz, intensities, '--', label=temperature)
plt.show()


frequencies_in_GHz = []
intensities = []
max_j = 20
temperature = 1000

mu = reduced_mass(m_H, m_Cl)
I = moment_of_inertia(mu, r_HCl)
B = rotational_constant(I)

for j_initial in range(0, max_j):
    for j_final in range(j_initial+1, max_j):
        if transition_probability(j_initial, j_final) == 1:
            frequencies_in_GHz.append(transition_energy(j_initial, j_final, B)/(h*1e9))
            intensities.append(intensity(j_final, B, temperature))

plt.stem(frequencies_in_GHz, intensities, basefmt='-')
plt.plot(frequencies_in_GHz, intensities, '--', label=temperature)
plt.show()


"""Plotting raman spectra"""


def transition_probability_raman(j_initial, j_final):
    if j_final - j_initial == +2.0:
        return 2.0
    elif j_final - j_initial == -2.0:
        return -2.0
    else:
        return 0


r_CO = 1.13e-10
m_C = 12 / Avogadro * 1e-3
m_O = 16 / Avogadro * 1e-3

frequencies_in_GHz = []
intensities = []
max_j = 21
temperature = 298  # Kelvin

mu = reduced_mass(m_C, m_O)
I = moment_of_inertia(mu, r_CO)
B = rotational_constant(I)

for j_initial in range(0, max_j):
    for j_final in range(0, max_j):
        if transition_probability(j_initial, j_final) == 1:
            frequencies_in_GHz.append(transition_energy(j_initial, j_final, B) / (h * 1e9))
            intensities.append(intensity(j_final, B, temperature))

plt.stem(frequencies_in_GHz, intensities, basefmt='-')

plt.show()
