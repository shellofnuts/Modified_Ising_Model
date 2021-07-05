# -*- coding: utf-8 -*-
"""
Created on Mon May 17 18:27:28 2021

@author: Robert Clampett
"""

import matplotlib.pyplot as plt

axisAni = [1.2, 1.13, 1.05, 0.98, 0.87, 0.8, 0.7, 0]
axisAni_Err = [0.01, 0.01, 0.02, 0.02, 0.03, 0.03, 0.03, 0]

axisExch = [0.9, 0.9, 0.9, 0.87, 0.82, 0.8, 0.72, 0]
axisExch_Err = [0.01, 0.01, 0.02, 0.03, 0.02, 0.03, 0.03, 0]

gamma = [1.0, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05, 0]

fig, ax = plt.subplots()

ax.grid()

ax.errorbar(gamma, axisAni, yerr=axisAni_Err,
				   capsize=5,
				   marker='o',
				   markerfacecolor = 'g',
				   label='Axis Anisotropy')

ax.errorbar(gamma, axisExch, yerr=axisExch_Err,
				   capsize=5,
				   marker='o',
				   markerfacecolor = 'g', label='Exchange Anisotropy')

ax.legend()
ax.set_ylabel(r'Critical Temperature $T_C$')
ax.set_xlabel(r'Gamma Factor $\gamma$')