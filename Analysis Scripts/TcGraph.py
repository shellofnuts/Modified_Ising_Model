import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np

plt.rcParams.update({'font.size' : 22})

axisAni = [0.0, 0.25, 0.38, 0.41, 0.50, 0.51, 0.52, 0.58, 0.61, 0.65, 0.68]
axisAni_Err = [0.01, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]

axisExch = [0.0, 0.25, 0.32, 0.35, 0.38, 0.42, 0.46, 0.47, 0.50, 0.50, 0.50]
axisExch_Err = [0.01, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04]

gamma = np.arange(0.0, 1.1, 0.1)


fig, ax = plt.subplots(figsize=(8, 6))

ax.grid()

#ax.hlines(y=2.269, xmin=0, xmax = 1, linestyle='dashed')

ax.errorbar(gamma, axisAni, yerr=axisAni_Err,
			color = 'blue',
				   capsize=5,
				   marker='o',
				   markerfacecolor = 'g',
				   label='Axis Anisotropy')

ax.errorbar(gamma, axisExch, yerr=axisExch_Err,
			color = 'green',
				   capsize=5,
				   marker='o',
				   markerfacecolor = 'g', label='Exchange Anisotropy')

ax.legend()


ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='minor', length=3, color='black', width = 0.8)
ax.tick_params(which='major', length=5, color='black', width = 1.2)

ax.set_title(r"Variation of $\mathrm{T_c}$ due to $\gamma$ at fixed J")

ax.set_ylabel(r'Critical Temperature $T_C$')
ax.set_xlabel(r'$\gamma$')