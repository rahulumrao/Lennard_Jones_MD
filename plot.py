#!/bin/bash

import numpy as np
import matplotlib.pyplot as plt
##############################################################################
# # Time, Kinetic Energy, Potential Energy, Total Energy, Temperature
# data = np.loadtxt('energy.dat', delimiter=',')
pe_data = np.loadtxt('potential.dat', delimiter=',')

dist, pe = pe_data[:,0], pe_data[:,1]
#
plt.plot(dist, pe, c='r', lw=2)
# axis ticks grid
plt.xticks(np.arange(1.5, max(dist), 0.5))
# axis tick size
plt.tick_params(axis='x', labelsize=16)
plt.tick_params(axis='y', labelsize=16)

# axis minimum and maximum limit
plt.xlim([1.5, None])
#
plt.grid(color='grey', linestyle='-.', linewidth='1.0')
# set axis label
plt.xlabel(r' R ${\rm{(\AA)}}$', fontsize='14')
plt.ylabel(r'Potential Energy (a.u)', fontsize='14')
# saving plot
plt.savefig('PE_plot.jpeg', bbox_inches='tight', transparent=True, dpi=300)
