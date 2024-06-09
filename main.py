'''
main function to run the simulation
#
-> pos = position of each atom ; vel = velocity of each atom ;
-> force = forces on each atom ; mass = mass of each atom ;
-> n = number of atoms ; dt = time step ; Kb = Boltzmann constant ;
-> sigma = sigma of the LJ potential ; epsilon = epsilon of the LJ potential ;
#
Default: READ XYZ-file FOR INTIAL POSITIONS (Units -> Angstrom)
Default: epsilon = 0.20834          # Electron Volt (eV)
Default: sigma = 2.2209             # Angstrom
Default: target_T = 0.0005          # Kelvin
Default: box = [10.0, 10.0, 10.0]   # Angstrom
Default: time_step = 20.0           # fs
Default: Kb = 1.0                   # atomic unit
Default: langevin = False 
Default: v_rescale = False
'''
import numpy as np
import sys
import matplotlib.pyplot as plt
from LJengo.dynamics import run_dynamics
##############################################################################
input_file = sys.argv[1]
##############################################################################
with open(input_file, 'r') as f:
    lines = f.readlines()
    natom = int(lines[0])
    dummy = str(lines[1].strip())
    pos = np.zeros([natom, 3])
    atom = []
    for i in range (0, natom):
        cord = lines[2 + i].split()
        atom_type = cord[0]
        atom.append(atom_type)
        pos[i, 0] = float(cord[1])
        pos[i, 1] = float(cord[2])
        pos[i, 2] = float(cord[3])
##############################################################################
# DO NOT MODIFY ABOVE THIS, UNLESS YOU WANT TO CHANGE THE INPUT FILE FORMAT
##############################################################################
time_step = 20.0 ; kB = 1.0
step = int(5e3) ; Temp = 0.0005
box = [10.0, 10.0, 10.0]
ang_to_au = 1.0/0.529177
amu_to_au = 1840.0
epsilon = 0.20834   # electron volt (eV)
sigma = 2.2209    # angstrom
#langevin = False
#v_rescale = False
mass = 40*np.ones(natom)*amu_to_au
##############################################################################
run_dynamics(step=step, T=Temp, n=natom, mass=mass, box=box, eps=epsilon, sig=sigma, 
        dt=time_step, Kb=kB, pos=pos, dyn_type='NVE')
# run_dynamics(T=Temp,dt=time_step, pos=pos) #, dyn_type='NVT: langevin')



##############################################################################
# # Time, Kinetic Energy, Potential Energy, Total Energy, Temperature
data = np.loadtxt('energy.dat', delimiter=',')
pe_data = np.loadtxt('potential.dat', delimiter=',')

i_step, KE, PE, TE, Temp  = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
dist, pot = pe_data[:,0], pe_data[:,1]
#

