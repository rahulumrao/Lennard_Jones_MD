'''
# main Molecular Dynamics loop
# 
The Lennard-Jones potential is given by the equation:

V(r) = 4 * ε * [(σ / r)^12 - (σ / r)^6]

Where:
- V(r) is the potential energy between two particles at distance r.
- ε is the depth of the potential energy well.
- σ is the finite distance at which the inter-particle potential is zero.
- r is the distance between the particles.

Here I am using A, B scheme for LJ potential.
Where:
-> A = 4 * ε * σ ^12 
-> B = 4 * ε * σ ^12 
#
-> T = temp ; n = number of atoms ; mass = mass of each atom ; 
-> # box = box size ; eps = epsilon ; sig = sigma ;
-> dt = time step ; Kb = Boltzmann constant ;
-> pos = position of each atom ; vel = velocity of each atom ;
-> force = forces on each atom ;
-> R = distance between two atoms ; pot = potential energy ;
-> ke = kinetic energy ; T = temperature ;
#
# Velocity-Verlet integration
# update_positions()
# update_velocities()
# calculate_forces()
# update_velocities()
'''
import numpy as np
import warnings
import inspect
from tqdm import tqdm
from .initial_condition import initial_vel
from .force import compute_force
from .update import update_vel, update_pos
from .thermostat import compute_ke, vel_rescale, langevin_dynamics
#####################################################################
def run_dynamics(step=100, 
                 T=0.0005,
                 n=2,
                 mass=[40,40],
                 box=[10.0, 10.0, 10.0],
                 eps=0.20834,
                 sig=2.2209,
                 dt=20.0,
                 Kb=1.0,
                 pos=None,
                 dyn_type=None):

    # Preprocess input based on the dynamics type (NVT or NVE)
    if dyn_type not in ["NVT: langevin", "NVT: v_rescale", "NVE"]: dyn_type = "NVE"
    if dyn_type == "NVT: langevin":
        langevin = True
        v_rescale = False
        print("# NVT Dynamics (Langevin).!")
    elif dyn_type == "NVT: v_rescale":
        langevin = False
        v_rescale = True
        print("#  NVT Dynamics (Velocity re-scale).!")
    elif dyn_type == "NVE" or dyn_type == None:
        dyn_type = "NVE"
        langevin = False
        v_rescale = False
        print("# Running NVE Dynamics.!")
    #        
    # Print the keyword arguments and their values as a table
    print("Keyword arguments:")
    print("------------------")
    print(f"MD STEPS             : {step}")
    print(f"Temperature          : {T}")
    print(f"Number of atoms      : {n}")
    print(f"Box                  : {box}")
    print(f"Epsilon              : {eps}")
    print(f"Sigma                : {sig}")
    print(f"MD Time Step         : {dt}")
    print(f"Boltzmann Constant   : {Kb}")
    print(f"Dynamics Type        : {dyn_type}")
    #
    # constants
    ang_to_au = 1.0/0.529177
    amu_to_au = 1840.0
    mass = 40*np.ones(n)*amu_to_au
    target_T = T
    box = np.array(box) * ang_to_au
    # box = [x * ang_to_au for x in box]  # convert: units from angstrom to a.u_
    #
    with open('traj.xyz', 'w') as f:
        f.write('')

    with open('energy.dat', 'w') as enfile:
        enfile.write('# Time, Kinetic Energy, Potential Energy, Total Energy, Temperature\n')

    with open('potential.dat', 'w') as ptfile:
        ptfile.write('# distance, Potential Energy\n')

    # Default: initial positions and velocities for 2 atoms
    # if pos is None: pos = [[0.0, 0.0, 0.0],[2.8, 0.0, 0.0]]
    vel = [[1.1088e-04, 1.1763e-05, 3.3963e-05],
        [6.1605e-05, 6.2051e-05, 7.7081e-05]]
    
    # Convert positions to atomic units
    # pos = np.array(pos)
    pos = pos * ang_to_au
    
    # Initialize velocities
    vel = initial_vel(T, vel, mass, n, Kb)

    # Compute intial force and potential energy
    R, pot, force = compute_force(pos, box, n, eps, sig)

    # main MD loop
    # for md in range (0, step):
    for md in tqdm(range(step), desc='Running MD simulation', unit="step"):

        pos = update_pos(pos, force, vel, dt, n, mass, box)
        vel = update_vel(vel, force, mass, dt, n)
        R, pot, force = compute_force(pos, box, n, eps, sig)
        if (langevin == True): 
            # print(f'Thermostat Called and Temp is {T:8.5f}')
            force = langevin_dynamics(force, vel, mass, n, target_T, dt, Kb, gamma=0.1)
        if (v_rescale == True):
            # print(f'Thermostat Called and Temp is {T:8.5f}')
            vel = vel_rescale(vel, mass, n, dt, target_T, Kb)
        vel = update_vel(vel, force, mass, dt, n)

        # force = np.array(F)
        ke , T = compute_ke(vel, mass, n, Kb)

        if (md % 100 == 0):
            with open('traj.xyz', 'a') as f:
                f.write(f'{n}\n')
                f.write(f'Lennard Jones MD\n')
                for i in range (0, n):
                    f.write(f"Ar {(pos[i,0]/ang_to_au):12.6f} {(pos[i,1]/ang_to_au):12.6f} {(pos[i,2]/ang_to_au):12.6f}\n")
            with open('energy.dat', 'a') as enfile:
                enfile.write(f'{md:6d}, {ke:12.6f}, {pot:12.6f}, {(ke+pot):12.6f}, {(T):12.5f}\n')
            with open('potential.dat', 'a') as ptfile:
                ptfile.write(f'{(R*0.529177):8.4f}, {pot:12.6f}\n')
############################ End Program #########################################
