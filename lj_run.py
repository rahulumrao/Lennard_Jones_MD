#!/usr/bin/python
import os
import sys
import numpy as np
import math as m
import matplotlib.pyplot as plt
#

def compute_force(pos, box, n, eps, sigma):
    '''
    # computing forces on each atom from cartesian positions
    # Using A, B scheme of Lennard Jones
    A = 12000.0 ; B = 100.0
    '''
    A = 4*eps*sigma**12
    B = 4*eps*sigma**6
    ev_au  = 27.21139
    ang_au = 0.529177
    A = A/(ev_au*(ang_au**12))      # converting to the units from eV to a.u_
    B = B/(ev_au*(ang_au**6))      # converting to the units from eV to a.u_
    force = np.zeros([n, 3])      # creating a numpy array of shape (atom_num , 3)
    pot = 0

    for i in range (0, n-1):
        for j in range (i+1, n):
            dx, dy, dz = 0.0, 0.0, 0.0 
            R = 0.0
            
            # XYZ components of i and j atom pairs
            dx = pos[i, 0] - pos[j, 0]
            dy = pos[i, 1] - pos[j, 1]
            dz = pos[i, 2] - pos[j, 2]
            
            # apply periodic boundary conditions
            dx = dx - box[0] * round(dx/box[0])
            dy = dy - box[1] * round(dy/box[1])
            dz = dz - box[2] * round(dz/box[2])
            
            # calculate distance between i and j atom pairs
            R = np.sqrt(dx*dx + dy*dy + dz*dz)

            # compute potential energy
            pot = pot + ((A/R**12) - (B/R**6))

            # compute force
            fx = -((-12*A/R**13) - (-6*B/R**7))*(dx/R)
            fy = -((-12*A/R**13) - (-6*B/R**7))*(dy/R)
            fz = -((-12*A/R**13) - (-6*B/R**7))*(dz/R)

            # appending force components to ith element
            force[i, 0] = force[i, 0] + fx
            force[i, 1] = force[i, 1] + fy
            force[i, 2] = force[i, 2] + fz

            # appending force components to jth element
            force[j, 0] = force[j, 0] - fx
            force[j, 1] = force[j, 1] - fy
            force[j, 2] = force[j, 2] - fz
        
        print(f'dist: {R*0.529177} pot: {pot}\n')
    return R, pot, force

################################

def update_pos(pos, force, vel, dt, n, mass, box):
    # updating positions of each n-atom
    for i in range (0, n):
        for j in range (3):
            pos[i, j] = pos[i, j] + (dt*vel[i,j]) + (0.5 * dt * dt * (force[i,j] / mass[i]))

        # apply periodic boundary conditions using Minimum Image Convention
        for j in range (3):
            if (pos[i, j] < -box[j] * 0.5): 
                pos[i, j] += box[j]
            elif (pos[i, j] >= box[j] * 0.5): 
                pos[i, j] -= box[j]
            
    return pos

################################

def update_vel(vel, force, mass, dt, n):
    # updating velocities of each n-atom
    for i in range (0, n):
        for j in range (3):
            vel[i, j] = vel[i, j] + (0.5 * dt * (force[i,j]) / mass[i])

    return vel
################################
def initial_vel(T, vel, mass, n, Kb=1.0):
    import math as m
    # intialize velocity from a Boltzmann distribution
    vel = np.zeros([n, 3])
    vel = np.random.rand(n, 3) # intialize random velocity for each atom
    for i in range (0, n):
        vel[i, 0] = vel[i,0]*(m.sqrt(mass[i]/(2.0*m.pi*Kb*T)))
        vel[i, 1] = vel[i,1]*(m.sqrt(mass[i]/(2.0*m.pi*Kb*T)))
        vel[i, 2] = vel[i,2]*(m.sqrt(mass[i]/(2.0*m.pi*Kb*T)))

        temp = (0.5*mass[i]*(vel[i,0]**2 + vel[i,1]**2 + vel[i,2]**2))/Kb
        # scale velocity to the target temperature T
        vel[i,0] = vel[i,0]*(m.sqrt(T/temp))
        vel[i,1] = vel[i,1]*(m.sqrt(T/temp))
        vel[i,2] = vel[i,2]*(m.sqrt(T/temp))

        temp = (0.5*mass[i]*(vel[i,0]**2 + vel[i,1]**2 + vel[i,2]**2))/Kb
        print(f'\n Velocities:,{vel[i,0], vel[i,1], vel[i,2]}, Temperature:,{temp}')
    return vel
################################
def compute_ke(vel, mass, n, Kb=1.0):
    # computing kinetic energies and temperature
    ke = 0.0
    temp = 0.0
    for i in range (0, n):
        ke += (0.5*mass[i]*(vel[i,0]**2 + vel[i,1]**2 + vel[i,2]**2))
    temp = (2.0 * ke)/(3.0 *n * Kb)
    return ke, temp
################################
def vel_rescale(vel, mass, n, temp, Kb=1.0):
    '''
    # velocity rescaling
    '''
    _ , temp_i = compute_ke(vel, mass, n, Kb)

    scale = np.sqrt(temp/temp_i)
    tau = 200 * dt  # Relaxation time constant for smoothing, adjust as needed

    # Calculate rescaled velocities
    rescaled_vel = vel * scale
    
    # Smooth the rescaling over time
    vel += (rescaled_vel - vel) * (1 - np.exp(-dt / tau))
    # vel = vel * scale
    # print(f'scale : {scale}')
    return vel

################################
def langevin_dynamics(force, vel, mass, n, temp, dt, Kb=1.0, gamma=0.01):
    '''
    Langevin modified forces to control system temperature. 
    This consists of adding a friction term to the force.
    The random force is calculated such that the average force is zero
    gamma : friction coefficient
    '''
    for i in range (0, n):
        for j in range (3):
            # r = np.random.randn()
            # r -= 0.5
            r = np.random.uniform(-0.5, 0.5)
            gamma2 = 2.0 * gamma * Kb * temp
            random_force =  np.sqrt((gamma2 * mass[i]) / dt )          
            noise = random_force * r
            friction = (- gamma * vel[i, j] * mass[i])
            force[i, j] += (friction + noise)
            #force[i, j] += friction + random_force
            # force[i,j] += np.sqrt(2 * gamma * Kb * temp / mass[i]) * r
    return force
################################
input_file = sys.argv[1]

with open(input_file, 'r') as f:
    lines = f.readlines()
    n = int(lines[0])
    dummy = str(lines[1].strip())
    pos = np.zeros([n, 3])
    atom = []
    for i in range (0, n):
        cord = lines[2 + i].split()
        atom_type = cord[0]
        atom.append(atom_type)
        pos[i, 0] = float(cord[1])
        pos[i, 1] = float(cord[2])
        pos[i, 2] = float(cord[3])
################################
print(f'num of atoms: {n}')
dt = 20.0 ; Kb = 1.0
step = 10000 ; T = 0.0005
target_T = T
ang_to_au = 1.0/0.529177
amu_to_au = 1840.0
box = [10.0, 10.0, 10.0]
vel = np.zeros([n, 3])
eps = 0.20834   # electron volt (eV)
sig = 2.2209    # angstrom
langevin = False 
v_rescale = False 
mass = 40*np.ones(n)*amu_to_au
vel = initial_vel(T, vel, mass, n, Kb)
#converting to the units from Angstrom to a.u.
box = [x * ang_to_au for x in box] 
######################################################################################################################
with open('traj.xyz', 'w') as f:
    f.write('')

with open('energy.dat', 'w') as enfile:
    enfile.write('# Time, Kinetic Energy, Potential Energy, Total Energy, Temperature\n')

with open('potential.dat', 'w') as ptfile:
    ptfile.write('# distance, Potential Energy\n')

pos = pos * ang_to_au
R, pot, force = compute_force(pos, box, n, eps, sig)
for md in range (0, step):

    pos = update_pos(pos, force, vel, dt, n, mass, box)
    vel = update_vel(vel, force, mass, dt, n)
    R, pot, force = compute_force(pos, box, n, eps, sig)
    if (langevin == True): 
        print(f'Thermostat Called and Temp is {T:8.5f}')
        force = langevin_dynamics(force, vel, mass, n, target_T, dt, Kb, gamma=0.1)
    if (v_rescale == True):
        print(f'Thermostat Called and Temp is {T:8.5f}')
        vel = vel_rescale(vel, mass, n, target_T, Kb)
    vel = update_vel(vel, force, mass, dt, n)

    # force = np.array(F)
    ke , T = compute_ke(vel, mass, n, Kb)

    with open('traj.xyz', 'a') as f:
        f.write(f'{n}\n')
        f.write(f'Lennard Jones MD\n')
        for i in range (0, n):
            f.write(f"Ar {(pos[i,0]/ang_to_au):12.6f} {(pos[i,1]/ang_to_au):12.6f} {(pos[i,2]/ang_to_au):12.6f}\n")
    with open('energy.dat', 'a') as enfile:
        enfile.write(f'{md:6d}, {ke:12.6f}, {pot:12.6f}, {(ke+pot):12.6f}, {(T):12.5f}\n')
    with open('potential.dat', 'a') as ptfile:
        ptfile.write(f'{(R*0.529177):8.4f}, {pot:12.6f}\n')

######################################################################################################################
