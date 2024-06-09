'''
updating positions and velocities of each atom using Velocity-Verlet integration
#
-> pos = position of each atom ; force = forces on each atom ;
-> vel = velocity of each atom ; dt = time step ; n = number of atoms ;
-> mass = mass of each atom ; box = size of the simulation box
#
update_pos(pos, force, vel, dt, n, mass, box)
#
'''
import numpy as np
################################################################
def update_pos(pos, force, vel, dt, n, mass, box):
    # updating positions of each n-atom
    for i in range (0, n):
        for j in range (3):
            # apply periodic boundary conditions using Minimum Image Convention
            if (pos[i, j] < (box[j] * 0.5)): 
                pos[i, j] += box[j]
            if (pos[i, j] >= (box[j] * 0.5)): 
                pos[i, j] -= box[j]

        # Position Update
        for j in range (3):
            pos[i, j] += (dt*vel[i,j]) + (0.5 * dt * dt * (force[i,j] / mass[i]))

        # # apply periodic boundary conditions using Minimum Image Convention
        # for j in range (3):
        #     if (pos[i, j] < -box[j] * 0.5): 
        #         pos[i, j] += box[j]
        #     elif (pos[i, j] >= box[j] * 0.5): 
        #         pos[i, j] -= box[j]
            
    return pos
################################################################
def update_vel(vel, force, mass, dt, n):
    # updating velocities of each n-atom
    for i in range (0, n):
        for j in range (3):
            vel[i, j] += (0.5 * dt * (force[i,j]) / mass[i])

    return vel
################################################################