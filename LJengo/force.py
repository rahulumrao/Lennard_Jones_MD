'''
The Lennard-Jones potential is given by the equation:

V(r) = 4 * ε * [(σ / r)^12 - (σ / r)^6]

Where:
- V(r) is the potential energy between two particles at distance r.
- ε is the depth of the potential energy well.
- σ is the finite distance at which the inter-particle potential is zero.
- r is the distance between the particles.

<a/> https://en.wikipedia.org/wiki/Lennard-Jones_potential

Here I am using A, B scheme for LJ potential.
Where:
-> A = 4 * ε * σ ^12 
-> B = 4 * ε * σ ^12 
#
-> pos = position of each atom ; vel = velocity of each atom ;
-> force = forces on each atom ; mass = mass of each atom ;
-> box = box size ; eps = epsilon ; sig = sigma ;
-> dt = time step ; n = number of atoms ;
-> Kb = Boltzmann constant ;
#
'''
#####################################################################
import numpy as np
def compute_force(pos, box, n, eps, sigma):

    # Using A, B scheme of Lennard Jones
    # A = 12000.0 ; B = 100.0

    A = 4*eps*sigma**12
    B = 4*eps*sigma**6
    ev_au  = 27.21139
    ang_au = 0.529177
    A = A/(ev_au*(ang_au**12))      # convert: the units from eV to a.u_
    B = B/(ev_au*(ang_au**6))      # convert:  the units from eV to a.u_
    box = [x / 0.529177 for x in box] #convert: units from Angstrom to a.u.
    force = np.zeros([n, 3])      # create: an empty numpy array of shape (atom_num , 3)
    pot = 0
    #
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
            force[i, 0] += fx
            force[i, 1] += fy
            force[i, 2] += fz

            # appending force components to jth element
            force[j, 0] -= fx
            force[j, 1] -= fy
            force[j, 2] -= fz
        
        # print(f'dist: {R*0.529177} pot: {pot}\n')
    return R, pot, force
#####################################################################
