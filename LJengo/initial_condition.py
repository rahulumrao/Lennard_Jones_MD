'''
Intialiaze velocity from a Boltzmann distribution
'''
import numpy as np
import math as m
#####################################################################
def initial_vel(T, vel, mass, n, Kb=1.0):
    '''
    intialize velocity from a Boltzmann distribution
    '''
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
        #print(f'\n Velocities:,{vel[i,0], vel[i,1], vel[i,2]}, Temperature:,{temp}')
    return vel
#####################################################################