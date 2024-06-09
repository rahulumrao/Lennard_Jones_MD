'''
If NVT is activated, the temperature of the system is kept constant.
Langevin thermostat, or Velocity rescaling is used to control 
the temperature of the system.
#
-> pos = position of each atom ; vel = velocity of each atom ;
-> force = forces on each atom ; mass = mass of each atom ;
-> n = number of atoms ; dt = time step ; Kb = Boltzmann constant ;
-> target_T = target temperature of the system ;
-> langevin = True/False ; v_rescale = True/False ;
-> gamma = friction coefficient ;
#
'''
import numpy as np
#####################################################################
def compute_ke(vel, mass, n, Kb=1.0):
    # computing kinetic energies and temperature
    ke = 0.0
    temp = 0.0
    for i in range (0, n):
        ke += (0.5*mass[i]*(vel[i,0]**2 + vel[i,1]**2 + vel[i,2]**2))
    temp = (2.0 * ke)/(3.0 *n * Kb)
    return ke, temp
#####################################################################
def vel_rescale(vel, mass, n, dt, temp, Kb=1.0):
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

    return vel
#####################################################################
def langevin_dynamics(force, vel, mass, n, temp, dt, Kb=1.0, gamma=0.01):
    '''
    Langevin modified forces to control system temperature. 
    This consists of adding a friction term to the force.
    The random force is calculated such that the average force is zero.
    gamma : friction coefficient
    '''
    for i in range (0, n):
        for j in range (3):
            r = np.random.randn()
            r -= 0.5
            # r = np.random.uniform(-0.5, 0.5)
            gamma2 = 2.0 * gamma * Kb * temp
            random_force =  np.sqrt((gamma2 * mass[i]) / dt )          
            noise = random_force * r
            friction = (- gamma * vel[i, j] * mass[i])
            force[i, j] += (friction + noise)

    return force
#####################################################################
