# Lennard-Jones Molecular Dynamics
```bash
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
```

## How to install
### download the package -
` git clone https://github.com/`

` pip install .`

## How to use:
_After instantiation the package use the following command:_ <br>

**from LJengo import run_dynamics**

_The use of this command:_ <br>

```bash
run_dynamics(pos = 'cartesian coordinates of each atom (dimension: N x 3)', 
                 step=100, 
                 T=0.0005,
                 n=2,
                 mass=[40,40],
                 box=[10.0, 10.0, 10.0],
                 eps=0.20834,
                 sig=2.2209,
                 dt=20.0,
                 Kb=1.0,
                 pos=None,
                 dyn_type=None)
```

# Arguments:
-> T = temp                 (float)
-> n = number of atoms      (integer)
-> eps = epsilon            (float)
-> sig = sigma              (float)
-> dt = time step           (float)
-> Kb = Boltzmann constant  (float)
-> mass = mass of each atom (np.array [dimension: N]) 
-> box = box size       (assuming Orhorhomic box [dimension: 3])
-> pos = position of each atom (np.array [dimension: N x 3])
-> vel = velocity of each atom (np.array [dimension: N x 3])
-> force = forces on each atom (np.array [dimension: N x 3])
