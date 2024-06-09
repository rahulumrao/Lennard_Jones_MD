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

## How to
### Download :
```bash
git clone https://github.com/rahulumrao/Lennard_Jones_MD`
```

### Install :
```bash
cd Lennard_Jones_MD

pip install .
```

## How to use:
_After installation, use the following command:_ <br>

```bash
from LJengo import run_dynamics
```

_use of `run_dynamics` :_ <br>

```bash
run_dynamics(pos = 'cartesian coordinates of each atom (dimension: N x 3)', 
                 step = 'number of molecular dynamics steps', 
                 T = 'temperature',
                 n = 'number of atoms',
                 mass = 'atomic mass of each atom',
                 box = 'box size',
                 eps = 'LJ - epsilon',
                 sig = 'LJ - sigma',
                 dt = 'md time step',
                 Kb = 'Boltzmann constant',
                 dyn_type = 'type of dynamics')
```

### Run:

``` bash
python main.py input.xyz
```

## Arguments:
-> T = temp                 (float)     <br>
-> n = number of atoms      (integer)   <br>
-> eps = epsilon            (float)     <br>
-> sig = sigma              (float)     <br>
-> dt = time step           (float)     <br>
-> Kb = Boltzmann constant  (float)     <br>
-> mass = mass of each atom (np.array [dimension: N])               <br>
-> box = box size       (assuming Orhorhomic box [dimension: 3])    <br>
-> pos = position of each atom (np.array [dimension: N x 3])        <br>
-> vel = velocity of each atom (np.array [dimension: N x 3])        <br>
-> force = forces on each atom (np.array [dimension: N x 3])        <br>
-> dn_type = MD type (str [NVE or NVT])


# Full code in a single `Python` file -> `lj_run.py`

_run_:

```bash
python lj_run.py input.xyz
```