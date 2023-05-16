# wCDM Training Simulations - Input Files
Necessary files to run wCDM training simulations. They include:
- The Latin Hypercube Sample of cosmological parameters
- The transfer functions for each cosmology
- `Transferinfo` files, i.e. files that tell COLA where the transfer functions are located and at which redshifts they are calculated.
- Lua files, which are the basic COLA inputs containing all simulation settings (e.g. cosmology, initial condition generation, number of particles, resolution settings...)


## Generating Transfer Functions
COLA needs the transfer functions of all species in the N-body gauge. They are calculated using Gui's modified [`hi_class`](https://github.com/SBU-Jonathan/hi_class_Nbody). Make 
sure to install it before running the generator script, `transfer_function_generator.py`.

The script reads the LHS, calculates transfer functions for each cosmology and saves them in `./transfer_functions/{i}/`, where `i` is the index (zero-based, going from 0 to 499) 
of the cosmology in the LHS. The script also saves a `transferinfo` file: a file which contains the names of the transfer function files and to which redshift they refer to. These 
`transferinfo` files are needed to run the simulations.

Additionally, the script will generate a `cosmo_params/cosmo_{i}.txt` file. CLASS takes as input `Omega_cdm`, whereas our LHS samples over `Omega_m`, the total matter density, 
defined as `Omega_cdm + Omega_b + Omega_massivenu`. The problem is that we don't fix `Omega_massivenu`, but `m_nu`. We don't know `Omega_massivenu` in advance. The workaround is: 
using an approximation `Omega_nu_approx * h^2 = m_nu / (93.14 eV)` to calculate `Omega_cdm = Omega_m - Omega_b - Omega_nu_approx`, run CLASS, obtain the exact value of 
`Omega_massivenu` from CLASS and save it (and other cosmological parameters) to this file. These cosmological parameters are then written in the Lua file.

**NOTE**: Before running the script, open it to adjust a few paths:
- `amypond_direc`: the path to save the transfer functions in the current system
- `cluster_path`

To generate the transfer functions from cosmology `NMIN` to cosmology `NMAX`, simply 
run the command:

```
$ python transfer_function_generator.py NMIN NMAX
```

**NOTE**: Since CLASS is being run with high precision, it takes around 1min (in AmyPond) to generate transfers for one cosmology.
**NOTE**: The script seems to be leaking memory: each subsequent transfer function calculation occupies more and more memory in the system. In AmyPond, I recommend running only 5 
or 6 cosmologies at time. I thought Python wasn't supposed to leak memory :(

The script
