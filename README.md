# wCDM Training Simulations - Input Files
Necessary files to run COLA wCDM training simulations. They include:
- The Latin Hypercube Sample of cosmological parameters
- The transfer functions for each cosmology
- `Transferinfo` files, i.e. files that tell COLA where the transfer functions are located and at which redshifts they are calculated. These files are passed as arguments in the lua files mentioned below.
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
`Omega_massivenu` from CLASS and save it (and other cosmological parameters) to this file. These cosmological parameters are then written in the Lua file. Once the lua files are generated, the `cosmo_{i}.txt` files are not needed anymore.

**NOTE**: Before running the script, open it to adjust a few paths:
- `amypond_direc`: the path to save the transfer functions in the current system
- `cluster_path`: the path where the transfer functions will be in the cluster that will run the simulations. I currently use `/scratch/decola/joao.reboucas2/COLA_projects/wCDM/transfer_functions/{i}`, where `i` is the index of the cosmology in the LHS.
- `transferinfo_path`: the path to save the `transferinfo` file in the current system. In the `transferinfo` file, the `cluster_path` will be written, telling COLA where the transfer functions are. I use `/scratch/decola/joao.reboucas2/COLA_projects/wCDM/transferinfo_files/transferinfo{i}.txt`
Adapt these paths to your cluster.

To generate the transfer functions from cosmology `NMIN` to cosmology `NMAX`, simply 
run the command:

```
$ python transfer_function_generator.py NMIN NMAX
```

**NOTE**: Since CLASS is being run with high precision, it takes around 1min (in AmyPond) to generate transfers for one cosmology.

**NOTE**: The script seems to be leaking memory: each subsequent transfer function calculation occupies more and more memory in the system. In AmyPond, I recommend running only 5 
or 6 cosmologies at time. I thought Python wasn't supposed to leak memory :(

**NOTE**: my workaround for this problem is to use a bash script to sequentially generate transfer functions in batches of 5 or 6. The script should be something like
```
#!/bin/bash
python transfer_function_generator.py 0 5
python transfer_function_generator.py 6 11
...
```
and it can be executed, for instance, using the `nohup` command. If you use this approach, keep in mind that cancelling the commands is not trivial, since `ctrl+c` will only cancel the current command 
and move to the next one.

Once you run the script for all cosmologies, you should have the transfer functions and transferinfo files.

**NOTE**: The transfer function generator script can also generate files for three test wCDM cosmologies, where the cosmological parameter are set to the EE2 reference values and w varies between -1.3, 
-0.7 and -0.9999 (these runs are marked as `test_0`, `test_1` and `test_2` respectively). To generate the files for test simulations, run the command:
```
$ python transfer_function_generator.py 0 2 --test
```

## Generating lua files
The lua file generation is much simpler. It is done via the notebook `lua_script_generator.ipynb`. It finds the cosmological parameters in `cosmo_params` and writes them in the lua.

Apart from the cosmological parameters, the lua files also contain:
- The path where the transferinfo files will be located. I use `/scratch/decola/joao.reboucas2/COLA_projects/wCDM/transferinfo_files/transferinfo{i}.txt`. Adapt this path for your cluster.
- The output path where the power spectra will be stored during the simulation. I use `/scratch/decola/joao.reboucas2/COLA_projects/wCDM/outputs/{i}`.
- The path where the lua files will be saved in the system.
- The box size of the simulation
- Number of particles in 1D (such that `N_total = Np1d^3`)
- Force resolution
- Power spectra resolution
- Whether to reverse phases (important for pair-fixing)

To generate the lua files for the LHS runs, just run the first three cells of the notebook. To generate lua files for the test runs described above, run the fourth cell.

## Running the simulations
Install [`FML`](https://github.com/SBU-Jonathan/FML_AUGUST_2020) in your cluster. I will assume the following structure:
```
FML/ (path where the simulations need to be submitted)
COLA_projects/
---- wCDM/
     ---- transfer_functions/
     ---- transferinfo_files/
     ---- lua_files/
     ---- outputs/
```
but `COLA_projects` can also be inside `FML` if the paths are correctly adjusted. The COLA binary is inside `FML/FML/COLASolver/nbody` (there in an FML directory inside the FML directory, be careful). The command to run COLA is (assuming you are in the root directory of FML where `start_cola` is):
```
$ ./FML/COLASolver/nbody ../COLA_projects/wCDM/lua_files/run_${SLURM_JOB_ARRAY_ID}.lua
```

An example script that works in Santos Dumont is given in `scripts/`. Adapt it to your cluster.
