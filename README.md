# About
Tutorial for using the usage of the dcTMD approach and Langevin simulations using the trypsin-benzamidine complex (see [Wolf et al., Nat. Commun. 2020, 11, 2918](https://www.nature.com/articles/s41467-020-16655-1)) as example.
The basis of this tutorial are constraint pulling simulations as described in [Wolf and Stock, JCTC 2018, 14, 6175](https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00835). For pathway separation (advanced topic, see below), we recommend the usage of a principal component analysis (PCA) as possible via our [fastPCA code](https://github.com/moldyn/FastPCA).

# Licensing

If you perform a dcTMD analysis of your own data for published works, please cite the appropriate literature:

- for the dcTMD analysis itself: _Wolf, S., & Stock, G. (2018). Targeted Molecular Dynamics Calculations of Free Energy Profiles Using a Nonequilibrium Friction Correction. Journal of Chemical Theory and Computation, 14, 6175–6182._
- for temperature-boosted Langevin equation simulations and/or in case you need to perform a pathway separation: _Wolf, S., Lickert, B., Bray, S., & Stock, G. (2020). Multisecond ligand dissociation dynamics from atomistic simulations. Nature Communications, 11(1), 2918. http://doi.org/10.1038/s41467-020-16655-1_


# Tutorial

## Getting started

You will need the Python analysis scripts from [dcTMD](https://github.com/floWneffetS/dcTMD) as well as the C++ Langevin equation integrator code and Jupyter notebooks from [here](https://github.com/floWneffetS/Langevin_T_boost). Simulations need to be carried out in [Gromacs v2018 and higher](https://manual.gromacs.org/documentation/). Older versions support constraint pulling, as well, but you need to adjust the syntax in the MDP files accordingly. 

**Note:** You may also use other simulation MD software packages for the simulations that support constraint pulling (e.g. CHARMM). The dcTMD analysis only requires an ASCII file with two colums: time in the first column, constraint forces in the second column. Please be advised that the analysis scripts all assume Gromacs-internal units, i.e., times in ps and forces in kJ / (mol nm)

In the following, we assume that you are familiar with the general usage of Gromacs. Good tutorials on how to use the program can be found [here](http://www.mdtutorials.com/gmx/).

## TMD simulation

### Input files
Download and unpack tutorial_files.tar.gz via
```
tar -xzvf ./tutorial_files.tar
```
You will find a folder with the following files:

Files for run input commands:
`3ptb_AMBER99SB_ben_pushEQUIBRUN.mdp 
3ptb_AMBER99SB_ben_pushRUN_v0.001.mdp
`
- The `pushEQUIBRUN.mdp` file is an initial equilibration file for generating start simulation files with different initial velocity distributions. 
- The `pushRUN_v0.001.mdp` file is the respective command input for the non-equilibrium pulling. 
- `3ptb` refers to the protein data base code of a trypsin structure, `AMBER99SB` to the employed force field, and `ben` to the benzamidine ligand.

Structure file: `3ptb_AMBER99SB_ben.gro` in which the trypsin-bezamidine complex is equilibrated in TIP3P water with a physiological NaCl concentration and a single Ca2+ ion.

Topologies and position restraint files:
`3ptb_AMBER99SB_ben.top
3ptb_AMBER99SB_ben_Protein_chain_A.itp
3ptb_AMBER99SB_ben_Ion_chain_B.itp
3ptb_ben_H2_GMX_RESP.itp
posre_Protein_chain_A.itp
posre_Ion_chain_B.itp
posre_ben.itp
`

Index file:
`3ptb_AMBER99SB_ben.ndx` 
**Important:** the index file needs to include an anchor group from whose center of mass the ligand is pulled away (in this case: the group `[sheet]` containing C-alpha atoms from the central beta-sheet) and the ligand itself (or better, the heavy atoms of the ligand, here group `[ BEN_heavy ]`). If you want to create a respective anchor index for your own simulation problem, choose an anchor group that is tighly connected to the remainder of the protein (such as C-alpha atoms in alpha-helices and beta-sheets). The vector connecting the centers of mass of anchor and ligand needs to roughly point into the direction of a putative unbinding path.


### Carrying out pulling MD simulations

For the generation of the input structure in your own project, we advise you to carry out an initial NPT equilibration simulation of at least 10 ns length. Here, we have done this already for you and generated an equilibrated structure.

You will need a number (between 100-200) of equilibrated trajectories with different initial velocity distributions. For this, generate an initial equilibration folder and the simulation start TPR files using:
```
mkdir equib
cd equib/

for i in {000..099}
do
gmx grompp -f ../3ptb_AMBER99SB_ben_pushEQUIBRUN.mdp -c ../3ptb_AMBER99SB_ben.gro -r ../3ptb_AMBER99SB_ben.gro -p ../3ptb_AMBER99SB_ben.top -n ../3ptb_AMBER99SB_ben.ndx -o 3ptb_AMBER99SB_ben_pushEQUIBRUN_"$i".tpr -maxwarn 1 
done
```
and run the individual simulations via e.g.
```
gmx mdrun -v -deffnm 3ptb_AMBER99SB_ben_pushEQUIBRUN_001
```
As these initial runs only require simulations of 0.1 ns length, they should be ready within a reasonably short time, i.e., some minutes.

When all equilibration simulations simulations have been carried out, prepare a separate directory for the pulling simulations and the individual pulling input TPR files via:
```
cd ..
mkdir v0.001
cd v0.001/

for i in {000..099}
do
gmx grompp -f ../3ptb_AMBER99SB_ben_pushRUN_v0.001.mdp -c ../equib/3ptb_AMBER99SB_ben_pushEQUIBRUN_"$i".gro -p ../3ptb_AMBER99SB_ben.top -n ../3ptb_AMBER99SB_ben.ndx -o 3ptb_AMBER99SB_ben_pushRUN_0.001_"$i".tpr
done
```
Note that the notation `*_0.001_*` stands for a velocity in Gromacs units of 0.001 nm/ps, i.e., 1 m/s. To our current experience, this is a sweet-spot velocity with the best trade-off between slow pulling and minimal computational effort. Run the simulations via e.g.:
```
gmx mdrun -v -deffnm 3ptb_AMBER99SB_ben_pushRUN_0.001_000
```
These simulations will each require 1-2 hours on a modern workstation, so you better run them in parallel on a HPC cluster of your choice.

For all further analysis, you require the `3ptb_AMBER99SB_ben_pushRUN_0.001_*_pullf.xvg` files (with `*` denoting the respective run number).


## dcTMD analysis

### Carrying out the analysis

Within the folder containing the `*pullf.xvg` pulling force files, use our dcTMD script as:
```
python3 NEQGamma.py -i 3ptb_AMBER99SB_ben_pushRUN_0.001_ -s _pullf -o 3ptb_AMBER99SB_ben_pushRUN_0.001_dG.dat -ofrict 3ptb_AMBER99SB_ben_pushRUN_0.001_frict.dat -vel 0.001 -T 290.15 -N 100 -av 40000 -sigma 40000
```
The respecitve flags read:
- `-i`: prefix of all force files (needs to be identic)
- `-s`: suffix of all force files excluding the `*.xvg` ending (needs to be identic, too)
- `-o`: output name for the file containing the dissipation-corrected free energy estimate with:
  - column #1: x axis in nm
  - column #2: non-equilibrium work <W>
  - column #3: fricition coefficient Gamma
  - column #4: dissipative work estimate 
  - column #5: free energy estimate 
- `-ofrict`: output name for the file containing a closer analysis of the friction:
  - column #1: x axis in nm
  - column #2: force autocorrelation function (with center time equal to the last time point)
  - column #3: friction coefficient Gamma
  - column #4: Gaussian filtered Gamma
  - column #5: running average window filtered Gamma
- `-vel`: pulling velocity used in the simulations in nm/ps
- `-T`: temperature used in the pulling simulations in K
- `-N`: number of input trajectories
- `-av`: width in data points of the running average window. We recommend to use a width of 40 to 200 per 1000 data points.
- `-sigma`: sigma width in data points of the Gaussian filter. We recommend to use a sigma of 40 per 1000 data points.

### Results of the analysis
  
A typical output of the free energy file for this tutorial would look like this:
![energy comparison](https://github.com/floWneffetS/tutorial_dcTMD/blob/main/figs/Tryp_100traj_energies.png)
  
The Gauss- and average window-filtered friction profiles (with a point width of 40000 for both `-av` and `-sigma`) would look like this:
![friction comparison](https://github.com/floWneffetS/tutorial_dcTMD/blob/main/figs/Tryp_100traj_friction.png)
  
Please note that the friction factors are very noisy and converge very badly. This is natural, as they technically are a measure of the force variance, whose estimator converges significantly slower than the one of the mean force. However, the only sources of this noise comes from thermal fluctuations and therefore can be removed via the two filtering functions. The exact width for the filter functions is a heuristic parameter that varies between investigated systems. As general guideline: remember that friction has to be always positive, so `-av` and `-sigma` should be chosen such that this constraint is fullfilled. We usually use the Gauss-filtered data for display and Langevin simulations.
  
### The friction overestimation artefact

In your analysis you may encounter a significat drop in free energy to unfeasibly low (and even negative) values up to several hundrets of kJ/mol. This clearly erroneous result stems from an overestimation of friction coming from the presence of different unbinding pathways. In this case you will need to cluster trajectories accordingly to the pathways they take out of the binding site and perform the friction correction for each cluster of trajectories separately. For more details, please consult [Wolf et al. Nat. Commun. 2020](https://www.nature.com/articles/s41467-020-16655-1) or the authors of this article directly.
  
  
## Temperature-boosted Langevin simulations

Take your self-generated free energy and friction files and generate two new files in which you 
- only contain the two columns with x and the free energy
- only contain the two columns with x and the Gauss-filtered friction.
You may want to reduce the resolution of both files down to 1 pm. 
  
Alternatively, you may download and unpack the respective TAR archive containing two such files via
```
tar -xzvf ./Tboost_tutorial_files.tar
```
and continue from there. Furthermore, download the C++ code, the CMake file and the two Jupyter notebooks from [the T-boosting repository](https://github.com/floWneffetS/Langevin_T_boost). Follow the instructions there to compile the Langevin integrator code.
  
You will need to decide on a suitable temperature range for your T-boosting calculations. This requires a bit of testing and is different for each system. In the case of Trypsin as used here, a temperature range of 400 to 900 K is a good start. As a rule of thumb, you should observe several 1000s of transitions over the main barrier for the highest temperture. The notebook `error_estimation_T_boosting.ipynb` helps you to estimate the number of temperature points and transitions for the convergence of transition rates. You need to generate the two files `start.dat` containing an arbitrary point along x to start the simulation (e.g., 1.0) and `Try_ligmass.dat` containing the ligand mass in kg/mol. Remember that in case of Trypsin as used here, transition rates are independent from the mass used, so we can here use an arbitrary high mass of 10 times benzamidine's mass and by this increase the Langevin time step to 5 fs (see again [Wolf et al. Nat. Commun. 2020](https://www.nature.com/articles/s41467-020-16655-1)). Start Langevin simulations using the Langevin integrator as:
```
for i in 400 450 500 550 600 650 700 750 800 850 900
do 
./LE_1dim_reflect -start start.dat -free Trypsin_all_1pm_dG.dat -gamma Trypsin_all_1pm_frict.dat -mass Try_ligmass.dat -o Trypsin_LE_5fsdt_5ms_1pmres_"$i"K.dat -t 0.000005 -T "$i" -I 89492"$i"90093 -L 1000000 -s 200000 -n 2001 -ngamma 2001 >& Trypsin_LE_5fsdt_5ms_1pmres_"$i"K.log & 
done
```
The program flags indicate:
- `-start`: name of file with starting point in nm
- `-free`: name of file with free energy (x,G(x)) in kJ/mol
- `-gamma`: name of file with friction (x,Gamma(x)) in units kg/(mol*ns), equals the kJ ps /(mol nm^2) from Gromacs
- `-mass`: name of file with mass in kg/mol
- `-o`: name of output trajectory with units of ns and nm
- `-t`: integration timestep in ns (so here: 5 fs)
- `-T`: temperature in K
- `-I`: seed for random number generator (we use the info from the bash varriable $i as a pseudo-random seed)
- `-s`: write out every sth point, so here a point in time is written out each nanosecond
- `-L`: length of output trajectory in points defined in `-s`, so here 1 ms
- `-n`: number of points of the free energy in the input free energy file
- `-ngamma`: number of points of Gamma in the input friction file
  
Here we gegnerate one single Langevin trajectory file of 1 ms length per temperature. On a single core of a modern PC, this should require ca. 5 hours. Alternatively, you may start 5 simulations with 0.2 ms length each ad different seed, which is actually beneficial for convergence.
  
After completion of the simulations, use the Pytho script `LE_bin_corer.py` to core the data (see e.g. [Nagel et al.](https://aip.scitation.org/doi/10.1063/1.5081767) for the background of this approach) and write out files conntaining the observed transition times. For trypsin, we decided to attribute distances <0.3 nm to the bound state coree, and distances >0.9 nm to the unnbound state core. Use the script simply as e.g.
```
./LE_bin_corer.py -i Trypsin_LE_5fsdt_5ms_1pmres_400K
```
ommiting the `.dat` ending. The output will be the two files `Trypsin_LE_5fsdt_5ms_1pmres_400K_left_to_right.dat` and `Trypsin_LE_5fsdt_5ms_1pmres_400K_right_to_left.dat`. `left` signifies here the "left" side of the x axis, i.e., the bound state, and `right` the "right" side, i.e., the unbound state. Hence, `*left_to_right.dat` contains statistics on the unbinding waiting times, and `*right_to_left.dat` on the binding waiting times. **Note:** the output "times" are strictly in the time scale of steps written out in `Trypsin_LE_5fsdt_5ms_1pmres_"$i"K.dat`, so if this file is in nanoseconds as in the example here, your waiting times will be in nanoseconds, as well. Calculate mean waiting times for the respective temperatures and processes from these files. Use then a fit of the mean waiting times and their errors to the T-boosting formula (4) in [Wolf et al. Nat. Commun. 2020](https://www.nature.com/articles/s41467-020-16655-1) to recover an estimate of kinetics at 290.15 K at which the initial MD pulling simulations were carried out. Remember that the binding rate calculated by you is a function of the ligand's concentration, which is 1 ligand per spherical volume based on the maximal pulling distance as radius. Alternatively, you can use our notebook `LE_Tboost_example.ipynb`available [here](https://github.com/floWneffetS/Langevin_T_boost) for this analysis. Here you can directly use the two types of output files from `LE_bin_corer.py` as input. 
