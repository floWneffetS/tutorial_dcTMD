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
These simulations will require 1-2 hours on a modern workstation, so you better run them in parallel on a HPC cluster of your choice.

For all further analysis, you require the `3ptb_AMBER99SB_ben_pushRUN_0.001_*_pullf.xvg` files (with `*` denoting the respective run number).


## dcTMD analysis

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
- `-T`: temperature used in the pulling simulations
- `-N`: number of input trajectories
- `-av`: width in data points of the running average window. We recommend to use a width of 40 to 200 per 1000 data points.
- `-sigma`: sigma width in data points of the Gaussian filter. We recommend to use a sigma of 40 per 1000 data points.
  
A typical output of the free enrgy file would look like this:
![energy comparison](https://github.com/floWneffetS/tutorial_dcTMD/blob/main/figs/Tryp_100traj_energies.png)
  
Please note that the friction factors are very noisy and converge very badly. This is natural, as they technically are a measure of the force variance, whose estimator converges significantly slower than the one of the mean force. However, the only sources of this noise comes from thermal fluctuations and therefore can be removed via the two filtering functions. The exact width for the filter functions is a heuristic parameter that varies between investigated systems. As general guideline: remember that friction has to be always positive, so `-av` and `-sigma` should be chosen such that this constraint is fullfilled. 
  
  
## Temperature-boosted Langevin simulations
