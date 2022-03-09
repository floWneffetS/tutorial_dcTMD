# About
Tutorial for using the usage of the dcTMD approach and Langevin simulations using the trypsin-benzamidine complex (see [Wolf et al., Nat. Commun. 2020, 11, 2918](https://www.nature.com/articles/s41467-020-16655-1)) as example. Required programs are Gromacs v2018 and higher, the script NEQGamma.py from the dcTMD repository, and LE_1dim_reflect.cpp from the Langevin_T_boost repository.

The basis of this tutorial are constraint pulling simulations as described in [Wolf and Stock, JCTC 2018, 14, 6175](https://pubs.acs.org/doi/full/10.1021/acs.jctc.8b00835). For pathway separation (advanced topic, see below), we recommend the usage of a principal component analysis (PCA) as possible via our [fastPCA code](https://github.com/moldyn/FastPCA).

# Licensing

If you perform a dcTMD analysis of your own data for published works, please cite the appropriate literature:

- for the dcTMD analysis itself: _Wolf, S., & Stock, G. (2018). Targeted Molecular Dynamics Calculations of Free Energy Profiles Using a Nonequilibrium Friction Correction. Journal of Chemical Theory and Computation, 14, 6175–6182._
- for temperature-boosted Langevin equation simulations and/or in case you need to perform a pathway separation: _Wolf, S., Lickert, B., Bray, S., & Stock, G. (2020). Multisecond ligand dissociation dynamics from atomistic simulations. Nature Communications, 11(1), 2918. http://doi.org/10.1038/s41467-020-16655-1_

# Tutorial

## TMD simulation

### Input files
Download and unpack tutorial_files.tar.gz. You will find a folder with the following files:

Files for run input commands:
`3ptb_AMBER99SB_ben_pushEQUIBRUN.mdp 
3ptb_AMBER99SB_ben_pushRUN_v0.001.mdp
`
The `pushEQUIBRUN.mdp` file is an initial equilibration file for generating start simulation files with different initial velocity distributions. The `pushRUN_v0.001.mdp` file is the respective command input for the non-equilibrium pulling. `3ptb` refers to the protein data base code of a trypsin structure, `AMBER99SB` to the employed force field, and `ben` to the benzamidine ligand.

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


### Carrying out simulations



## dcTMD analysis

## Langevin simulations
