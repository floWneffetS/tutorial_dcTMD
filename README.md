# tutorial_dcTMD
Tutorial for using the usage of the dcTMD approach and Langevin simulations using the trypsin-benzamidine complex (see Wolf et al., Nat. Commun. 2020, 11, 2918) as example. Required programs are Gromacs v2018 and higher, the script NEQGamma.py from the dcTMD repository, and LE_1dim_reflect.cpp from the Langevin_T_boost repository.

The basis of this tutorial are constraint pulling simulations as described in 

## TMD simulation

### Input files
Download and unpack tutorial_files.tar.gz. You will find a folder with the following files:

Files for run input commands:
`3ptb_AMBER99SB_ben_pushEQUIBRUN.mdp 
3ptb_AMBER99SB_ben_pushRUN_v0.001.mdp
`
The `pushEQUIBRUN.mdp` file is an initial equilibration file for generating start simulation files with different initial velocity distributions. The `pushRUN_v0.001.mdp` file is the respective command input for the non-equilibrium pulling. `3ptb` refers to the protein data base code of a trypsin structure, `AMBER99SB` to the employed force field, and `ben` to the benzamidine ligand.

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
Important: the index file needs to include an anchor group from where the ligand is pulled away (in this case: the group `[sheet]` containing C-alpha atoms from the central beta-sheet) and the ligand itself (or better, the heavy atoms of the ligand, here group `[ BEN_heavy ]`).

Structure file: `3ptb_AMBER99SB_ben.gro` in which the trypsin-bezamidine complex is equilibrated in TIP3P water with a physiological NaCl concentration and a single Ca2+ ion.

### Carrying out simulations

## dcTMD analysis

## Langevin simulations
