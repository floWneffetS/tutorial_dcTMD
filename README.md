# tutorial_dcTMD
Tutorial for using the usage of the dcTMD approach and Langevin simulations using the trypsin-benzamidine complex as example. Required programs are Gromacs v2018 and higher, the script NEQGamma.py from the dcTMD repository, and LE_1dim_reflect.cpp from the Langevin_T_boost repository.

The basis of this tutorial are constraint pulling simulations as described in 

## Usage
Download and unpack tutorial_files.tar.gz. You will find a folder with the following files:

Files for run input commands:
`3ptb_AMBER99SB_ben_pushEQUIBRUN.mdp 
3ptb_AMBER99SB_ben_pushRUN_v0.001.mdp
`

Topologies:
`3ptb_AMBER99SB_ben.top
3ptb_AMBER99SB_ben_Protein_chain_A.itp
3ptb_AMBER99SB_ben_Ion_chain_B.itp
3ptb_ben_H2_GMX_RESP.itp
posre_Protein_chain_A.itp
posre_Ion_chain_B.itp
posre_ben.itp
`
