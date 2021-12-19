# Spin_fluctuations

Codes for comm_phys paper:

The following are the codes necessary to replicate the results in "Coupled charge and spin excitations in the cuprates and their role in the superconducting transition". Since the calculations were done for different data sets (at different temperatures and dopings), only the codes for the  OP80 sample at 70 K are given. This set was chosen because it represents the more interesting case where the eigenvalue is close to unity. The codes were written in MATLAB. To run them, the OP80 ARPES data (at 70K) must be used as input. 

Here is the suggested order (it is recommended that the codes be run in succession, and that the output variables are kept in the workspace):

import_data: imports data and relevant scales (energy, angles, etc.) from the raw data files as they come from the electron detector software. The input data come in the form of 16-19 "slices", corresponding to dispersions along the ky direction for 16-19 kx values.   

BZ_quadrant: takes the ARPES data and organizes them into a 3D array, for ease of manipulations. This script also subtracts background, superlattice, and "symmetrizes" the data. Requires ARPES data, stored in variable called "data", plus other variables created in the code "import_data". The main output variable is the spectral function sampled over one BZ quadrant, and covering energies in the range [-640,640] meV. This output is called "bsces70". 

realw_chi0
Calculates the bare spin susceptibility starting from the variable "bsces70". 
The output is chi_0 (called "chi0_70"), sampled over the entire BZ and in the energy range [0,400] meV.

chi_w_wm:
Calculates the real and Matsubara spin response functions. Requires "chi0_70" as input.
The outputs are 
i) Im(chi) throughout the BZ and in the energy range [-400,400]  meV (called "chi") and
ii) The Matsubara chi throughout the BZ in the energy range [-3.5,3.5] eV, called ("chiM")
***Note: in this code, one must also choose the appropriate spin-fermion coupling "U" to reproduce Im(chi) seen in INS experiments. For this data set, U=670 meV.***

matsubara_G_BSE"
Calculates the Matsubara Green's function and solves the BSE via the POWER METHOD.
Requires "bsces70", "chiM", and "U" as inputs. The output is the leading eigenvalue and eigenfunction in Matsubara frequencies. 
