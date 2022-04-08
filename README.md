# four-flux-method

based on 
Four-flux models to solve the scattering transfer equation in terms of Lorenz-Mie parameters
by
B. Maheu, J. N. Letoulouzan, and G. Gouesbet
https://doi.org/10.1364/AO.23.003353

the code works from validation.m file

it uses C Mätzler Lorenz Mie code to solve particle characteristics

Currently it uses Find_R_T_size_dist to take into accound particle size distribution
for monodisperse particles Find_R_T can be used

in the validation.m file;
lamda is wavelength spectrum in meters
f_v is particle volume fraction
h is coating thickness in meters


in the Find_R_T_size_dist.m file;
r_avg is averge particle radius in meter
