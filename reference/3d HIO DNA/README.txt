%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D RECONSTRUCTION BY ITERATIVE HYBRID INPUT OUTPUT ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Citation for this code/algorithm or any of its parts:
% Tatiana Latychevskaia and Hans-Werner Fink
% "Three-dimensional double helical DNA structure directly revealed 
% from its X-ray fiber diffraction pattern by iterative phase retrieval",
% Optics Express 26(23), 30991 - 31017 (2018) 
% https://doi.org/10.1364/OE.26.030991
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Tatiana Latychevskaia, 2015
% The version of Matlab for this code is R2010b

1. Program "a_simulate_3d_dp_no_phase.m" simulates 3d angular-averaged 
   diffraction pattern of discontinuous double helix. The simulated
   diffraction pattern is stored as 2D slices in BIN format.

2. To reconstruct the 3D diffraction patter, start "b_HIO_3d_sequence.m" and
   select the 2D slices simulated in step (1). Obtain several reconstructions, 
   300 iterations for each reconstruction is sufficient. The obtained 3D 
   reconstructions are saved as sequence of 2D slices, the name contains the
   calculated error.

3. The final reconstruction can be obtained either as the reconstruction with 
   the smallest error or as an average over several reconstructions with the 
   smallest error.  

4. To view 3D reconstruction, run "c_reconstruction_show_3d.m"