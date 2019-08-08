%{
---------------------------------------------------------------------------
    Name: Shift hologram of a makred shaped on two layered galsses 
    
    Author: Peng Xia
    Date: 2019
---------------------------------------------------------------------------
%}
holo_type = 'complex';  


lambda = 532e-9;         % the laser wavelength
pps = 3.45e-6;    % CCD pixel size; 3272x2469 

z = 10e-3+[0 9.5]*1e-3;   % the distance between the object and the CCD camera

Nz = length(z);

tau = 0.2;   % This effects, need further investigation
tau_psi = 0.2;
%% ========================================== Resize =================================================
load('mark.mat');

[holo, pps] = holoResize(holo, pps, 512, 5);
[Ny, Nx] = size(holo);

