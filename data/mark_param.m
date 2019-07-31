%{
---------------------------------------------------------------------------
    Name: Experiment of phase shift hologram.
    
    Author: Min-Ok Jeong
    Date: 2008  
    
    Object: coin
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

[holo, pps] = holoResize(holo, pps, 256, 0);
[Ny, Nx] = size(holo);

