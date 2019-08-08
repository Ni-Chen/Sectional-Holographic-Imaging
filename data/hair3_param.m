%{
---------------------------------------------------------------------------
    Name: Shift hologram of a three overlapped hair 
    
    Author: Peng Xia
    Date: 2019 
    
---------------------------------------------------------------------------
%}
holo_type = 'phaseshift';

lambda = 532e-9;         % the laser wavelength
pps = 3.45e-6;    % CCD pixel size; 3272x2469 

z = [2 5 14]*1e-3;   % the distance between the object and the CCD camera


load('hair3.mat');

holo = holo - max(abs(holo(:)));
[holo, pps] = holoResize(holo, pps, 512, 5);
[Ny, Nx] = size(holo);