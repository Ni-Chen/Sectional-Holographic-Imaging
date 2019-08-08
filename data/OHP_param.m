%{
---------------------------------------------------------------------------
    Name:  off-axis transmittant digital hologram of OHP
    
    Author: Peng Xia 
    Date: Nov. 2018    
    
    Object: 
---------------------------------------------------------------------------
%}

holo_type = 'complex';  

lambda = 473e-9;    % wavelength of illumination beam
pps = 3.45e-6;      % pixel pitch of CCD camera

% z = [132 -4 3]*1e-3;
% z = [136 3 10]*1e-3;
z = [173 179 192]*1e-3;
% z = [132 133 134 135 136]*1e-3;
Nz = length(z);

% offsetZ = 68e-3;        % Distance between object and camera
% Nz = 10;
% z = zd - ((1:Nz)- round(Nz/2)).*dz; 

tau = 0.2;   % This effects, need further investigation
tau_psi = 0.2;
%% ========================================== Resize =================================================
load('OHP.mat');

[holo, pps] = holoResize(holo, pps, 512, 0);
[Ny, Nx] = size(holo);


