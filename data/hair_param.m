%{
---------------------------------------------------------------------------
    Name: Experiment of on-axis transmittant digital holography.
    
    Author: Gao  
    Date: Nov. 2018    
    
    Object: hair on glass
---------------------------------------------------------------------------
%}

holo_type = 'complex';  

lambda = 532e-9;    % wavelength of illumination beam
pps = 3.45e-6;      % pixel pitch of CCD camera

% z = [132 -4 3]*1e-3;
% z = [0 8 135]*1e-3;
z = [-16 -8 0 8 16 24 32 40 48 56 64 72 80 88 96 104 120 128 135 143]*1e-3;
% z = [85 -50 -42]*1e-3;
% z = [132 133 134 135 136]*1e-3;
Nz = length(z);

% offsetZ = 68e-3;        % Distance between object and camera
% Nz = 10;
% z = zd - ((1:Nz)- round(Nz/2)).*dz; 

tau = 0.2;   % This effects, need further investigation
tau_psi = 0.2;
%% ========================================== Resize =================================================
load('hair.mat');

[holo, pps] = holoResize(holo, pps, 512, 5);
[Ny, Nx] = size(holo);


