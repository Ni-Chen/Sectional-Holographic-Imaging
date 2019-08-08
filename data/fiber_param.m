%{
---------------------------------------------------------------------------
    Name: Experiment of off-axis transmittant digital holography.
    
    Author: Xia 
    Date: Nov. 2018    
    
    Object: 
---------------------------------------------------------------------------
%}

holo_type = 'complex';  

lambda = 473e-9;    % wavelength of illumination beam
pps = 3.45e-6;      % pixel pitch of CCD camera

% z = [138 155]*1e-3;
z = [138 141 155]*1e-3;   % can not smaller than 3
% z = [135 138 141 155 159]*1e-3;   % can not smaller than 3

Nz = length(z);

tau = 0.2;   % This effects, need further investigation
tau_psi = 0.2;
%% ========================================== Resize =================================================
load('fiber.mat');

[holo, pps] = holoResize(holo, pps, 512, 5);
[Ny, Nx] = size(holo);


