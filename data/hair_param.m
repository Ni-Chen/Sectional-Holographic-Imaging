%{
---------------------------------------------------------------------------
    Name: Experiment of on-axis transmittant digital holography.
    
    Author: Gao  
    Date: Nov. 2018    
    
    Object: hair on glass
---------------------------------------------------------------------------
%}

load('hair.mat');
% data size
[Ny, Nx] = size(holo);

lambda = 632.8e-9;    % wavelength of illumination beam
pps = 3.45e-6;      % pixel pitch of CCD camera
offsetZ = 180e-3;        % Distance between object and camera

deltaX = pps;
deltaY = pps;


deltaZ = 10e-3;
Nz = 11;

% size of detector (um)
sensor_size = Nx*pps;


z_scope = linspace(offsetZ-50e-3, offsetZ+50e-3, Nz);

%% ========================================== Crop =================================================
N_crop = Ny;

% size of detector (um)
sensor_size = N_crop*pps;

holo_crop = holo(Ny/2-N_crop/2+1:Ny/2+N_crop/2, Nx/2-N_crop/2+1:Nx/2+N_crop/2);
clear holo;
holo = holo_crop;

[Ny, Nx] = size(holo);

