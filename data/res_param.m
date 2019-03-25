%{
---------------------------------------------------------------------------
    Name: Experiment of on-axis transmittant digital holography.
    
    Author: Gao  
    Date: Nov. 2018    
    
    Object: hair on glass
---------------------------------------------------------------------------
%}

load('res.mat');
% data size
[Ny, Nx] = size(holo);

lambda = 632.8e-9;    % wavelength of illumination beam
pps = 3.45e-6;      % pixel pitch of CCD camera
offsetZ = 123e-3;        % Distance between object and camera

deltaX = pps;
deltaY = pps;


deltaZ = 10e-3;
Nz = 10;

% % size of detector (um)
% sensor_size = Nx*pps;


% z_scope = linspace(offsetZ-5e-3, offsetZ+5e-3, Nz);
z_scope = offsetZ - ((1:Nz)- round(Nz/2))*deltaZ

%% ========================================== Crop =================================================
N_crop = Ny;

holo_crop = holo(Ny/2-N_crop/2+1:Ny/2+N_crop/2, Nx/2-N_crop/2+1:Nx/2+N_crop/2);
clear holo;
holo = holo_crop;

[Ny, Nx] = size(holo);

% size of detector (um)
sensor_size = N_crop*pps;

