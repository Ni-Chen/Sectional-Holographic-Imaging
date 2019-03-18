%{
Gabor hologram of particles.

partical diameter = 10um

- magnification = 10x
- FOV = 
%}

% z = 72e3;          % um
pps = 3.45;          % um
lambda = 0.6328;     % wavelength in um

holo = double(imread('particle10x2.bmp'));


[Ny_ori, Nx_ori] = size(holo);
Nz = 3;
deltaX = pps;
deltaY = pps;

% distance between each axial plane (um)
deltaZ = 1000;
offsetZ = 72000;

z_scope = ((1:Nz)- round(Nz/2))*deltaZ + offsetZ

%% Crop hologram

N_crop = 1000;
Nx = N_crop;
Ny = N_crop;

% size of detector (um)
sensor_size = Nx*pps;

holo_crop = holo(Ny_ori/2-N_crop/2:Ny_ori/2+N_crop/2-1, Nx_ori/2-N_crop/2:Nx_ori/2+N_crop/2-1);
clear holo;
holo = holo_crop;