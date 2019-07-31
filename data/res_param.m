%{
---------------------------------------------------------------------------
    Name: Experiment of on-axis transmittant digital holography.
    
    Author: Gao  
    Date: Nov. 2018    
    
    Object: hair on glass
---------------------------------------------------------------------------
%}

% load('res.mat');
holo = double(imread('resthu.png'));

holo_type = 'inline';  

% data size
[Ny, Nx] = size(holo);
Nz = 1;

lambda = 532e-9;    % wavelength of illumination beam
pps = 3.8e-6;      % pixel pitch of CCD camera
z0 = 14e-3;        % Distance between object and camera

Lx = Nx*pps;  % size of detector (um)
NA = Lx/2/sqrt(z0^2+(Lx/2)^2)

% Resolution of diffraction limits
res_z = 2*lambda/(NA^2)
% res_x = lambda/2/NA

% dz = res_z*1
dz = 14e-3;  % axial spacing (um)

%% Indirect parameters
% Center of the object is located at the origin of the coordinates, z0 is the location of the
% hologram
% z = z0 - ((1:Nz)- round(Nz/2))*dz;
z = z0 + [0:(Nz-1)]*dz

holo = max(holo(:)) - holo; % holo = holo*(-1)+abs(min(min(holo*(-1))));

holo = holo./max(abs(holo(:)));

[holo, pps] = holoResize(holo, pps, 512, 0);
[Ny, Nx] = size(holo);
