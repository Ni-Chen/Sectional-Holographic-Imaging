%{
OSH of 'S' and 'H'.

Parameters can be found from:
- Xin Zhang, "Blind sectional image reconstruction for optical scanning holography," Opt. Lett. 34, 
  3098-3100 (2009)

z1=87, z2=107cm.
in HKU, relative k of k001 = 555 = 89.452cm; k002 = 776 = 63.977cm; since rk=pi/lambda/z, these don't 
meet the original z, different deltax may be used in the reconstruction
%}

holo_type = 'complex';  

load('sh.mat');
[Ny, Nx] = size(holo);

% wavelength (um), 632nm
lambda = 632.8e-9;  % Need to confirm

sensor_size = 10e-3; % size of detector (um), 1cm used in HKU reconstruction
pps = sensor_size/Nx;   

% deltaZ = 8e-2;  % distance between each axial plane (um)
% offsetZ = 77e-2;
% z = offsetZ - ((1:Nz)- round(Nz/2))*deltaZ 

% z = [65 89]*1e-2;
z = [61 65 69 73 77 81 85 89 93 97]*1e-2;  % Even number
Nz = length(z);

tau = 0.01;   % This effects, need further investigation
tau_psi = 0.15;

%% Resize 
load('sh.mat');
[holo, pps] = holoResize(holo, pps, 256, 10);
[Ny, Nx] = size(holo);