%{
OSH of 'S' and 'H'.

Parameters can be found from:
- Xin Zhang, "Blind sectional image reconstruction for optical scanning holography," Opt. Lett. 34, 
  3098-3100 (2009)

z1=87, z2=107cm.
in HKU, relative k of k001 = 555 = 89.452cm; k002 = 776 = 63.977cm; since rk=pi/lambda/z, these don't 
meet the original z, different deltax may be used in the reconstruction
%}

% data size
Nx = 500;
Ny = 500;
Nz = 8;  % 85~110

% wavelength (um), 632nm
lambda = 632.8e-9;  % Need to confirm

sensor_size = 10e-3; % size of detector (um), 1cm used in HKU reconstruction

pps = sensor_size/Nx;   
deltaX = pps;
deltaY = pps;

deltaZ = 8e-2;  % distance between each axial plane (um)

offsetZ = 77e-2;
z_scope = offsetZ - ((1:Nz)- round(Nz/2))*deltaZ 


%% Pad
N = 512;
pad_size = (N-Nx)/2;
holo_pad = padarray(holo, [pad_size pad_size]);  

[Nx, Ny] = size(holo_pad);  
holo = holo_pad;

sensor_size = Nx*pps;