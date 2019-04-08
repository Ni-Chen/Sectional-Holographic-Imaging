% data size
Nx = 256;
Ny = 256;
Nz = 1;

%% System specification
lambda  = 532e-9;  % wavelength 
pps     = 5e-6;    % pixel pitch of the hologram 
z0      = 100e-3;  % 40, distance from detector to center of the object 

Lx = Nx*pps;  % size of detector (um)
NA = Lx/2/sqrt(z0^2+(Lx/2)^2)

% Resolution of diffraction limits
res_z = 2*lambda/(NA^2)
res_x = 0.61*lambda/NA

rp = round(res_x/pps/2)   % radius of a point

% Shift 
x_shift = round(res_x/pps/1.5);

dz = res_z;
%% Indirect parameters
% Center of the object is located at the origin, z0 is the location of the hologram
z = z0 - ((1:Nz)- round(Nz/2))*dz;
% z_axis = (z - z0)*1e3;

x_axis = ((1:Nx)- round(Nx/2))*pps;
