% data size
Nx = 64;
Ny = 64;
Nz = 5;

%% System specification
lambda  = 633e-9;  % wavelength 
pps     = 30e-6;   % pixel pitch of the hologram 
z0      = 0e-3;   % distance from detector to center of the object 

Lx = Nx*pps;  % size of detector (um)
NA = Lx/2/sqrt(z0^2+(Lx/2)^2)

% Resolution of diffraction limits
res_z = 2*lambda/(NA^2)
% res_x = lambda/2/NA

% dz = res_z*1
dz = 10e-3;  % axial spacing (um)


%% Indirect parameters
% Center of the object is located at the origin of the coordinates, z0 is the location of the
% hologram
% z = z0 - ((1:Nz)- round(Nz/2))*dz;
z = z0 + [0:(Nz-1)]*dz;