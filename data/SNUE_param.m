% data size
Nx = 64;
Ny = 64;
Nz = 5;

%% System specification
lambda  = 633e-9;  % wavelength 
pps     = 30e-6;   % pixel pitch of the hologram 
offsetZ = 120e-3;   % distance from detector to center of the object 

deltaZ  = 15e-3;   % distance between each axial plane

%% Indirect parameters
sensor_size = Nx*pps;  % size of detector (um)

deltaX = pps;
deltaY = pps;

% Center of the object is located at the origin of the coordinates, offsetZ is the location of the
% hologram
z = offsetZ - ((1:Nz)- round(Nz/2))*deltaZ;
% z = offsetZ + ((0:Nz-1))*deltaZ;

