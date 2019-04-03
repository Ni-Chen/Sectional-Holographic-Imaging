% data size
Nx = 64;
Ny = 64;
Nz = 60;

%% System specification
lambda  = 632.8e-9;  % wavelength (um)
pps     = 5e-6;         % pixel pitch of the hologram (um)
offsetZ = 80e-3;       % distance from detector to center of the object (um)

deltaZ  = 0.5e-3;       % distance between each axial plane (um)

%% Indirect parameters
sensor_size = Nx*pps;  % size of detector (um)

deltaX = pps;
deltaY = pps;

% Center of the object is located at the origin of the coordinates, offsetZ is the location of the
% hologram
z = offsetZ - ((1:Nz)- round(Nz/2))*deltaZ