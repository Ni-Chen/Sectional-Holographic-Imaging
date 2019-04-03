% data size
Nx = 64;
Ny = 64;
Nz = 60;

%% System specification
lambda  = 632.8e-3;  % wavelength (um)
pps     = 5;         % pixel pitch of the hologram (um)
offsetZ = 8e4;       % distance from detector to center of the object (um)

deltaZ  = 500;       % distance between each axial plane (um)


%% Indirect parameters
sensor_size = Nx*pps;  % size of detector (um)

deltaX = pps;
deltaY = pps;

z = offsetZ - ((1:Nz)- round(Nz/2))*deltaZ;

% z_depth = Nz*deltaZ