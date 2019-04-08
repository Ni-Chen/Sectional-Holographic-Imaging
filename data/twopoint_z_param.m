% data size
Nx = 256;
Ny = 1;
Nz = 31;

% Shift 
z_shift = 1;
x_shift = 0;

%% System specification
lambda  = 532e-9;  % wavelength 
pps     = 8e-6;   % pixel pitch of the hologram 
z0      = 80e-3;   % distance from detector to center of the object 

Lx = Nx*pps;  % size of detector (um)
D = sqrt(2)*Lx/2;
% D = 1*Lx/2;
NA = Lx/2/sqrt(z0^2 + D^2)
% NA = Lx/2/sqrt(z0^2)

% Abbe’s diffraction formula for axial resolution
% res_z = 2*lambda/(NA^2)  
% res_x = lambda/2/NA

% Rayleigh’s diffraction formula for axial resolution
% https://www.microscopyu.com/techniques/super-resolution/the-diffraction-barrier-in-optical-microscopy
% https://micro.magnet.fsu.edu/primer/digitalimaging/deconvolution/deconresolution.html
res_z = 2*lambda/(NA^2)  
res_x = 0.61*lambda/NA

% res_x = z0*lambda/Lx
% res_z = 8*lambda*z0^2/Lx

dz_factor = 4;
dz = res_z/(z_shift*2)/dz_factor;  % depth interval of the object, make the resolution as the reference

%% Indirect parameters
% Center of the object is located at the origin, z0 is the location of the hologram
z = z0 - ((1:Nz)- round(Nz/2))*dz
% z = (0:Nz-1)*dz
z_axis = (z-z0)*1e3;




