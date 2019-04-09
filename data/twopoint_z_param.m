% Sampling of the two points setting
Nx = 256;
Ny = 1;
Nz = 31;

% location of the two points
z_shift = 1;
x_shift = 0;

%% specification of the holographic System 
lambda  = 532e-9;  % wavelength 
pps     = 8e-6;   % pixel pitch of the hologram 
z0      = 80e-3;   % distance from detector to center of the object 

Lx = Nx*pps;  % size of detector (um)
% D = sqrt(2)*Lx/2;
D = 1*Lx/2;
NA = Lx/2/sqrt(z0^2 + D^2)
% NA = Lx/2/sqrt(z0^2)

%% Rayleigh’s diffraction resolution
% <https://micro.magnet.fsu.edu/primer/digitalimaging/deconvolution/deconresolution.html>
%
% <https://www.microscopyu.com/techniques/super-resolution/the-diffraction-barrier-in-optical-microscopy>
%
% $\delta z = 2\frac{\lambda}{NA^2}$, 
% $\delta x = 0.61\frac{\lambda}{NA}$
% 
res_z = 2*lambda/(NA^2)  
res_x = 0.61*lambda/NA

dz_factor = 2;
dz = res_z/(z_shift*2)/dz_factor;  % depth interval of the object, make the resolution as the reference

%% Indirect parameters
% Center of the object is located at the origin, z0 is the location of the hologram
z = z0 - ((1:Nz)- round(Nz/2))*dz
z_axis = (z-z0)*1e3;


