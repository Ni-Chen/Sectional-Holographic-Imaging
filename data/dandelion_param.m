%{ 
Gabor hologram made at  Duke University.

z1 = 1.5cm, z2=5.5cm

%}

load('dandelion.mat');

holo_type = 'inline';  

% number of pixels of the detected images
pixel_num = 1024;

% size of detector pixels (um)
pps = 5.2e-6;

lambda = 632.8e-9;  % wavelength (um)

deltaZ = 5e-3;  % distance between each axial plane (um)
offsetZ = 35e-3;  % distance from detector to first reconstructed plane (um)

Nz = 12;   % number of axial planes, should cover the depth range

z = offsetZ -((1:Nz)- round(Nz/2))*deltaZ 

sigma = 0.0277;
%% ==================================================================
pad_size = 0;   % 100, % number of zeros to pad matrix by in each direction

shrinkage_factor = pixel_num/size(holo, 1)
sensor_size = pixel_num*pps; 
deltaX = pps*shrinkage_factor; 
deltaY = pps*shrinkage_factor; 

holo_pad = padarray(holo, [pad_size pad_size]);  
range = pad_size*2 + pixel_num; 

[Nx, Ny] = size(holo_pad);  
holo = holo_pad;
