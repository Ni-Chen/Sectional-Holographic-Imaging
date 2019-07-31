%{
OSH of fluorescent beads(DukeR0200, 2um in diameter, excitation around 542nm, emission around 612nm).

The beads mainly assemble at the top and bottom surfaces. The distance between the two planes was 
around 35um. 
  - z0 = 140um
  - z1 = z0-15um, z2= z0+20um

The plane of the Fabry-Perot is imaged in the pupil plane of the objective with a magnification of 
20x, and the Fabry-Perot fringes localized in the focal plane of lens L1 are imaged in the focal 
plane of the objective with a magnification of 1/60x. Therfore the total manification is:
  - magnification = 20*(1/60)
  - FOV = 140um x 140um

- Guy Indebetouw, "Scanning holographic microscopy with spatially incoherent sources: reconciling the 
holographic advantage with the sectioning advantage," J. Opt. Soc. Am. A 26, 252-258 (2009)
- Edmund Y. Lam, "Three-dimensional microscopy and sectional image reconstruction using optical 
scanning holography," Appl. Opt. 48, H113-H119 (2009)
  - z0 = 140um
  - z1 = 85um, z2= 120um

- expected axial resolution is thus ~10um
%}

load('beads.mat');

holo_type = 'complex';  

% data size
[Ny, Nx] = size(holo);

lambda = 532e-9;  % wavelength (um)

% sensor_size = 140e-6/(20*(1/60));  % size of detector (um)
sensor_size = 105e-6/(20*(1/60));    % size of detector (um), FOV/M
pps = sensor_size/Nx;

% Nz = 14;
% deltaZ = 5e-6;  % distance between each axial plane (um)
% offsetZ = 105e-6;  % distance from detector to center of the object plane (um)
% z_scope = offsetZ - ((1:Nz)- round(Nz/2))*deltaZ

z = [85 125]*1e-6;
% z = [84 85 86 124 125 126 127 128]*1e-6;
Nz = length(z);

tau = 0.1;   % This effects, need further investigation
tau_psi = 0.2;

%% Resize 
[holo, pps] = holoResize(holo, pps, 512, 0);
[Ny, Nx] = size(holo);

