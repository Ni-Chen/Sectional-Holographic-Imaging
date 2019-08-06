close all; clear; clc;

addpath(genpath('../function/'));  % Add funtion path with sub-folders

lambda = 532e-9;    % wavelength of illumination beam
pps = 3.45e-6;      % pixel pitch of CCD camera

holo_amp = importdata('Amp_0.013m.txt');
holo_phase = importdata('Phase_1e-05m.txt');

holo = holo_amp.*exp(1i*holo_phase);

[Ny, Nx] = size(holo);

holo = max(holo(:)) - holo; 
% holo = holo*(-1)+abs(min(min(holo*(-1))));

% holo = (holo)./(max(abs(holo(:))));
holo = holoNorm(holo);
holo = imadjust(abs(holo)).*exp(1i*angle(holo));
holo = holoNorm(holo);
holo = (holo)./(max(abs(holo(:))));

save(['../beads.mat'], 'holo');

z = linspace(-20e-3, 20e-3, 41);
[otf3d, psf3d, pupil3d] = OTF3D(Ny, Nx, lambda, pps, z);
holo_reobj = iMatProp3D(holo, otf3d, pupil3d);
write3d(abs(holo_reobj), z*1e3, './', 'beads');

