close all; clear; clc;

addpath(genpath('.../function/'));  % Add funtion path with sub-folders

lambda = 532e-9;    % wavelength of illumination beam
pps = 3.45e-6;      % pixel pitch of CCD camera

load('Amp_0.003m_dz.mat');
holo_3_amp = double(Amp);
holo_3_phase = importdata('Phase_0.003m.txt');

holo_3 = holo_3_amp.*exp(1i*holo_3_phase);

[Ny, Nx] = size(holo_3);

%% ======================================= Hologram ================================================
[otf3d, ~, pupil3d] = OTF3D_z(Ny, Nx, lambda, pps, -3e-3);

%% ========================== Reconstruction with back-propagation =================================
holo_capture = iMatProp3D(holo_3, otf3d, pupil3d, 'complex');

holo = holo_capture;
holo = max(holo(:)) - holo; % holo = holo*(-1)+abs(min(min(holo*(-1))));

holo = holoNorm(holo);
holo = imadjust(abs(holo)).*exp(1i*angle(holo));
holo = holoNorm(holo);

save(['../hair.mat'], 'holo');
 
[otf3d, psf3d, pupil3d] = OTF3D_z(Ny, Nx, lambda, pps, [132e-3 3e-3 10e-3]);
holo_reobj = iMatProp3D(holo, otf3d, pupil3d, 'complex');

% Back vecProp reconstruction
figure; imagesc(plotdatacube(abs(holo))); title('Hologram'); axis image; drawnow; colormap(hot); colorbar; axis off;
figure; imagesc(plotdatacube(abs(holo_reobj))); title('BackPropagation'); axis image; drawnow; colormap(hot); colorbar; axis off;



