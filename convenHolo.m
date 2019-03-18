%{
---------------------------------------------------------------------------
Name: Hologram of 3D objects.

Author:   Ni Chen (chenni@snu.ac.kr)
Date:     
Modified:
---------------------------------------------------------------------------
%}

close all;
clear;
clc;

format short;

addpath('./function');

indir = './data/';
outdir = './output/';

% random_scatter, geo_overlap, helix_circular, helix_conical
obj_name = 'helix_conical'; 
run([indir, obj_name, '_param.m']);
load([indir, obj_name, '.mat']);

%% ========================================= Propagation ===========================================
[otf3d, psf3d, pupil3d] = OTF3D(Nx, Ny, Nz, lambda, deltaX, deltaY, deltaZ, offsetZ, sensor_size);
prop_field = Propagation3D(obj3d, otf3d, pupil3d);  % Field at the center plane of the 3D object
holo = prop_field;
% holo = prop_field.^2;  % Gabor hologram

% save([indir, obj_name, '.mat'], 'holo', 'obj3d');
save([indir, obj_name, '_conv_holo', '.mat'], 'holo');

