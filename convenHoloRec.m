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

% geo, random_scatter, geo_overlap, helix_circular, helix_conical, SNUE
obj_name = 'SNUE'; 
run([indir, obj_name, '_param.m']);
load([indir, obj_name, '.mat']);

%% ========================================= Propagation ===========================================
[otf3d, psf3d, pupil3d] = OTF3D(Ny, Nx, Nz, lambda, deltaX, deltaY, deltaZ, offsetZ, sensor_size);
prop_field = Propagation3D(obj3d, otf3d, pupil3d);  % Field at the center plane of the 3D object
holo = prop_field;

% save([indir, obj_name, '.mat'], 'holo', 'obj3d');
save([indir, obj_name, '_conv_holo', '.mat'], 'holo');

