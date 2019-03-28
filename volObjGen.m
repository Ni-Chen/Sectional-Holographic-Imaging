%{
----------------------------------------------------------------------------------------------------
Name: Code for generate 3D object.

Author   : Ni Chen (chenni@snu.ac.kr)
Date     :     
Modified :

----------------------------------------------------------------------------------------------------
%}

close all;
clear;
clc;

addpath('./function');
addpath('./data/');

indir = './data/';
outdir = './output/';

% random, geo, overlap, cirhelix, conhelix, SNUE
obj_name = 'random'; 
run([indir, obj_name, '_param.m']);

%% ====================================== 3D object ================================================
switch obj_name
     case 'SNUE'
        % Transmitance should be 0~1
        S = mat2gray((imread('S.tif')));   
        N = mat2gray((imread('N.tif')));
        U = mat2gray((imread('U.tif')));
        E = mat2gray((imread('E.tif')));
         
        obj3d = zeros(Nx, Ny, Nz);
        obj3d(:,:,2)=0.6*S;
        obj3d(:,:,3)=0.8*N;
        obj3d(:,:,4)=0.9*U;
        obj3d(:,:,5)=1*E;
        
     case 'geo'
        obj3d = zeros(Nx, Ny, Nz);
        
        obj3d(40:60,40:50, 15) = 0.6;
        obj3d(22:42,22:42, 30) = 0.8;
        obj3d(5:25,5:22, 55) = 1;
        
    case 'random'
        obj3d = zeros(Nx, Ny, Nz);
        for iz = 1:Nz
            ix = ceil(rand(1)*Nx);
            iy = ceil(rand(1)*Ny);
            
            r=0;
            
            ix_range = ix:ix+r;
            iy_range = iy:iy+r;
            
            ix_range(ix_range>Nx) = Nx;
            iy_range(iy_range>Ny) = Ny;
            
            obj3d(iy_range, ix_range, iz) = 1;     % rand(1)        
        end
        
    case 'overlap'
        obj3d = zeros(Nx, Ny, Nz);
        obj3d(20:50,20:50, 55) = 0.8;
        obj3d(30:35,10:55, 30) = 0.9;
        obj3d(10:45,15:25, 15) = 1.0;
   
    case 'cirhelix' %  circular helix
        obj3d = zeros(Ny,Nx,Nz);
        a = 1.5;
        c = 0.8;
        t = 0:0.01:12*pi;

        ix = 0.1;
        iy = ix;
        iz = ix;

        xx = round(a.*sin(t)/ix) + Nx/2;
        yy = round(a.*cos(t)/iy) + Ny/2;
        zz = round(t/(2*pi*c)/iz) + 1;

        xyz = sub2ind(size(obj3d), xx, yy, zz);

        obj3d(xyz) = 1;
        
    case 'conhelix'  % conical helix 
        a = 0.25;
        c = 1.2;
        t = 0:0.005:12*pi;     
        
        obj3d = zeros(Ny,Nx,Nz);
        
        ix = 0.025*2;
        iy = ix;
        iz = 0.035*2;
        
        xx = round(a*t.*sin(t)/(2*pi*c)/ix) + Nx/2;
        yy = round(a*t.*cos(t)/(2*pi*c)/iy) + Ny/2;
        zz = round(t/(2*pi*c)/iz) + 1;
        
        xyz = sub2ind(size(obj3d),xx,yy,zz);
        
        obj3d(xyz) = 1;
        
end

%% =================================================================================================
figure;
show3d(obj3d, 0.01);
title('Object');
print('-dpng', [outdir, obj_name, '.png']);
% print('-dsvg', [outdir, obj_name, '.svg']);
% fig2svg([outdir, obj_name, '.svg']);

save([indir, obj_name, '_3d.mat'], 'obj3d');

