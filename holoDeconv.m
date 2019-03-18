%{
----------------------------------------------------------------------------------------------------
Name:

Author:   Ni Chen (chenni@snu.ac.kr)
Date:     
Modified:

Reference:
- 
----------------------------------------------------------------------------------------------------
%}

close all; clear; clc; 

%% ====================================== Parameters =============================================== 
addpath('./function');  

indir = './data/';
outdir = './output/';

% Simulation: geo_overlap, random_scatter, helix_conical, helix_circular
% Experiment: dandelion, sh, beads
obj_name = 'helix_circular';  

deconv_type = 'TwISTL1';  % 'TwISTTV', 'TwISTL1'
holo_type = 'conv';             % 'LFFR', 'LFFT', 'conv'

iter_num = 400;  

isWrite2DRecon = 0;
isWriteEPS = 0;

%% ======================================= Hologram ================================================ 
% load([indir, obj_name, holo_type, '.mat']);
load([indir, obj_name, '_conv_holo', '.mat']);
run([indir, obj_name, '_param.m']);

% holo = abs(holo);

%% =================================================================================================  
[otf3d, psf3d, Pupil] = OTF3D(Nx, Ny, Nz, lambda, deltaX, deltaY, deltaZ, offsetZ, sensor_size);  

% Since the TV-based denoising term is not defined in a complex domain, the real and imaginary 
% parts are independently processed in a vectorized form.
field2dvec = C2V(holo(:));  

A_fun = @(vol3d) Propagation(vol3d, otf3d, Pupil);  
AT_fun = @(field2d) iPropagation(field2d, otf3d, Pupil);  

%% ========================== Reconstruction with back-propagation ================================= 
reobj_raw = AT_fun(field2dvec);  
reobj_raw = reshape(abs(V2C(reobj_raw)), Nx, Ny, Nz);  
% write3d(abs(reobj_raw), z_scope*1000, outdir, obj_name); 
%% ============================= Reconstruction with deconvolution ================================= 
Nx_prime = Nx; 
Ny_prime = Ny*Nz*2; 
Nz_prime = 1; 

tau = 0.01; 
piter = 4; 
tolA = 1e-6; 

switch deconv_type
    case 'Wiener'
        reobj_Wiener = 1;
        reobj_deconv = reobj_Wiener;
        
    case 'TwISTL1'  % Fast calculation and convergence, Better for no overlapping?
        Phi_fun = @(vol3d, weight, epsi) L1phi(vol3d);
        [reobj_TwIST, ~, obj_twist, times_twist, dummy, mse_twist] ...
            = TwIST(field2dvec, A_fun, tau, ...
                'AT', AT_fun, ...
                'Phi', Phi_fun, ...
                'Initialization', 2, ...
                'Monotone', 1, ...
                'StopCriterion', 1, ...
                'MaxIterA', iter_num, ...
                'MinIterA', iter_num, ...
                'ToleranceA', tolA, ...
                'Verbose', 1);  
         
         reobj_TwIST = reshape(V2C(reobj_TwIST), Nx, Ny, Nz);  
         reobj_deconv = reobj_TwIST;
         
    case 'TwISTTV'  
        Psi_fun = @(vol3d, th) TVpsi(vol3d, th, 0.05, piter, Nx_prime, Ny_prime, Nz_prime);  
        Phi_fun = @(vol3d) TVphi(vol3d, Nx_prime, Ny_prime, Nz_prime);
%         vol3d_vec = C2V(vol3d(:)); 
        [reobj_TwIST, ~, obj_twist, times_twist, ~, mse_twist] ...
            = TwIST(field2dvec, A_fun, tau, ...
                'AT', AT_fun, ...
                'Psi', Psi_fun, ...
                'Phi', Phi_fun, ...
                'Initialization', 2, ...
                'Monotone', 1, ...
                'StopCriterion', 1, ...
                'MaxIterA', iter_num, ...
                'MinIterA', iter_num, ...
                'ToleranceA', tolA, ...
                'Verbose', 1); ...
%                 'TRUE_X', vol3d_vec);

        reobj_TwIST = reshape(V2C(reobj_TwIST), Nx, Ny, Nz);  
        reobj_deconv = reobj_TwIST;
        
    case 'ADMM'
        reobj_deconv = reobj_ADMM;
end

save([outdir, obj_name, '_', deconv_type ,'_result.mat'], 'reobj_raw', 'reobj_deconv', 'iter_num', 'mse_twist');
 
%% ======================================= Show the images =========================================
% figure; imagesc(plotdatacube(real(otf3d)));title('OTF');axis image;drawnow;colorbar;
% print('-dpng', [outdir, obj_name, '_', 'otf', num2str(deltaZ), '.png']);
% figure; imagesc(plotdatacube(abs(psf3d)));title('PSF');axis image;drawnow;colorbar;
% figure; imagesc(plotdatacube(real(E)));title('E');axis image;drawnow;colorbar;
% figure; plot(mse_twist); xlabel('iteration number'); ylabel('MSE');

% Back Propagation reconstruction
figure; imagesc(plotdatacube(abs(reobj_raw))); title('BackPropagation'); axis image; drawnow; colorbar; axis off;
print('-dpng', [outdir, obj_name, '_', 'BP', '.png']);

% TwIST reconstruction
figure; imagesc(plotdatacube(abs(reobj_deconv))); title('Reconstruction'); axis image; drawnow; colorbar; axis off;
print('-dpng', [outdir, obj_name, '_', deconv_type, '_', num2str(iter_num), '.png']);


% % Show 3D OTF
% figure; show3d(abs(psf3d), 0.001); title('PSF');
% print('-dpng', [outdir, obj_name, '_PSF', num2str(deltaZ),'.png']);

% Show 3D volume
figure; show3d(abs(reobj_raw), 0.01); title('BackPropagation');
print('-dpng', [outdir, obj_name, '_BP_3d.png']);

figure; show3d(abs(reobj_deconv), 0.01); title('Reconstruction');
print('-dpng', [outdir, obj_name, '_', deconv_type, '_3d_', num2str(iter_num), '.png']);

% Write 2D images of the reconstruction
if isWrite2DRecon ~= 0
    % figure; imagesc(abs(holo)); title('Hologram'); axis image; colorbar;
    write3d(abs(reobj_raw), outdir, [obj_name, '_', deconv_type, '.png']); 
%     write3d(abs(reobj_deconv), outdir, [obj_name, '_', deconv_type, '.png']);
end