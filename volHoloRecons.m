%{
----------------------------------------------------------------------------------------------------
Name: Volume hologram reconstruction with deconvolution

Author:   Ni Chen (chenni@snu.ac.kr)
Date:
Modified:

Reference:
-
----------------------------------------------------------------------------------------------------
%}

close all; clear; clc;

FT2 = @(x) ifftshift(fft2(fftshift(x)));
iFT2 = @(x) ifftshift(ifft2(fftshift(x)));
%% ====================================== Parameters ===============================================
addpath('./function/');

indir = './data/';
outdir = './output/';

% Simulation: geo, geo_overlap, random_scatter, helix_conical, helix_circular, SNUE
% Experiment: dandelion, sh, beads, res
obj_name = 'SNUE';
issim = 1; % Simulation

holo_type = 'inline';  % complex (OSH, phase-shifting); inline; offline;

regu = 'TV'; % 'TV', 'L1'
deconv_type = 'TwIST';  % 'TwIST', 'GPSR', SALSA

iter_num = 200;

isWrite2DRecon = 0;
isWriteEPS = 0;

%% ======================================= Hologram ================================================
load([indir, obj_name, '.mat']);
if(issim)
    load([indir, obj_name, '_conv_holo.mat']); 
end
run([indir, obj_name, '_param.m']);

holo = holo./max(max(abs(holo)));

% Filter
% sigma = 1;
% gausFilter = fspecial('gaussian', [3,3], sigma);
% holo = imfilter(holo, gausFilter, 'replicate');

switch holo_type
    case 'inline'
        %{
          For inline hologram, I = |R|^2 + OR* + O*R + |O|^2, |R|^2 can be captured directly, and
          eliminated by I_prime = iFT(FT(I) - FT(|R|^2 )).
        %}
%         holo = holo + conj(holo) + sum(obj3d.^2,3)/Nz;  % hologram without |R|^2, suppose |R|=1
        holo = holo + conj(holo);  % in weak object case, = 2real(O),
    case 'offline'
end

%% ============================== Construst Minimization problem ===================================
[otf3d, psf3d, Pupil] = OTF3D(Ny, Nx, Nz, lambda, deltaX, deltaY, deltaZ, offsetZ, sensor_size);
A_fun = @(vol3d) Propagation(vol3d, otf3d, Pupil, holo_type);
AT_fun = @(field2d) iPropagation(field2d, otf3d, Pupil, holo_type);

%% ========================== Reconstruction with back-propagation =================================
reobj_raw = iPropagation3D(holo, otf3d, Pupil, holo_type);
% write3d(abs(reobj_raw), z_scope*1e9, outdir, obj_name);

%% ============================= Reconstruction with deconvolution =================================
Nyv = Ny;
Nxv = Nx*Nz*2;  %!!!!
Nzv = 1;

% Since the TV-based denoising term is not defined in a complex domain, the real and imaginary parts
% are independently processed in a vectorized form.
holo_vec = C2V(holo(:));

if issim
    vol3d_vec = C2V(obj3d(:));
else
    temp = zeros(Ny, Nx, Nz);
    vol3d_vec = C2V(temp(:));
end

switch deconv_type
    case 'TwIST'  % gold but slowly
        tau = 0.05;   % need an estamation function
        tolA = 1e-6;
        
        if strcmp(regu, 'L1')
            Phi_fun = @(vol3d, weight, epsi) L1phi(vol3d);
            
            [reobj_TwIST, ~, objective, times, ~, mses] = TwIST(holo_vec, A_fun, tau, ...
                'AT', AT_fun, ...
                'Phi', Phi_fun, ...
                'MaxIterA', iter_num, ...
                'ToleranceA', tolA, ...
                'TRUE_X', vol3d_vec);
            
        elseif strcmp(regu, 'TV')
            Phi_fun = @(vol3d) TVphi(vol3d, Nyv, Nxv, Nzv);
            piter = 5;
            Psi_fun = @(vol3d, th) TVpsi(vol3d, th, tau, piter, Nyv, Nxv, Nzv);
            
            [reobj_TwIST, ~, objective, times, ~, mses]= TwIST(holo_vec, A_fun, tau, ...
                'AT', AT_fun, ...
                'Phi', Phi_fun, ...
                'Psi', Psi_fun,...
                'MaxIterA', iter_num, ...
                'ToleranceA', tolA, ...
                'TRUE_X', vol3d_vec);
        end
        
        reobj_TwIST = reshape(V2C(reobj_TwIST), Ny, Nx, Nz);
        reobj_deconv = reobj_TwIST;
        
    case 'GPSR' % works but not for overlapping
        tau = 0.01;
        tolA = 1e-6;
        [reobj_SALSA, objective, obj_gpsr, times, debias_s, mses] = GPSR_Basic(holo_vec, A_fun, tau,...
            'Debias', 1, ...
            'AT', AT_fun,...
            'ToleranceA', tolA, ...
            'MaxIterA', iter_num, ...
            'TRUE_X', vol3d_vec);
        
        reobj_SALSA = reshape(V2C(reobj_SALSA), Ny, Nx, Nz);
        reobj_deconv = reobj_SALSA;
    case 'SALSA'  % not working
        A_fun = @(vol3d) Propagation(vol3d, otf3d, Pupil, holo_type);
        AT_fun = @(field2d) iPropagation(field2d, otf3d, Pupil, holo_type);
        
        %{
        reg parameter, extreme eigenvalues, rule of thumb:
        1e-4: severyly ill-conditioned problems
        1e-1: mildly  ill-conditioned problems
        1:    when A = Unitary matrix
        %}
        tau = 1e-0;
        mu = tau/5;
        
        tolA = 1e-3; % Smoothing parameter (empirical setting)
        
        otf3dvec = C2V(otf3d(:));
        filter_FFT = 1./(abs(otf3dvec) + mu);
        
        invLS = @(x) real(ifftshift(ifft(fftshift(filter_FFT.*ifftshift(fft(fftshift(x)))))));  % ????
        invLS = @(x) callcounter(invLS, x);
        
        Psi_fun = @(vol3d, th) TVpsi(vol3d, tau, tolA, 10, Nyv, Nxv, Nzv);
        Phi_fun = @(vol3d) TVphi(vol3d, Nyv, Nxv, Nzv);
        
        holo_vec = C2V(holo(:));
        [reobj_salsa, numA, numAt, objective, distance, times] = SALSA(holo_vec, A_fun, tau, ...
            'MU', mu, ...
            'AT', AT_fun, ...
            'ToleranceA', tolA,...
            'MaxIterA', iter_num, ...
            'Psi', Psi_fun, ...
            'Phi', Phi_fun, ...
            'LS', invLS, ...
            'True_x', vol3d_vec);
        
        reobj_salsa = reshape(V2C(reobj_salsa), Ny, Nx, Nz);
        reobj_deconv = reobj_salsa;
        
    case 'CSALSA' % not working
        A_fun = @(vol3d) Propagation(vol3d, otf3d, Pupil, holo_type);
        AT_fun = @(field2d) iPropagation(field2d, otf3d, Pupil, holo_type);
        
        tolA = 1e-10;
        
        mu1 = 1;
        mu2 = mu1;
        %         epsilon = sqrt(N^2+8*sqrt(N^2))*sigma;
        epsilon = 1e-6;
        %         Pnum = norm(x(:)-y(:),2)^2;
        Pnum = 0.1;
        tol = 1e-6;
        sigma = 0.1;
        
        otf3dvec = C2V(otf3d(:));
        
        H_FFT = ifftshift( fft(fftshift(otf3dvec)));
        H2 = abs(H_FFT).^2;
        tau = mu1/mu2;
        filter_FFT = H2./(H2 + tau);
        %         invLS = @(x, mu1)(1/mu1)*( x - real(ifft( filter_FFT.*fft( x ) ) ) );
        invLS = @(x, mu1)(1/mu1)*( x - real(ifftshift(ifft(fftshift(filter_FFT.*ifftshift(fft(fftshift(x)))) ) )) );
        LS = @(x,mu) callcounter(invLS,x,mu);
        
        piter = 10;
        Psi_fun = @(vol3d, th) TVpsi(vol3d, th, tau, piter, Nyv, Nxv, Nzv);
        Phi_fun = @(vol3d) TVphi(vol3d, Nyv, Nxv, Nzv);
        
        holo_vec = C2V(holo(:));
        [reobj_csalsa, numA, numAt, objective, distance1, distance2, criterion, times, mses] = ...
            CSALSA(holo_vec, A_fun, mu1, mu2, sigma,...
            'AT', AT_fun, ...
            'PHI', Phi_fun, ...
            'PSI', Psi_fun, ...
            'StopCriterion', 3, ...
            'True_x', vol3d_vec, ...
            'ToleranceA', tol,...
            'MaxIterA', iter_num, ...
            'LS', LS, ...
            'VERBOSE', 1, ...
            'EPSILON', epsilon, ...
            'INITIALIZATION', 2, ...
            'TVINITIALIZATION', 0, ...
            'CONTINUATIONFACTOR', 0);
        
        reobj_csalsa = reshape(V2C(reobj_csalsa), Ny, Nx, Nz);
        reobj_deconv = reobj_csalsa;
        
    otherwise
        
end

% save([outdir, obj_name, '_', deconv_type ,'_result.mat'], 'reobj_raw', 'reobj_deconv', 'iter_num', 'mse_twist');

%% ======================================= Show the images =========================================
% figure; imagesc(plotdatacube(real(otf3d)));title('OTF');axis image;drawnow;colorbar;
% print('-dpng', [outdir, obj_name, '_', 'otf', num2str(deltaZ), '.png']);
% figure; imagesc(plotdatacube(abs(psf3d)));title('PSF');axis image;drawnow;colorbar;
% figure; imagesc(plotdatacube(real(E)));title('E');axis image;drawnow;colorbar;
figure; semilogy(times, objective, 'LineWidth',2); ylabel('objective distance'); xlabel('CPU time (sec)');
figure; plot(times, mses, 'LineWidth',2); ylabel('MSE'); xlabel('CPU time (sec)');

pause(1);

% Back Propagation reconstruction
figure; imagesc(plotdatacube(abs(reobj_raw))); title('BackPropagation'); axis image; drawnow; colormap(hot); colorbar; axis off;
print('-dpng', [outdir, obj_name, '_', 'BP', '.png']);

% TwIST reconstruction
figure; imagesc(plotdatacube(abs(reobj_deconv))); title('Reconstruction'); axis image; drawnow; colormap(hot); colorbar; axis off;
print('-dpng', [outdir, obj_name, '_', deconv_type, '_', num2str(iter_num), '.png']);

% % Show 3D OTF
% figure; show3d(abs(psf3d), 0.001); title('PSF');
% print('-dpng', [outdir, obj_name, '_PSF', num2str(deltaZ),'.png']);

% Show 3D volume
% figure; show3d(abs(reobj_raw), 0.01); title('BackPropagation');
% print('-dpng', [outdir, obj_name, '_BP_3d.png']);

figure; show3d(abs(reobj_deconv), 0.01); title('Reconstruction');
print('-dpng', [outdir, obj_name, '_', deconv_type, '_3d_', num2str(iter_num), '.png']);
% write3d(abs(reobj_deconv), z_scope*1e9, outdir, [obj_name, '_' ,deconv_type]);

% Write 2D images of the reconstruction
if isWrite2DRecon ~= 0
    % figure; imagesc(abs(holo)); title('Hologram'); axis image; colorbar;
    write3d(abs(reobj_raw), outdir, [obj_name, '_', deconv_type, '.png']);
    %     write3d(abs(reobj_deconv), outdir, [obj_name, '_', deconv_type, '.png']);
end