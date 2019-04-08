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

FT = @(x) ifftshift(fft(fftshift(x)));
iFT = @(x) ifftshift(ifft(fftshift(x)));

%% ====================================== Parameters ===============================================
addpath(genpath('./function/'));  % Add funtion path with sub-folders
indir = './data/';  % Hologram data

obj_name = 'twopoint_x';
holo_type = 'complex';  % complex; inline; offline;

% Deconvolution setting
iter_num = 200;
regu_type = 'TV';  % 'TV', 'L1'
deconv_type = 'TwIST';  % 'TwIST','GPSR', TVAL3, SALSA, NESTA, TVPD

outdir = './output/';  % Output files

%% ======================================= Hologram ================================================
run([indir, obj_name, '_param.m']);  % Parameters of the object and hologram
[otf3d, psf3d, pupil3d] = OTF3D_z(Ny, Nx, lambda, pps, z);

% figure; imagesc(plotdatacube(angle(otf3d))); title('OTF'); axis image; drawnow; colormap(hot); colorbar; axis off;
% figure; imagesc(plotdatacube(abs(psf3d))); title('PSF'); axis image; drawnow; colormap(hot); colorbar; axis off;
% figure; imagesc(plotdatacube(angle(pupil3d))); title('Pupil in Frequency'); axis image; drawnow; colormap(hot); colorbar; axis off;

% Simulation data
issim = 1;
load([indir, obj_name, '_3d.mat']);

vol3d_vec = C2V(obj3d(:));  % For calculating MSE
prop_field = MatProp3D(obj3d, otf3d, pupil3d);  % Field at the center plane of the 3D object

holo = prop_field;
tau = 0.1;   % This effects, need further investigation
tau_psi = 0.25;

holo = holoNorm(holo);
% holo = holo./max(abs(holo(:)));

% add noise
% holo = awgn(holo, 40);  % holo = imnoise(holo, 'gaussian', 0, 0.001);
out_filename = [outdir, obj_name, '_', holo_type , '_'];

%% ========================== Reconstruction with back-propagation =================================
reobj_raw = iMatProp3D(holo, otf3d, pupil3d, holo_type);
%% ============= Construct Minimization problem and Reconstruction with deconvolution ==============
A_fun = @(field3d_vec) VecProp3D(field3d_vec, otf3d, pupil3d, holo_type);
AT_fun = @(field2d_vec) iVecProp3D(field2d_vec, otf3d, pupil3d, holo_type);

Nyv = Ny;
Nxv = Nx*Nz*2;  %!!!!
Nzv = 1;

% Since the TV-based denoising term is not defined in a complex domain, the real and imaginary parts
% are independently processed in a vectorized form.
holo_vec = C2V(holo(:));

tolA = 1e-6;

Phi_fun = @(vol3d) TVphi(vol3d, Nyv, Nxv, Nzv);
piter = 5;
Psi_fun = @(vol3d, th) TVpsi(vol3d, th, tau_psi, piter, Nyv, Nxv, Nzv);

[reobj_TwIST, ~, objective, times, ~, mses]= TwIST(holo_vec, A_fun, tau, ...
    'AT', AT_fun, ...
    'Phi', Phi_fun, ...
    'Psi', Psi_fun,...
    'MaxIterA', iter_num, ...
    'ToleranceA', tolA, ...
    'TRUE_X', vol3d_vec);

reobj_TwIST = reshape(V2C(reobj_TwIST), Ny, Nx, Nz);
reobj_deconv = reobj_TwIST;


%% ======================================= Show the images =========================================
if strcmp(obj_name , 'twopoint_z')
    % Back vecProp reconstruction
    raw_temp = abs(reobj_raw);
    raw_temp = (selectline(raw_temp,x_shift,z_shift));
    raw_temp = mat2gray(raw_temp);
    
    % TwIST reconstruction
    deconv_temp = (abs(reobj_deconv));
    deconv_temp = selectline(deconv_temp,x_shift,z_shift);
    deconv_temp = mat2gray(deconv_temp);
    
    figure;
    plot(z_axis, selectline(obj3d, x_shift, z_shift), 'k', 'LineWidth', 1);
    hold on; plot(z_axis, raw_temp, 'b', 'LineWidth', 1.5);
    hold on; plot(z_axis, deconv_temp, 'r-.', 'LineWidth',1.5);
    
    zf = round(Nz/2);
    plot([z_axis(zf-z_shift),z_axis(zf-z_shift)],[0,1],'k--','LineWidth',1);
    plot([z_axis(zf+z_shift),z_axis(zf+z_shift)],[0,1],'k--','LineWidth',1);
    
    plot([z_axis(zf+z_shift),z_axis(zf-z_shift)],[0,1],'g','LineWidth',1);
    % plot([z_axis(27),z_axis(27)],[0,1],'k:','LineWidth',1);
    
    lgnd = legend('Original points', 'Back propagation', 'Proposed method', 'FontSize',13);
    xlim([min(z_axis) max(z_axis)]);
    set(lgnd, 'Box', 'off', 'color', 'none'); xlabel('z (mm)'); ylabel('Normalized intensity');
    set(gcf,'Position',[0,0,500,300]); set(gca,'FontSize',14)
    print('-depsc', [outdir, '2point_zplot', '.eps']);
    
    %======================================================================
    figure;
    plot(z_axis, selectline(obj3d, x_shift, z_shift), 'b', 'LineWidth', 1);
    hold on;
    plot(z_axis, deconv_temp, 'r-.', 'LineWidth',1.5);
    
    zf = round(Nz/2);
    plot([z_axis(zf-z_shift),z_axis(zf-z_shift)],[0,1],'k--','LineWidth',1);
    plot([z_axis(zf+z_shift),z_axis(zf+z_shift)],[0,1],'k--','LineWidth',1);
    
    plot([z_axis(zf+z_shift),z_axis(zf-z_shift)],[0,1],'g','LineWidth',1);
    % plot([z_axis(27),z_axis(27)],[0,1],'k:','LineWidth',1);
    
    lgnd = legend('Original points', 'Proposed method', 'FontSize',13);
    xlim([min(z_axis) max(z_axis)]);
    set(lgnd, 'Box', 'off', 'color', 'none');
    xlabel('z (mm)');ylabel('Normalized intensity');
    set(gcf,'Position',[0,0,500,300]); set(gca,'FontSize',14)
    print('-depsc', [outdir, '2point_zplot_dz', num2str(dz_factor), '.eps'])
       
%     eval(['2point_dz', num2str(dz_factor), '=deconv_temp;']);
%     save('2point_zplot_deconv.mat', eval(['2point_dz', num2str(dz_factor)]),'-append');
    
    %======================================================================
%     raw_temp = squeeze(abs(reobj_raw(Ny, :, :)));
%     raw_temp = mat2gray(raw_temp);
%     deconv_temp = squeeze(abs(reobj_deconv(Ny, :, :)));
%     raw_temp = mat2gray(raw_temp);
%     
%     figure;
%     imshow(raw_temp,[],'Border','tight', 'InitialMagnification', 100); colormap(hot); 
% %     set(gcf,'Position',[0,0,200,200]); 
%     set(gca,'LooseInset', get(gca,'TightInset')); print([outdir,'2point_raw.png'],'-dpng','-r300');
%     figure;
%     imshow(deconv_temp,[],'Border','tight', 'InitialMagnification', 100); colormap(hot); 
% %     set(gcf,'Position',[0,0,200,200]); 
%     set(gca,'LooseInset', get(gca,'TightInset')); print([outdir,'2point_twist.png'],'-dpng','-r300');
    
else
    raw_temp = abs(reobj_raw);
    raw_temp = mat2gray(raw_temp);
    figure;
    imshow(raw_temp,[],'Border','tight', 'InitialMagnification', 100);colormap(hot);    
    
    deconv_temp = abs(reobj_deconv);
    deconv_temp = mat2gray(deconv_temp);
    figure;
    imshow(deconv_temp,[],'Border','tight', 'InitialMagnification', 100); colormap(hot); 
    
    figure;
    plot(x_axis, obj3d(Ny/2,:), 'k');hold on;
    plot(x_axis, raw_temp(Ny/2,:), 'b--');hold on;
    plot(x_axis, deconv_temp(Ny/2,:), 'r--');
    
    xlim([min(x_axis) max(x_axis)]);
    legend('Original', 'Back propagation', 'Proposed method');
    xlabel('x'); ylabel('y');    
    print('-depsc', [outdir, '2point_xplot', '.eps']);
end