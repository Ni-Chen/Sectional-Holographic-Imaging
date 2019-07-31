
%{
----------------------------------------------------------------------------------------------------
Name: Volume hologram reconstruction with deconvolution

Author:   Ni Chen (chenni@snu.ac.kr)
Date:
Modified:

Reference:
-`
----------------------------------------------------------------------------------------------------
%}

%% Initialization
close all; clear; clc;

set(0, 'DefaultAxesFontWeight', 'normal', ...
      'DefaultAxesFontAngle', 'normal', ... % Not sure the difference here
      'DefaultAxesFontWeight', 'normal', ... % Not sure the difference here
      'DefaultAxesTitleFontWeight', 'normal', ...
      'DefaultAxesTitleFontSizeMultiplier', 1,...
      'DefaultAxesLooseInset',[0,0,0,0]);

FT2 = @(x) ifftshift(fft2(fftshift(x)));
iFT2 = @(x) ifftshift(ifft2(fftshift(x)));

FT = @(x) ifftshift(fft(fftshift(x)));
iFT = @(x) ifftshift(ifft(fftshift(x)));

addpath(genpath('./function/'));  % Add funtion path with sub-folders

indir = './data/';  % Hologram data

% Simulations: random, geo, overlap, cirhelix, conhelix, SNUE, DISP
% Experiments: dandelion, sh, beads, res, hair
obj_name = 'mark';
holo_type = 'complex';  % complex; inline; offline;

% Output setting
isDebug = 1;

% Deconvolution setting
iter_num = 200;
regu_type = 'TV';  % 'TV', 'L1'
deconv_type = 'TwIST';  % 'TwIST','GPSR', TVAL3, SALSA, NESTA, TVPD

if(isDebug)
    outdir = './output/';  % Output files
else
    outdir = './Complex/figure/';  % Output images for paper writing
%     outdir = './Sectional/figure/';  % Output images for paper writing
end
%% 
%% Hologram Generation / loading (experimental) 

run([indir, obj_name, '_param.m']);  % Parameters of the object and hologram 
% [otf3d, psf3d, pupil3d] = OTF3D(Ny, Nx, Nz, lambda, deltaX, deltaY, deltaZ, offsetZ, sensor_size);
[otf3d, psf3d, pupil3d] = OTF3D(Ny, Nx, lambda, pps, z);


if any(strcmp(obj_name, {'geo', 'overlap', 'random', 'conhelix', 'cirhelix', 'SNUE', 'finger', 'DISP'}))
% Simulation data
    issim = 1;    
    load([indir, obj_name, '_3d.mat']);
%     figure; imagesc(plotdatacube(abs(obj3d))); 

    vol3d_vec = C2V(obj3d(:));  % For calculating MSE
    
    prop_field = MatProp3D(obj3d, otf3d, pupil3d);  % Field at the center plane of the 3D object
     
    switch holo_type
        case 'complex'
            holo = prop_field;
            tau = 0.1;   % This effects, need further investigation
            tau_psi = 0.25;           
        
            holo = holoNorm(holo);
            
%         case 'inline'
%             %{
%               Inline hologram, I = |R|^2 + OR* + O*R + |O|^2, |R|^2 can be captured directly, 
%               OR* + O*R + |O|^2 = iFT(FT(I) - FT(|R|^2 )), suppose |R|=1, real(OR* + O*R)~=2real(O)
%             %}
%             holo = prop_field + conj(prop_field) + sum(obj3d.^2,3)/Nz;  % Nz effects much???
% %             holo = prop_field + conj(prop_field) + sum(obj3d.^2,3); 
% %             holo = prop_field  + sum(obj3d.^2,3);    % Compressive holography
% 
%             holo = holo./max(abs(holo(:)));
% %             holo = holoNorm(holo);
%             tau = 0.01;   % This effects, need further investigation
%             tau_psi = 0.1;
%         case 'offline'
    end
    
    % add noise
%     holo = awgn(holo, 40);  % holo = imnoise(holo, 'gaussian', 0, 0.001);
    
    out_filename = [outdir, obj_name, '_', holo_type , '_'];
else
% Experimental data
    issim = 0;
    tau = 0.005;   % This effects, need further investigation
    tau_psi = 0.05;
    temp = zeros(Ny, Nx, Nz);
    vol3d_vec = C2V(temp(:));  % Experiments have no reference object, just for coding convenience
%     vol3d_vec = C2V(temp);  % Experiments have no reference object, just for coding convenience
   
    out_filename = [outdir, obj_name, '_'];
    
%     holo = holodenoise(holo);
    
%     holo = holo./max(abs(holo(:)));
%     holo = max(holo(:))-holo;
    

     holo = holoNorm(holo);
      
    
    clear C C_F;
end

% [sigma, mu]= noiseAnalysis(holo)
% holo = holodenoise(holo); % Filter noise


%% Reconstruction with back-propagation 
reobj_raw = iMatProp3D(holo, otf3d, pupil3d);
temp = reobj_raw;

% write3d(abs(temp), z*1e6, outdir, obj_name);
% figure; imagesc(plotdatacube(abs(temp))); title('BackPropagation'); axis image; colormap(hot); axis off;
% figure; imagesc(plotdatacube(angle(pupil3d))); title('BackPropagation'); axis image; drawnow; colormap(hot); colorbar; axis off;
%% 
%% Reconstruction with deconvolution
% Construct Minimization problem 

A_fun = @(field3d_vec) VecProp3D(field3d_vec, otf3d, pupil3d);
AT_fun = @(field2d_vec) iVecProp3D(field2d_vec, otf3d, pupil3d);

Nyv = Ny;
Nxv = Nx*Nz*2;  %!!!!
Nzv = 1;

% Nyv = Ny*2;
% Nxv = Nx;  %!!!!
% Nzv = Nz;


% Since the TV-based denoising term is not defined in a complex domain, the real and imaginary parts

% %% Use spacial filter to get noise out of the measurement
% N = Nx;
% C = zeros(N);
% % N = size(C);
% %C is a diagonal regularization matrix, a 'laplacian' filter
% C(ceil(N/2)-1:ceil(N/2)+1,ceil(N/2)-1:ceil(N/2)+1) = fspecial('laplacian',0.1);
% fC = fft2(C); %two-dimensional discrete Fourier transform for filter C
% fC = fC/max(max(abs(fC)));
% fC_holo = im2col(1-fC,[N N],'distinct');
% % d_C = im2col(conj(fC).*fC,[N N],'distinct');
% d_C = conj(fC).*fC;
% holo = holo + tau*d_C;

% are independently processed in a vectorized form.
holo_vec = C2V(holo(:));
% holo_vec = C2V(holo);

% Solving Minimization problem withdifferent deconvolution algorithms

switch deconv_type
    case 'TwIST'  % works for both complex and inline holograms, but very slowly
%         tau = 0.01;   % This effects, need further investigation
%         tau_psi = 0.25;
        tolA = 1e-6;
        
%         if strcmp(regu_type, 'L1')
%             Phi_fun = @(vol3d, weight, epsi) L1phi(vol3d);
%             
%             [reobj_TwIST, ~, objective, times, ~, mses] = TwIST(holo_vec, A_fun, tau, ...
%                 'AT', AT_fun, ...
%                 'Phi', Phi_fun, ...
%                 'MaxIterA', iter_num, ...
%                 'ToleranceA', tolA, ...
%                 'TRUE_X', vol3d_vec);
%         
%         elseif strcmp(regu_type, 'TV')
            Phi_fun = @(vol3d) TVphi(vol3d, Nyv, Nxv, Nzv);
            piter = 5;
            Psi_fun = @(vol3d, th) TVpsi(vol3d, th, tau_psi, piter, Nyv, Nxv, Nzv);
            
            [reobj_TwIST, ~, objective, times, ~, mses]= TwIST(holo_vec, A_fun, tau, ...
                'AT', AT_fun, ...
                'Phi', Phi_fun, ...
                'Psi', Psi_fun,...
                'MaxIterA', iter_num, ...
                'ToleranceA', tolA);
%                 'TRUE_X', vol3d_vec);
%         end
        
        reobj_TwIST = reshape(V2C(reobj_TwIST), Ny, Nx, Nz);
        reobj_deconv = reobj_TwIST;
        
    case 'TwIST_GPSR'  % works for both complex and inline holograms, but very slowly
        tau = 0.005;  % regularization parameter
        tolA = 1e-6; % stopping threshold
        [reobj_GPSR, reobj_debias, obj_gpsr, times, debias_s, mses] = GPSR_BB(holo_vec, A_fun, tau,...
            'Debias', 1, ...
            'AT', AT_fun,...
            'Continuation', 0, ...
            'ToleranceA', tolA, ...
            'MaxIterA', 20, ...   % 15
            'TRUE_X', vol3d_vec);
        
        tau = 0.01;   % This effects, need further investigation
        tau_psi = 0.2;
        Phi_fun = @(vol3d) TVphi(vol3d, Nyv, Nxv, Nzv);
        piter = 5;
        Psi_fun = @(vol3d, th) TVpsi(vol3d, th, tau_psi, piter, Nyv, Nxv, Nzv);

        [reobj_TwIST, ~, objective, times, ~, mses]= TwIST(holo_vec, A_fun, tau, ...
            'AT', AT_fun, ...
            'Phi', Phi_fun, ...
            'Psi', Psi_fun,...
            'Initialization', reobj_GPSR,...
            'MaxIterA', iter_num, ...
            'ToleranceA', tolA, ...
            'TRUE_X', vol3d_vec);
        
        reobj_TwIST = reshape(V2C(reobj_TwIST), Ny, Nx, Nz);
        reobj_deconv = reobj_TwIST;
        
    case 'GPSR' % works but not for overlapping
        tau = 0.005;  % regularization parameter
        tolA = 1e-6; % stopping threshold        
        
        [reobj_SALSA, reobj_debias, obj_gpsr, times, debias_s, mses] = GPSR_BB(holo_vec, A_fun, tau,...
            'Debias', 1, ...
            'AT', AT_fun,...
            'Continuation', 0, ...
            'ToleranceA', tolA, ...
            'MaxIterA', iter_num, ...
            'TRUE_X', vol3d_vec);
        
        reobj_SALSA = reshape(V2C(reobj_SALSA), Ny, Nx, Nz);
        reobj_deconv = reobj_SALSA;
        
  
    otherwise        
end
%% 
%% Show the images 

if isDebug
    figure; semilogy(objective, 'LineWidth',2); ylabel('objective distance'); xlabel('Iteration');
%     print('-dpng', [out_filename, 'objdis.png']);
    
%     figure; plot(mses, 'LineWidth',2); ylabel('MSE'); 
        
    % Back vecProp reconstruction
    figure; imagesc(plotdatacube(abs(reobj_raw))); title('BackPropagation'); axis image; drawnow; colormap(hot);  axis off;
    print('-dpng', [out_filename, 'BP', '.png']);
    
    % TwIST reconstruction
    temp = abs(reobj_deconv);
    temp = (temp-min(temp(:)))/(max(temp(:))-min(temp(:)));
    figure; imagesc(plotdatacube(abs(temp))); title('Deconvolution'); axis image; drawnow; colormap(hot);  axis off;
    print('-dpng', [out_filename, deconv_type, '_tau', num2str(tau), '_psi' , num2str(tau_psi), '_', num2str(iter_num), '.png']);
    
%     figure; show3d(temp, 0.01); axis normal; set(gcf,'Position',[0,0,500,500]);
    
%     figure; show3d(abs(reobj_deconv), 0.01); axis normal; set(gcf,'Position',[0,0,500,500]);
%     print('-dpng', [out_filename, deconv_type,  '_', num2str(iter_num), '_3d.png']);
%     view(0, 0); colorbar off; print('-dpng', [out_filename, deconv_type,  '_', num2str(iter_num), '_side.png']);
%     view(0, 90); colorbar off; print('-dpng', [out_filename, deconv_type, '_', num2str(iter_num), '_top.png']);
    
%     write3d(abs(reobj_deconv), z_scope*1e6, outdir, [obj_name, '_' ,deconv_type]);
   
    % For turnning tau
%     figure; show3d(abs(reobj_deconv), 0.01); axis normal; set(gcf,'Position',[0,0,500,500]);
%     saveas(gca, [outdir, obj_name, '_', num2str(tau), '_' , num2str(tau_psi), '_3d.png'], 'png');
%     view(0, 0); colorbar off; 
%     saveas(gca, [outdir, obj_name, '_', num2str(tau), '_' , num2str(tau_psi),'_side.png'], 'png'); 
%     view(0, 90); colorbar off; 
%     saveas(gca, [outdir, obj_name, '_', num2str(tau), '_' , num2str(tau_psi),'_top.png'], 'png');
    
else
    save([outdir, obj_name, '_', deconv_type ,'_result.mat'], 'reobj_raw', 'reobj_deconv', 'iter_num');

    figure; plot(mses(10:end), 'LineWidth',2); ylabel('MSE'); xlabel('Iteration');
      
    % write3d(abs(reobj_deconv), z_scope*1e9, outdir, obj_name);
    
    imwrite(mat2gray(real(holo)), [out_filename, 'holo_real.png']);
    imwrite(mat2gray(imag(holo)), [out_filename, 'holo_img.png']);
    
    % Show 3D volume
%     figure; show3d(abs(reobj_raw), 0.01); axis normal;set(gcf,'Position',[0,0,500,500]);
% %     caxis([0 1]);
%     saveas(gca, [out_filename, 'BP_3d.png'], 'png');
%     view(0, 0); colorbar off; axis equal; 
%     saveas(gca, [out_filename, 'BP_side.png'], 'png');   
%     view(0, 90); colorbar off; axis equal; 
%     saveas(gca, [out_filename, 'BP_top.png'], 'png'); 

%     temp = (reobj_deconv-min(reobj_deconv(:)))/(max(reobj_deconv(:))-min(reobj_deconv(:)));
    temp = abs(reobj_deconv);
    temp = (temp-min(temp(:)))/(max(temp(:))-min(temp(:)));
    figure; show3d(temp, 0.01); axis normal; set(gcf,'Position',[0,0,500,500]);
    caxis([0 1]);
    saveas(gca, [out_filename, deconv_type, '_3d.png'], 'png');
%     view(0, 0); colorbar off; 
%     saveas(gca, [out_filename, deconv_type, '_side.png'], 'png'); 
%     view(0, 90); colorbar off; 
%     saveas(gca, [out_filename, deconv_type, '_top.png'], 'png');
    
%     figure; show3d(abs(obj3d), 0.01); axis normal; set(gcf,'Position',[0,0,500,500]);
%     caxis([0 1]);
%     saveas(gca,[outdir, obj_name, '.png'], 'png');
%     view(0, 0); colorbar off; 
%     saveas(gca,[outdir, obj_name, '_side.png'], 'png');
%     view(0, 90); colorbar off;
%     saveas(gca,[outdir, obj_name, '_top.png'], 'png');
    mses(end)
%     set(gcf,'paperpositionmode','auto');
%     print('-depsc', [out_filename, deconv_type, '_3d.eps']);
end