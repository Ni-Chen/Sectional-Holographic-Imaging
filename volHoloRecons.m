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

% Simulations: random, geo, overlap, cirhelix, conhelix, SNUE
% Experiments: dandelion, sh, beads, res, hair
obj_name = 'beads';
holo_type = 'complex';  % complex; inline; offline;

% Output setting
isDebug = 1;

% Deconvolution setting
iter_num = 1000;
regu_type = 'TV';  % 'TV', 'L1'
deconv_type = 'TwIST';  % 'TwIST','GPSR', TVAL3, SALSA, NESTA, TVPD

if(isDebug)
    outdir = './output/';  % Output files
else
    outdir = './Complex/figure/';  % Output images for paper writing
%     outdir = './Sectional/figure/';  % Output images for paper writing
end

%% ======================================= Hologram ================================================
run([indir, obj_name, '_param.m']);  % Parameters of the object and hologram 
% [otf3d, psf3d, pupil3d] = OTF3D(Ny, Nx, Nz, lambda, deltaX, deltaY, deltaZ, offsetZ, sensor_size);
[otf3d, psf3d, pupil3d] = OTF3D_z(Ny, Nx, lambda, pps, z);


if any(strcmp(obj_name, {'geo', 'overlap', 'random', 'conhelix', 'cirhelix', 'SNUE'}))
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
        
        case 'inline'
            %{
              Inline hologram, I = |R|^2 + OR* + O*R + |O|^2, |R|^2 can be captured directly, 
              OR* + O*R + |O|^2 = iFT(FT(I) - FT(|R|^2 )), suppose |R|=1, real(OR* + O*R)=2real(O)
            %}
            holo = prop_field + conj(prop_field) + sum(obj3d.^2,3)/Nz; 

            holo = holo./max(abs(holo(:)));
%             holo = abs(holo)./max(abs(holo(:))).*exp(1i*angle(holo));
        case 'offline'
    end
    
    % add noise
    holo = awgn(holo, 40);  % holo = imnoise(holo, 'gaussian', 0, 0.001);
    
    out_filename = [outdir, obj_name, '_', holo_type , '_'];
else
% Experimental data
    issim = 0;
    temp = zeros(Ny, Nx, Nz);
    vol3d_vec = C2V(temp(:));  % Experiments have no reference object, just for coding convenience
   
    out_filename = [outdir, obj_name, '_'];
    
    holo = holoNorm(holo);
end

% [sigma, mu]= noiseAnalysis(holo)
% holo = holodenoise(holo); % Filter noise    

%% ========================== Reconstruction with back-propagation =================================

reobj_raw = iMatProp3D(holo, otf3d, pupil3d, holo_type);
temp = reobj_raw;
write3d(abs(temp), z*1e6, outdir, obj_name);
% write3d(abs(temp), z_scope*1e6, outdir, obj_name);
figure; imagesc(plotdatacube(abs(temp))); title('BackPropagation'); axis image; drawnow; colormap(hot); colorbar; axis off;

%% ============= Construct Minimization problem and Reconstruction with deconvolution ==============
A_fun = @(field3d_vec) VecProp3D(field3d_vec, otf3d, pupil3d, holo_type);
AT_fun = @(field2d_vec) iVecProp3D(field2d_vec, otf3d, pupil3d, holo_type);

Nyv = Ny;
Nxv = Nx*Nz*2;  %!!!!
Nzv = 1;

% Since the TV-based denoising term is not defined in a complex domain, the real and imaginary parts
% are independently processed in a vectorized form.
holo_vec = C2V(holo(:));

switch deconv_type
    case 'TwIST'  % works for both complex and inline holograms, but very slowly
%         tau = 0.01;   % This effects, need further investigation
%         tau_psi = 0.15;
        tolA = 1e-6;
        
        if strcmp(regu_type, 'L1')
            Phi_fun = @(vol3d, weight, epsi) L1phi(vol3d);
            
            [reobj_TwIST, ~, objective, times, ~, mses] = TwIST(holo_vec, A_fun, tau, ...
                'AT', AT_fun, ...
                'Phi', Phi_fun, ...
                'MaxIterA', iter_num, ...
                'ToleranceA', tolA, ...
                'TRUE_X', vol3d_vec);
            
        elseif strcmp(regu_type, 'TV')
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
        end
        
        reobj_TwIST = reshape(V2C(reobj_TwIST), Ny, Nx, Nz);
        reobj_deconv = reobj_TwIST;
        
    case 'GPSR' % works but not for overlapping
        tau = 0.02;  % regularization parameter
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
        
    case 'TVPD' % works but not for overlapping
        tau = 0.01;  % regularization parameter
        tolA = 1e-6; % stopping threshold
        
        otf3d_vec = C2V(otf3d(:));
        reobj_TVPD = tvpd(otf3d_vec, holo_vec, Nx, vol3d_vec, 2, 10,10, 20); %total variation regularization
%         [reobj_SALSA, reobj_debias, obj_gpsr, times, debias_s, mses] = GPSR_BB(holo_vec, A_fun, tau,...
%             'Debias', 1, ...
%             'AT', AT_fun,...
%             'Continuation', 1, ...
%             'ToleranceA', tolA, ...
%             'MaxIterA', iter_num, ...
%             'TRUE_X', vol3d_vec);
        
        reobj_TVPD = reshape(V2C(reobj_TVPD), Ny, Nx, Nz);
        reobj_deconv = reobj_TVPD;
     case 'TVAL3' % works but not for overlapping
        opts.mu = 2^5;
        opts.beta = 2^5;
        opts.tol = 1e-3;
        opts.maxit = iter_num;
        opts.TVnorm = 1;
        opts.nonneg = true;

        A_fun = @(vol3d) vecProp(vol3d, otf3d, pupil3d, holo_type);
        AT_fun = @(field2d) ivecProp(field2d, otf3d, pupil3d, holo_type);

        [reobj_TVAL, out]= TVAL3(A_fun, holo_vec,  Nyv*Nxv*2, 1, opts);
                
        reobj_TVAL = reshape(V2C(reobj_TVAL), Ny, Nx, Nz);
        reobj_deconv = reobj_TVAL;
        
    case 'NESTA' % works but not for overlapping
        delta = 0;
        mu = 1e-3;
        cg_tol = 1e-6; 
        cg_maxit = 40;
        CGwrapper(); % (first, zero-out the CGwrapper counters)

        A_fun = @(vol3d) vecProp(vol3d, otf3d, pupil3d, holo_type);
        AT_fun = @(field2d) ivecProp(field2d, otf3d, pupil3d, holo_type);
        
        opts = [];
        opts.Verbose = 10;
        opts.TolVar = 1e-4;
        opts.AAtinv = @(b) CGwrapper(A_fun, b, cg_tol, cg_maxit);
        opts.typemin = 'tv';

        [reobj_NESTA, niter, resid, mses] = NESTA(A_fun, AT_fun, holo_vec, mu, delta, opts);
                
        reobj_NESTA = reshape(V2C(reobj_NESTA), Ny, Nx, Nz);
        reobj_deconv = reobj_NESTA;
        
    case 'SALSA'  % not working
        A_fun = @(vol3d) vecProp(vol3d, otf3d, pupil3d, holo_type);
        AT_fun = @(field2d) ivecProp(field2d, otf3d, pupil3d, holo_type);
        
        %{
        reg parameter, extreme eigenvalues, rule of thumb:
        1e-4: severely ill-conditioned problems
        1e-1: mildly  ill-conditioned problems
        1:    when A = Unitary matrix
        %}
        tau = 1e-0;
        mu = tau/5;
        
        tolA = 1e-3; % Smoothing parameter (empirical setting)
        
        otf3dvec = C2V(otf3d(:));
        filter_FFT = 1./(abs(otf3dvec) + mu);
        
        invLS = @(x) real(iFT(filter_FFT.*FT(x)));  % ????
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
        A_fun = @(vol3d) vecProp(vol3d, otf3d, pupil3d, holo_type);
        AT_fun = @(field2d) ivecProp(field2d, otf3d, pupil3d, holo_type);
        
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

%% ======================================= Show the images =========================================
if isDebug
%     figure; semilogy(objective, 'LineWidth',2); ylabel('objective distance'); xlabel('Iteration');
%     print('-dpng', [out_filename, 'objdis.png']);
    
    figure; plot(mses, 'LineWidth',2); ylabel('MSE'); 
        
    % Back vecProp reconstruction
    figure; imagesc(plotdatacube(abs(reobj_raw))); title('BackPropagation'); axis image; drawnow; colormap(hot); colorbar; axis off;
    print('-dpng', [out_filename, 'BP', '.png']);
    
    % TwIST reconstruction
    temp = abs(reobj_deconv);
    temp = (temp-min(temp(:)))/(max(temp(:))-min(temp(:)));
    figure; imagesc(plotdatacube(abs(temp))); title('Reconstruction'); axis image; drawnow; colormap(hot); colorbar; axis off;
    print('-dpng', [out_filename, deconv_type, '_tau', num2str(tau), '_psi' , num2str(tau_psi), '_', num2str(iter_num), '.png']);

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
    view(0, 0); colorbar off; 
    saveas(gca, [out_filename, deconv_type, '_side.png'], 'png'); 
    view(0, 90); colorbar off; 
    saveas(gca, [out_filename, deconv_type, '_top.png'], 'png');
    
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