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
% Experiments: dandelion, sh, beads, res
obj_name = 'SNUE';
holo_type = 'complex';  % complex; inline; offline;

% Output setting
isDebug = 1;

% Deconvolution setting
iter_num = 100;
regu_type = 'TV';  % 'TV', 'L1'
deconv_type = 'TwIST';  % 'TwIST','GPSR', TVAL3, SALSA, NESTA, TVPD

if(isDebug)
    outdir = './output/';  % Output files
else
    outdir = './Complex/figure/';  % Output images for paper writing
end

%% ======================================= Hologram ================================================
run([indir, obj_name, '_param.m']);  % Parameters of the object and hologram 
[otf3d, psf3d, pupil3d] = OTF3D(Ny, Nx, Nz, lambda, deltaX, deltaY, deltaZ, offsetZ, sensor_size);

if any(strcmp(obj_name, {'geo', 'overlap', 'random', 'conhelix', 'cirhelix', 'SNUE'}))
% Simulation data
    issim = 1;    
    load([indir, obj_name, '_3d.mat']);
    vol3d_vec = C2V(obj3d(:));  % For calculating MSE
    
    prop_field = MatProp3D(obj3d, otf3d, pupil3d);  % Field at the center plane of the 3D object
    
    %% Data
    H = LinOpVecProp3D(otf3d);
%     vol_field_ft = ifftshift(fft2(fftshift(obj3d)));
%     plane_field_ft = sum(otf3d.*vol_field_ft, 3);  % intergration along z axis
%     prop_field = ifftshift(ifft2(fftshift(plane_field_ft)));
    prop_field = H*obj3d;
    prop_field = sum(prop_field, 3);  % intergration along z axis
    
    
    switch holo_type
        case 'complex'
            holo = prop_field;
            tau = 0.02;   % This effects, need further investigation
            tau_psi = 0.05;
            
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
end

% holo = holodenoise(holo); % Filter noise    
% [sigma, mu]= noiseAnalysis(holo)

%% ========================== Reconstruction with back-propagation =================================
% reobj_raw = iMatProp3D(holo, otf3d, pupil3d, holo_type);
reobj_raw = conj(H)*holo;
write3d(abs(reobj_raw), z_scope*1e9, outdir, obj_name);

%% ============= Construct Minimization problem and Reconstruction with deconvolution ==============
A_fun = @(field3d_vec) VecProp3D(field3d_vec, otf3d, pupil3d, holo_type);
AT_fun = @(field2d_vec) iVecProp3D(field2d_vec, otf3d, pupil3d, holo_type);
sz = size(holo);
% H = LinOpConv(otf3d);   
H = LinOpConv('MTF', otf3d, false);
% H = LinOpVecProp3D(field3d_vec, otf3d, pupil3d, holo_type); 
L2 = CostL2(H.sizeout, holo);                            % L2 cost function
LS = L2*H;                                               % Least-Sqaures data term
pos = CostNonNeg(sz);                                    % Non-Negativity: Indicator function
Id = LinOpIdentity(sz);                                  % Identity Operator
KL = CostKullLeib(sz, holo, 1e-6);                       % Kullback-Leibler divergence data term

%% Parameters
lamb = 1e-4;        % Hyperparameter for initial deconvolution
maxIt = 30;         % Max iterations
Reg = 1;            % 1 for TV,  2 for Hessian-Schatten 
DataTerm = 1;       % 1 for LS,  2 for KL

if DataTerm==1
    dt = '0.5||Hx-y||_2^2'; 
else
    dt = 'KullLeib(Hx)';
end
if Reg==1
    rt = 'lamb*||x||_TV'; 
else
    rt = 'lamb*||x||_{Hess-Schat}';
end
fprintf('Minimized Cost function : \n');
fprintf(['\t F(x)  =  ', dt, ' + ', rt, ' + i_+(x) \n \n']);

%% Deconvolution
% -- Functions definition
if Reg==1
%     Freg = CostMixNorm21([sz, 3], 4);         % TV regularizer: Mixed Norm 2-1
    Freg = CostMixNorm21([sz, 2], 4);         % TV regularizer: Mixed Norm 2-1
    Opreg = LinOpGrad(sz);                   % TV regularizer: Operator Gradient
elseif Reg==2
    Freg = CostMixNormSchatt1([sz, 6], 1);     % Hessian-Shatten: Mixed Norm 1-Schatten (p=1)
    Opreg = LinOpHess(sz);                   % Hessian-Shatten: Hessian Operator
end

Opreg.useRFT = 0;

% -- ADMM
if DataTerm==1
    Fn = {lamb*Freg, pos};           % Functionals F_n
    Hn = {Opreg, Id};                    % Associated operators H_n
    rho_n = [1e-1, 1e-1];                % Multipliers rho_n
    ADMM = OptiADMM(LS, Fn, Hn, rho_n);   
else
    Fn = {KL, lamb*Freg, pos};          % Functionals F_n
    Hn = {H, Opreg, Id};                    % Associated operators H_n
    rho_n = [1e-3, 1e-3, 1e-3];             % Multipliers rho_n
    ADMM = OptiADMM([], Fn, Hn, rho_n);
end

ADMM.OutOp = OutputOpti(1, im, round(maxIt/10), [1 2]);
ADMM.CvOp = TestCvgCombine(TestCvgCostRelative(1e-4), 'StepRelative', 1e-4);  
ADMM.ItUpOut = round(maxIt/10);
ADMM.maxiter = maxIt;
ADMM.run(holo);

%% Displays
Orthoviews(ADMM.xopt, [], 'Deconvolved image');

figure; show3d(obj3d, 0.01); axis normal; title('Input Image');
figure; show3d(holo, 0.01); axis normal; title('Convolved and noisy data');
figure; show3d(ADMM.xopt, 0.01); axis normal; caxis([0 1]);title('Deconvolved image');

% -- Plot Evolution SNR,  cost  and Running Time for TV-Reg-Pos methods
figure;subplot(1, 3, 1); grid; hold all;
plot(ADMM.OutOp.iternum, ADMM.OutOp.evolcost, 'LineWidth', 1.5);  
set(gca, 'FontSize', 12);xlabel('Iterations');ylabel('Cost');
legend('ADMM');title('Cost evolution');
subplot(1, 3, 2); grid; hold all; title('Evolution SNR');set(gca, 'FontSize', 12);
semilogy(ADMM.OutOp.iternum, ADMM.OutOp.evolsnr, 'LineWidth', 1.5);
legend('ADMM', 'Location', 'southeast');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1, 3, 3);hold on; grid; title('Runing Time');set(gca, 'FontSize', 12);
orderCol = get(gca, 'ColorOrder');
bar(1, [ADMM.time], 'FaceColor', orderCol(1, :), 'EdgeColor', 'k');
set(gca, 'xtick', [1 2 3 4]);ylabel('Time (s)');
set(gca, 'xticklabels', {'ADMM'});set(gca, 'XTickLabelRotation', 45);




%% ======================================= Show the images =========================================
if isDebug
%     figure; semilogy(objective, 'LineWidth',2); ylabel('objective distance'); xlabel('Iteration');
%     print('-dpng', [out_filename, 'objdis.png']);
    
    figure; plot(times, mses, 'LineWidth',2); ylabel('MSE'); xlabel('CPU time (sec)');
        
    % Back vecProp reconstruction
    figure; imagesc(plotdatacube(abs(reobj_raw))); title('BackPropagation'); axis image; drawnow; colormap(hot); colorbar; axis off;
    print('-dpng', [out_filename, 'BP', '.png']);
    
    % TwIST reconstruction
    figure; imagesc(plotdatacube(abs(reobj_deconv))); title('Reconstruction'); axis image; drawnow; colormap(hot); colorbar; axis off;
    print('-dpng', [out_filename, deconv_type, '_', num2str(iter_num), '.png']);

    figure; show3d(abs(reobj_deconv), 0.01); axis normal; set(gcf,'Position',[0,0,500,500]);
    print('-dpng', [out_filename, deconv_type,  '_', num2str(iter_num), '_3d.png']);
    view(0, 0); colorbar off; print('-dpng', [out_filename, deconv_type,  '_', num2str(iter_num), '_side.png']);
    view(0, 90); colorbar off; print('-dpng', [out_filename, deconv_type, '_', num2str(iter_num), '_top.png']);
    
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