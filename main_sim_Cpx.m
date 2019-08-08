%{
----------------------------------------------------------------------------------------------------
Name: Complex deconvolution with GlobalBioIm

Author:   Ni Chen (chenni@snu.ac.kr)
Date:
Modified:

Reference:
- https://biomedical-imaging-group.github.io/GlobalBioIm/
----------------------------------------------------------------------------------------------------
%}

close all; clear; clc;
addpath(genpath('../function/'));

%% --------------------------------------- Initilization -------------------------------------------
global isSim;
isGPU = 0;
maxit = 100;       % Max iterations
case_type = 3;  %
isCpx = 1;
is3D = 0;

obj_name = 'star';
data_dir = './output/sim/';

if any(size(dir([data_dir, '*.mat']),1))
%     movefile([data_dir, '*.mat'], [data_dir, 'Cpx/'])
    delete([data_dir, '*.png']);
    delete([data_dir, '*.mat']);
end

%% -------------------------------------------------------------------------------------------------
rng(1);
useGPU(isGPU);
% if isGPU; reset(gpuDevice(1)); end

if ismember(case_type,[7, 8, 9])
    [im_ori,otf_ori,y_ori] = setData(obj_name, 'Poisson', 100, isCpx, is3D);
else
    [im_ori,otf_ori,y_ori] = setData(obj_name, 'Gaussian', 35, isCpx, is3D);
end

obj_name = strcat(obj_name, ~is3D*('2d'));

if ~is3D
    imwrite(uint8(255*mat2gray(abs(im_ori))), [data_dir,  obj_name,'_ori_abs.png']);
    imwrite(uint8(255*mat2gray(angle(im_ori))), [data_dir, obj_name, '_ori_phase.png']);
else
    figure('Name', 'Original'); show3d(im_ori, 0.01); axis normal; set(gcf,'Position',[0,0,500,500]); saveas(gca, [data_dir, obj_name, '_ori.png'], 'png');
end

imwrite(uint8(255*mat2gray(abs(y_ori))), [data_dir, obj_name, '_holo_abs.png']);
imwrite(uint8(255*mat2gray(angle(y_ori))), [data_dir,obj_name, '_holo_phase.png']);

if isGPU
    %     reset(gpuDevice());
    otf = gpuCpuConverter(otf_ori); im = gpuCpuConverter(im_ori); y = gpuCpuConverter(y_ori);
else
    otf = (otf_ori); im = (im_ori); y = (y_ori);
end

%% ----------------------------- Convolution Operator definition -----------------------------------
sz = size(otf);
Hp = LinOpConv(otf, 0, [1 2]);

if ~is3D
    H = Hp;
    
    % Back-propagation reconstruction
    im_bp = LinOpAdjoint(H)*y;
    if isGPU; temp = gather(im_bp); else  temp = im_bp;  end
    imwrite(uint8(255*mat2gray(abs(temp))), [data_dir, obj_name, '_BP_abs.png']);
    imwrite(uint8(255*mat2gray(angle(temp))), [data_dir, obj_name, '_BP_phase.png']);
else
    Hs = LinOpSum(sz,3);
    H = Hs*Hp;
    
    im_bp = LinOpAdjoint(H)*y;
    if isGPU; temp = gather(im_bp); else  temp = im_bp;  end
    figure('Name', 'BP'); show3d(abs(temp), 0.01); axis normal; set(gcf,'Position',[0,0,500,500]); saveas(gca, [data_dir, obj_name, '_BP.png'], 'png');
end

switch case_type
    case 1 % Not work
        % LS:               0.5 ||Hx - y||^2
        isNonNeg = 0; cost_type = 'LS';  reg_type = '';    solv_types = {'GD', 'CG'};     % GD, VMLMB, CG
        
    case 2 % work wrongly
        % LS + NonNeg:      0.5 ||Hx - y||^2  + i_{>0}(x)
        isNonNeg = 1; cost_type = 'LS';  reg_type = '';    solv_types = {'FISTA', 'DR'};  % 'VMLMB'
        gam_fista = 1;
        
    case 3 % gold
        % LS + TV:          0.5 ||Hx - y||^2  + lamb*TV(x)
        isNonNeg = 0; cost_type = 'LS';  reg_type = 'TV';  solv_types = {'ADMM', 'CP'};
        lamb = 1.5e-2; rho_admm = 0.1; tau_cp = 1;   tau_sig = 0.02;  % for 2D
%         lamb = 0.5e-3; rho_admm = 1e-1; tau_cp = 1;   % tau_sig = 0.02;  for 3D
        
    case 4 %
        % LS + TV + NonNeg: 0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*TV(x)
        isNonNeg = 1; cost_type = 'LS';  reg_type = 'TV';  solv_types = {'ADMM'}; %'VMLMB', PD
        lamb = 10e-3; rho_admm = 0.01; tau_pd = 0.1;  rho_pd = 1.1;
        
    case 5
        % LS + HS:          0.5 ||Hx - y||^2  + lamb*||Hess*x||_{1,S_p}
        isNonNeg = 0; cost_type = 'LS';  reg_type = 'HS';  solv_types = {'ADMM'}; % {'ADMM', 'CP'}
        lamb = 1e-2; rho_admm = 0.1; tau_cp = 0.2; sig_cp = 0.01;
        
    case 6
        % LS + HS + NonNeg: 0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*||Hess*x||_{S_p,1}
        isNonNeg = 1; cost_type = 'LS';  reg_type = 'HS';  solv_types = {'PD', 'ADMM'};
        lamb = 1e-2; rho_admm = 0.1; tau_pd = 1; sig_pd = 0.1; rho_pd = 1.7;
        
    case 7
        % KL + NonNeg:      \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x)
        isNonNeg = 1; cost_type = 'KL';  reg_type = '';    solv_types = {'FISTA', 'RL'};
        gam_fista = 5;
        
    case 8
        % KL + TV + NonNeg: \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x) + Tv(x)
        isNonNeg = 1; cost_type = 'KL';  reg_type = 'TV';  solv_types ={'PD', 'ADMM', 'RL'} ; %{'PD', 'ADMM', 'RL', 'VMLMB'};
        lamb = 1e-2; rho_admm = 1e-2; tau_pd = 100; sig_pd=1e-2;  rho_pd=1.2;
        
    case 9
        % KL + HS + NonNeg: \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x) + lamb*||Hess*x||_{1,S_p}
        isNonNeg = 1; cost_type = 'KL';  reg_type = 'HS';  solv_types = {'PD', 'ADMM'};
        lamb = 5e-3; rho_admm = 0.01; tau_pd = 100; sig_pd=1e-2; rho_pd=1.2;
end

method_name_head = strcat(obj_name, '_', cost_type, ~isempty(reg_type)*strcat('+', reg_type));
out_name = strcat(method_name_head, isNonNeg*('(NonNeg)'));

for n = 1:length(solv_types)
    solv_type = solv_types{n};
    method_name = strcat(method_name_head, isNonNeg*('(NonNeg)'), '+', solv_type);
    run('DeconvCpx.m');
    
    save([data_dir, method_name, '.mat'], 'optSolve');
    
    if isGPU; temp = gather(optSolve.xopt); else temp = optSolve.xopt; end
    if ~is3D
        abs(temp); imwrite(uint8(255*mat2gray(abs(temp))), [data_dir, method_name, '_abs.png']);
        imwrite(uint8(255*mat2gray(angle(temp))), [data_dir, method_name, '_phase.png']);
    else
        figure('Name', 'Deconv'); show3d(abs(temp), 0.01); axis normal; set(gcf,'Position',[0,0,500,500]); saveas(gca, [data_dir, obj_name, method_name, '.png'], 'png');
    end
    

    solve_lst = dir([data_dir, obj_name, '*.mat']);  run('PlotMult.m');    
    close all;
    
    if isGPU; reset(gpuDevice(1)); end
end
