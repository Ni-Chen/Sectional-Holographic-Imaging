%--------------------------------------------------------------
% This script performs 3D deconvolution using CP with
%    - Data-term: Least-Squares or Kullback-Leibler
%    - regul: TV or Hessian-Schatten norm
%--------------------------------------------------------------

close all; clear; clc;

addpath(genpath('./function/'));

%------------------------ Parameters --------------------------------------------
obj_name = 'star';  %'random', 'conhelix', 'circhelix', 'star'

isGPU = 0;
isNonNeg = 0;
cost_type = 'LS';  % LS, KL
reg_type = '3DTV';   % TV:, HS:Hessian-Shatten
solv_type = 'ADMM';  % CP, ADMM, CG, RL, FISTA, VMLMB

maxit = 1000;       % Max iterations

%% fix the random seed (for reproductibility)
rng(1);
useGPU(isGPU);

%% -------------------------------------- Image generation -----------------------------------------
[im_ori, otf_ori, y_ori] = setHoloData(obj_name, 'Gaussian', 50);

% CP: choose lambad and tau
lamb_try = [0.5e-3];
tau_try = [0.5 1];   % for CP
rho_try = [1e-1 1e-2];   % for ADMM

% lamb_try = 0.5e-3;
% tau_try = 1;

% lamb_try = 1e-3;
% tau_try = 1.5;

% lamb_try = 0.5e-3;
% tau_try = 0.5;

lamb_try =1e-3;
rho_try = 0.1;

n = 0;
if strcmp(solv_type, 'CP')
    OptPara_try = tau_try;
elseif strcmp(solv_type, 'ADMM')
    OptPara_try = rho_try;
end

for ilamb = 1:length(lamb_try)
    lamb = lamb_try(ilamb);
    
    for idxOptPara = 1:length(OptPara_try)
        OptPara = OptPara_try(idxOptPara);
        
        n = n + 1;
        disp(['Case ', num2str(n), ' of ', num2str(length(lamb_try)*length(OptPara_try))]);
       
        close all;
        if isGPU
            reset(gpuDevice(1));
            otf = gpuCpuConverter(otf_ori); im = gpuCpuConverter(im_ori); y = gpuCpuConverter(y_ori);
        else
             otf = (otf_ori); im = (im_ori); y = (y_ori);
        end
        
        %% ----------------------------------------- Foward Model ------------------------------------------
        sz = size(otf);
        H1 = LinOpConv(otf, 0, [1 2]);
        S = LinOpSum(sz,3);
        H = S*H1;
        if isNonNeg
            file_name = [obj_name, '_', cost_type, '+', reg_type, '(NonNeg)+', solv_type];
        else
            file_name = [obj_name, '_', cost_type, '+', reg_type, '+', solv_type];
        end
        file_name = [file_name, '(L', num2str(lamb), ',T', num2str(OptPara) ,')'];
        
        run('holoDeconv.m');
    end
end
if isGPU; reset(gpuDevice(1)); end
% solve_lst = dir(['./output/', obj_name, '*', reg_type, '*', solv_type, '*.mat']);
% solve_lst = dir(['./output/',obj_name, '*', solv_type, '*.mat']);
solve_lst = dir(['./output/',obj_name, '*.mat']);
run('dispPlotCmp.m');

