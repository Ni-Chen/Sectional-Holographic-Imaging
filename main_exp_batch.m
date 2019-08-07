%--------------------------------------------------------------
% This script performs 3D deconvolution using CP with
%    - Data-term: Least-Squares or Kullback-Leibler
%    - regul: TV or Hessian-Schatten norm
%--------------------------------------------------------------

close all; clear; clc;
addpath(genpath('./function/'));

%% ------------------------ Parameters --------------------------------------------
obj_name = 'fiber';  %'random', 'conhelix', 'circhelix', 'star'

isGPU = 1;
isNonNeg = 1;
cost_type = 'KL';  % LS, KL
reg_type = 'TV';   % TV:, HS:Hessian-Shatten
solv_type = 'ADMM';  % CP, ADMM, CG, RL, FISTA, VMLMB

maxit = 100;       % Max iterations

% CP: choose lambad and tau
% lamb_try = [0.5e-3 1e-3 5e-3];
% tau_try = [0.5 1 1.5];   % for CP
% rho_try = [1e-1 1e-2];   % for ADMM

lamb_try = 0.005;
rho_try = 0.1;
tau_try = 0.1;

%% fix the random seed (for reproductibility)
rng(1);
useGPU(isGPU);

%% -------------------------------------- Image loading -----------------------------------------
[im_ori, otf_ori, y_ori] = setHoloData(obj_name, 'Gaussian', 50);

n = 0;
if strcmp(solv_type, 'CP')
    OptPara_try = tau_try;
elseif strcmp(solv_type, 'ADMM')
    OptPara_try = rho_try;
else
    OptPara_try = tau_try;
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
        save(['./output/', file_name, '.mat'], 'optSolve');

    end
end

% % Back-propagation reconstruction
% im_bp = LinOpAdjoint(H)*y;
% temp = abs(im_bp);  %temp = abs(gather(optSolve.xopt));
% if isGPU
%     temp = gather(temp);
% end
% figure('Name', 'BP'); imagesc(plotdatacube(temp));  axis image; drawnow; colormap(hot);  axis off;
% print('-dpng', ['./output/', obj_name,'_BP.png']);

temp = abs((optSolve.xopt));  %temp = abs(gather(optSolve.xopt));
if isGPU
    temp = gather(temp);
end
figure('Name', 'Deconv'); imagesc(plotdatacube(temp));  axis image; drawnow; colormap(hot);  axis off;
print('-dpng', ['./output/', file_name,'_', num2str(maxit) ,'.png']);
clear otf im y optSolve

if isGPU; reset(gpuDevice(1)); end
% solve_lst = dir(['./output/', obj_name, '*', reg_type, '*', solv_type, '*.mat']);
% solve_lst = dir(['./output/', obj_name, '*', solv_type, '*.mat']);

solve_lst = dir(['./output/', obj_name, '*.mat']);
img_num = length(solve_lst);
if img_num > 0    
    legend_name = {};
    method_name = {};
    
    for imidx = 1:img_num 
        solve_name = solve_lst(imidx).name;
        load(['./output/', solve_name]);
        solve_result{imidx} = optSolve;
        
        temp = strrep(solve_name, [obj_name, '_'], '');
        method_name{imidx} = strrep(temp, '.mat', '');        
        legend_name = [legend_name, method_name{imidx}];        
    end
          
    figure('Name', 'Cost evolution'); 
    grid;title('Cost evolution'); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
    for imidx = 1:img_num 
        lineStyle = randLineStyle(); 
        plot(solve_result{imidx}.OutOp.iternum,solve_result{imidx}.OutOp.evolcost, lineStyle{1}, 'LineWidth', 1.5);     
        hold all;
    end
    
    legend(legend_name); 
    set(gcf,'paperpositionmode','auto');
    print('-dpng', ['./output/', obj_name, '_cost.png']);    
end

