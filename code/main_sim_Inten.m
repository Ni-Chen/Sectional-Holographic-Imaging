%{
----------------------------------------------------------------------------------------------------
Name: Gathering of the examples in GlobalBioIm

Author:   Ni Chen (chenni@snu.ac.kr)
Date:
Modified:

Reference:
- https://biomedical-imaging-group.github.io/GlobalBioIm/
----------------------------------------------------------------------------------------------------
%}


close all; clear; clc;
addpath(genpath('./function/'));

global isSim;

data_dir = './output/sim/real';
if any(size(dir([data_dir, '*.mat']),1))
%     movefile([data_dir, '*.mat'], [data_dir, 'real/'])
    % delete([data_dir, '*.mat']);
end

%% ------------------------ Parameters --------------------------------------------
isGPU = 0;
maxit = 200;       % Max iterations

%% fix the random seed (for reproductibility)
rng(1);
useGPU(isGPU);
case_type = 9;
isCpx = 0;
is3D = 0;

if ismember(case_type,[7, 8, 9])
    [im,psf,y] = setData('star','Poisson', 100, isCpx, is3D);
else
    [im,psf,y] = setData('star', 'Gaussian', 35, isCpx, is3D);
end


imwrite(uint8(255*mat2gray(abs(y))), ['./output/','conv_abs.png']);
imwrite(uint8(255*mat2gray(angle(y))), ['./output/','conv_phase.png']);

sz = size(y);
    
%% -- Convolution Operator definition
H = LinOpConv(fft2(psf));
  
switch case_type
    case 1  
        % LS:               0.5 ||Hx - y||^2
        isNonNeg = 0; cost_type = 'LS';  reg_type = '';    solv_types = {'GD', 'CG', 'VMLMB'};     % GD, VMLMB, CG
        
    case 2  
        % LS + NonNeg:      0.5 ||Hx - y||^2  + i_{>0}(x)
        isNonNeg = 1; cost_type = 'LS';  reg_type = '';    solv_types = {'FISTA', 'DR', 'VMLMB'};  % FISTA, Douglas-Rachford and VMLMB
        
    case 3  
        % LS + TV:          0.5 ||Hx - y||^2  + lamb*TV(x)
        isNonNeg = 0; cost_type = 'LS';  reg_type = 'TV';  solv_types = {'CP', 'ADMM'};  
        lamb = 1e-3; rho_admm = 0.1; tau_cp = 15; 
        
    case 4  
        % LS + TV + NonNeg: 0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*TV(x)
        isNonNeg = 1; cost_type = 'LS';  reg_type = 'TV';  solv_types = {'PD', 'ADMM', 'VMLMB'};   
        lamb = 1e-3; rho_admm = 0.1; tau_pd = 1;  rho_pd=1.95;  
        
    case 5  
        % LS + HS:          0.5 ||Hx - y||^2  + lamb*||Hess*x||_{1,S_p}
        isNonNeg = 0; cost_type = 'LS';  reg_type = 'HS';  solv_types = {'CP', 'ADMM'};   
        lamb = 5e-3; rho_admm = 0.1; tau_cp = 1; sig_cp = 0.02;
        
    case 6  
        % LS + HS + NonNeg: 0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*||Hess*x||_{S_p,1}
        isNonNeg = 1; cost_type = 'LS';  reg_type = 'HS';  solv_types = {'PD', 'ADMM'};   
        lamb = 5e-3; rho_admm = 0.1; tau_pd = 1; sig_pd = 0.01; rho_pd = 1.7;
        
    case 7  
        % KL + NonNeg:      \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x)
        isNonNeg = 1; cost_type = 'KL';  reg_type = '';    solv_types = {'FISTA', 'RL'};  
        gam_fista = 5;
        
    case 8  
        % KL + TV + NonNeg: \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x) + Tv(x)
        isNonNeg = 1; cost_type = 'KL';  reg_type = 'TV';  solv_types = {'PD', 'ADMM', 'RL', 'VMLMB'};  
        lamb = 5e-3; rho_admm = 1e-3; tau_pd = 100; sig_pd=1e-2;  rho_pd=1.2;
        
    case 9  
        % KL + HS + NonNeg: \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x) + lamb*||Hess*x||_{1,S_p}
        isNonNeg = 1; cost_type = 'KL';  reg_type = 'HS';  solv_types = {'PD', 'ADMM'}; 
        lamb = 5e-3; rho_admm = 0.01; tau_pd = 100; sig_pd=1e-2; rho_pd=1.2; 
end

method_name_head = strcat(cost_type, ~isempty(reg_type)*strcat('+', reg_type));
out_name = strcat(method_name_head, isNonNeg*('(NonNeg)'));
    
for n = 1:length(solv_types)
    solv_type = solv_types{n};

    method_name = strcat(method_name_head, isNonNeg*('(NonNeg)'), '+', solv_type); 
    
    %% --------------------------------------------- Cost ----------------------------------------------
    switch cost_type
        case 'LS'   % Least-Squares
            LS = CostL2([],y);
            F = LS*H;
            
        case 'KL'  % Kullback-Leibler divergence
            KL = CostKullLeib([],y,1e-6);     % Kullback-Leibler divergence data term
            F = KL*H;
    end
    F.doPrecomputation = 1;
    
    %% ----------------------------------------- regularizer -------------------------------------------
    switch reg_type
        case 'TV'  % TV regularizer
            G = LinOpGrad(sz);       % Operator Gradient
            R_N12 = CostMixNorm21([sz,2],3);   % Mixed Norm 2-1
            
        case 'HS'  % Hessian-Shatten
            G = LinOpHess(sz);                 % Hessian Operator
            R_N12 = CostMixNormSchatt1([sz,3],1); % Mixed Norm 1-Schatten (p = 1)            
    end
    if isNonNeg
        R_POS = CostNonNeg(sz);           % Non-Negativity
        Id = LinOpIdentity(sz);
    end    
    
    %% ---------------------------------------- Optimization --------------------------------------------
    switch solv_type
        case 'CP'
            optSolve = OptiChambPock(lamb*R_N12,G,F);
            optSolve.OutOp=OutputOpti(1,im,20);

            optSolve.tau = tau_cp;  
            
            if  strcmp(reg_type, 'TV') 
                disp('CP + TV');
                optSolve.sig = 1/(optSolve.tau*G.norm^2)*0.99;  % sig x tau x ?H?2 <= 1, 1/(optSolve.tau*G.norm^2)*0.99
            elseif strcmp(reg_type, 'HS') 
                disp('CP + HS');
                optSolve.sig = sig_cp;  % sig x tau x ?H?2 <= 1, 1/(optSolve.tau*G.norm^2)*0.99
            end
            
            method_name = [method_name, '(L', num2str(lamb), ',T', num2str(optSolve.tau) ,')'];
            
        case 'ADMM'
            %% ----------------------------------- ADMM LS + TV ----------------------------------------
            if ~isNonNeg % ADMM LS + TV/HS
                disp('ADMM + LS + TV/HS');
                Fn = {lamb*R_N12}; % Functionals F_n constituting the cost
                Hn = {G}; % Associated operators H_n
                rho_n = [rho_admm]; % Multipliers rho_n, [1e-1];

                optSolve = OptiADMM(F,Fn,Hn,rho_n); % Declare optimizer
                
            else
                if  strcmp(cost_type, 'LS') % ADMM LS + TV + NonNeg
                    disp('ADMM + LS + NonNeg');
                    Fn = {lamb*R_N12, R_POS}; % Functionals F_n constituting the cost
                    Hn = {G, LinOpIdentity(size(im))}; % Associated operators H_n
                    rho_n = [rho_admm, rho_admm]; % Multipliers rho_n

                    optSolve = OptiADMM(F,Fn,Hn,rho_n); % Declare optimizer
                    
                elseif strcmp(cost_type, 'KL') % ADMM KL + TV + NonNeg
                    disp('ADMM + KL + NonNeg');
                    Fn={KL,lamb*R_N12,R_POS};
                    Hn={H,G,LinOpDiag(sz)};
                    rho_n=[rho_admm,rho_admm,rho_admm];

                    optSolve=OptiADMM([],Fn,Hn,rho_n);
                end
            end
            optSolve.OutOp=OutputOpti(1,im,round(maxit/10),[1 2]);
            
            method_name = [method_name, '(L', num2str(lamb), ',R', num2str(rho_n(1)) ,')'];
            
        case 'FISTA' % Forward-Backward Splitting optimization
            if strcmp(cost_type, 'LS')
                disp('FISTA + LS');
                optSolve= OptiFBS(F,R_POS);
                
            elseif strcmp(cost_type, 'KL')
                disp('FISTA + KL');
                optSolve = OptiFBS(F,R_POS);
                
                optSolve.gam = gam_fista;     % descent step
                optSolve.momRestart  = false; % true if the moment restart strategy is used
                method_name = [method_name, '(G', num2str(optSolve.gam),')'];
            end
            optSolve.OutOp=OutputOpti(1,im,round(maxit/10));
            optSolve.fista = true;   % true if the accelerated version FISTA is used            
            
        case 'RL' % Richardson-Lucy algorithm
            if  strcmp(cost_type, 'LS') 
                disp('RL + LS');
                optSolve = OptiRichLucy(F);
            elseif  strcmp(cost_type, 'KL')    
                if isempty(reg_type)
                    disp('RL + KL');
                    optSolve = OptiRichLucy(F);
                else
                    disp('RL + KL +TV');
                    optSolve = OptiRichLucy(F,1,lamb);
                    method_name = [method_name, '(L', num2str(lamb),')'];
                end                
            end
            optSolve.OutOp=OutputOpti(1,im,round(maxit/10));
            
        case 'DR' % Douglas-Rachford
            disp('DR + LS + NonNeg');
            optSolve = OptiDouglasRachford(F,R_POS,[],10,1.5);
            optSolve.OutOp=OutputOpti(1,im,round(maxit/10));
            
        case 'PD' % PrimalDual Condat KL
            if ~isNonNeg
               
            else  % PD + LS + TV + NonNeg
                if  strcmp(cost_type, 'LS') % ADMM LS + TV + NonNeg
                    disp('PD + LS + NonNeg');
                    Fn = {lamb*R_N12};
                    Hn = {G};
                    optSolve = OptiPrimalDualCondat(F,R_POS,Fn,Hn);
                    optSolve.OutOp=OutputOpti(1,im,round(maxit/5),[1 3]);
                    optSolve.CvOp=TestCvgCombine(TestCvgCostRelative(1e-8,[1 3]), 'StepRelative', 1e-8);

                    optSolve.tau = tau_pd;          % set algorithm parameters
                    if  strcmp(reg_type, 'TV')
                        optSolve.sig = (1/optSolve.tau-F.lip/2)/G.norm^2*0.9;   % sig x tau x ?H?2 <= 1, 1/(optSolve.tau*G.norm^2)*0.99
                    elseif strcmp(reg_type, 'HS')
                        optSolve.sig = sig_pd;  % sig x tau x ?H?2 <= 1, 1/(optSolve.tau*G.norm^2)*0.99
                    end
            
                    optSolve.rho = rho_pd;
                elseif strcmp(cost_type, 'KL') % ADMM LS + TV + NonNeg
                    disp('PD + KL + NonNeg');
                    Fn = {lamb*R_N12, KL};
                    Hn = {G,H};
                    optSolve = OptiPrimalDualCondat([],R_POS,Fn,Hn);
                    optSolve.OutOp=OutputOpti(1,im,round(maxit/5),[2 3]);
                    optSolve.tau = tau_pd;          % set algorithm parameters
                    optSolve.sig = sig_pd;    %
                    optSolve.rho = rho_pd;
                end                
            end
            method_name = [method_name, '(L', num2str(lamb), ',T', num2str(optSolve.tau) ,',S', num2str(optSolve.sig) ,',R', num2str(optSolve.rho) ,')'];
             
        case 'VMLMB' % optSolve LS
            if ~isNonNeg  % VMLMB + NonNeg 
               disp('VMLMB');
                H.memoizeOpts.applyHtH=true;
                optSolve=OptiVMLMB(F,[],[]);
%                 optSolve.OutOp=OutputOpti(1,im,10);
                optSolve.m = 2;  % number of memorized step in hessian approximation (one step is enough for quadratic function)
            else
                if strcmp(reg_type, 'TV') && strcmp(cost_type, 'LS') % VMLMB LS + TV NonNeg 
                    disp('VMLMB + KL + TV + NonNeg');
                    hyperB = CostHyperBolic(G.sizeout, 1e-7, 3)*G;
                    C = F + lamb*hyperB;
                    C.memoizeOpts.apply=true;
                    optSolve=OptiVMLMB(C,0.,[]);
%                     optSolve.OutOp=OutputOpti(1,im,10);
                    optSolve.m=3;                                     % number of memorized step in hessian approximation
                    method_name = [method_name, '(L', num2str(lamb),')'];
                    
                elseif strcmp(reg_type, 'TV') && strcmp(cost_type, 'KL')  % VMLMB KL + TV NonNeg 
                    disp('VMLMB+ KL + TV + NonNeg');
                    H.memoizeOpts.applyHtH = true;
                    hyperB = CostHyperBolic(G.sizeout, 1e-7, 3)*G;
                    C = F + lamb*hyperB; 
                    C.memoizeOpts.apply=true;
                    
                    optSolve=OptiVMLMB(C,0.,[]);                   
                    optSolve.m=3;                                     % number of memorized step in hessian approximation
                    method_name = [method_name, '(L', num2str(lamb),')'];
                else % VMLMB LS +  NonNeg 
                    disp('VMLMB + LS + NonNeg');
                    H.memoizeOpts.applyHtH=true;
                    optSolve=OptiVMLMB(F,0.,[]);              
                    optSolve.m=3;                                     % number of memorized step in hessian approximation
                end                 
            end
            optSolve.OutOp=OutputOpti(1,im,10);   
            
        case 'CG'  % ConjGrad LS
            disp('CG');
            A = H.makeHtH();
            b = H'*y;
            optSolve = OptiConjGrad(A,b);
            optSolve.OutOp=OutputOptiConjGrad(1,dot(y(:),y(:)),im,10);
            
        case 'GD'  %  Gradient Descent LS
            disp('GD');
            optSolve = OptiGradDsct(F);
            optSolve.OutOp = OutputOpti(1,im,round(maxit/10));  % for simulation 
            
    end
    
    optSolve.maxiter = maxit;                             % max number of iterations
%     optSolve.OutOp = OutputOpti(1,im,round(maxit/10));  % for simulation
    optSolve.ItUpOut = 10;         % call OutputOpti update every ItUpOut iterations
%     optSolve.CvOp = TestCvgCombine(TestCvgCostRelative(1e-10), 'StepRelative', 1e-10);
%     optSolve.run(zeros(size(y)));             % run the algorithm
    optSolve.run(y);             % run the algorithm
    
    save([data_dir, method_name, '.mat'], 'optSolve');
    
    temp = abs((optSolve.xopt)); if isGPU; temp = gather(temp); end; imwrite(uint8(255*mat2gray(temp)), [data_dir,method_name, '.png']);
    imdisp(optSolve.xopt,method_name,1);
    
    

    
end

if isGPU; reset(gpuDevice(1)); end
solve_lst = dir([data_dir, method_name, '*.mat']);  run('PlotMult.m');
