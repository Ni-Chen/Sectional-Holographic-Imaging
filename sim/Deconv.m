close all; clear; clc;
addpath(genpath('./function/'));

% delete('./output/*');

%% ------------------------ Parameters --------------------------------------------
isGPU = 1;
noise_type = 'Gaussian'; % Poisson
maxit = 200;       % Max iterations

if strcmp(noise_type, 'Gaussian')
    [im,psf,y] = GenerateData('Gaussian',20);
elseif strcmp(noise_type, 'Poisson')
    [im,psf,y] = GenerateData('Poisson',100);
end
imwrite(uint8(255*mat2gray(abs(y))), ['./output/','conv_abs.png']);
imwrite(uint8(255*mat2gray(angle(y))), ['./output/','conv_phase.png']);

sz = size(y);

    
%% fix the random seed (for reproductibility)
rng(1);
useGPU(isGPU);

%% -- Convolution Operator definition
H=LinOpConv(fft2(psf));
LS=CostL2([],y);            % Least-Squares data term
F=LS*H;
F.doPrecomputation=1;
% H.memoizeOpts.applyHtH = true;
    
% Minimizing the Least-Squares function without any regularizer:
%     0.5 ||Hx - y||^2
isNonNeg = 0; cost_type = 'LS';  reg_type = '';  solv_types = {'GD', 'CG'};  % GD, VMLMB, CG

% Minimizing the Least-Squares function plus the NonNegativity constraint  without any regularizer:
%     0.5 ||Hx - y||^2  + i_{>0}(x)
% isNonNeg = 1; cost_type = 'LS';  reg_type = '';  solv_types = {'FISTA', 'DR'};  % FISTA, Douglas-Rachford and VMLMB

% % Minimizing the Least-Squares function plus the TV regularizer:
% %     0.5 ||Hx - y||^2  + lamb*TV(x)
% isNonNeg = 0; cost_type = 'LS';  reg_type = 'TV';  solv_types = {'CP', 'ADMM'};  % CP, ADMM
% lamb = 1e-3; rho_try = 0.1; tau_try = 15;

% Minimizing the Least-Squares function plus the NonNegativity constraint with TV regularizer:
%     0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*TV(x)
% isNonNeg = 1; cost_type = 'LS';  reg_type = 'TV';  solv_types = {'PD', 'ADMM'};  % PD, ADMM, VMLMB
% lamb = 1e-3; rho_try = 0.1; tau_try = 1;

% % Minimizing the Least-Squares function plus the Hessian-Schatten regularizer:
% %     0.5 ||Hx - y||^2  + lamb*||Hess*x||_{1,S_p}
% isNonNeg = 0; cost_type = 'LS';  reg_type = 'HS';  solv_types = {'CP', 'ADMM'};  % CP/ ADMM
% lamb = 5e-3; rho_try = 0.1; tau_try = 1; sig = 0.02;

% % Minimizing the Least-Squares function plus the NonNegativity constraint with Hessian-Schatten regularizer:
% %     0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*||Hess*x||_{S_p,1}
% isNonNeg = 1; cost_type = 'LS';  reg_type = 'HS';  solv_types = {'PD', 'ADMM'};  % PD
% lamb = 5e-3; rho_try = 0.1; tau_try = 1; sig = 0.01;


%%
% Minimizing the Kullback-Leibler divergence  plus the NonNegativity constraint without any
% regularizer: ???????
%     \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x)
% isNonNeg = 1; cost_type = 'KL';  reg_type = '';  solv_types = {'FISTA', 'RL'};  % GD, VMLMB, CG

%-----------------------------------------------------------
% Minimizing the Kullback-Leibler divergence plus the NonNegativity constraint with TV  regularizer:
%     \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x) + Tv(x)
% isNonNeg = 1; cost_type = 'KL';  reg_type = 'TV';  solv_types = {'PD', 'ADMM', 'RL'};  % PD, ADMM, RL, VMLMB
% lamb = 5e-3; rho_try = 0.1; tau_try = 100;

%-----------------------------------------------------------
% Minimizing the Kullback-Leibler divergence plus the NonNegativity constraint with Hessian-Schatten regularizer:
%    \sum_n -y_n log((H*x)_n + bet) + (H*x)_n + i_{>0}(x) + lamb*||Hess*x||_{1,S_p}
% isNonNeg = 1; cost_type = 'KL';  reg_type = 'HS';  solv_types = {'PD', 'ADMM'};  % PD, ADMM
% lamb = 5e-3; rho_try = 0.01; tau_try = 100;

for n = 1:length(solv_types)
    solv_type = solv_types{n};
    file_name = [];
    if ~isempty(reg_type)
        file_name = [cost_type, '+', reg_type];
    else
        file_name = [cost_type];
    end
    if isNonNeg
        file_name = [file_name, '(NonNeg)+', solv_type];
    else
        file_name = [file_name, '+', solv_type];
    end
    
    if any(strcmp(solv_type, {'CP', 'PD'}))
        OptPara = tau_try;
        file_name = [file_name, '(L', num2str(lamb), ',T', num2str(OptPara) ,')'];
    elseif any(strcmp(solv_type, {'ADMM'}))
        OptPara = rho_try;
        file_name = [file_name, '(L', num2str(lamb), ',T', num2str(OptPara) ,')'];
    else
        file_name = [file_name];
    end   
    
    %% --------------------------------------------- Cost ----------------------------------------------
    switch cost_type
        case 'LS'   % Least-Squares
            LS = CostL2([],y);
            F = LS*H;
            
            F.doPrecomputation = 1;
            %         C = LinOpCpx(sz);
            
        case 'KL'  % Kullback-Leibler divergence
            KL = CostKullLeib([],y,1e-6);     % Kullback-Leibler divergence data term
%             F = KL*H;
            
%             F.doPrecomputation = 1;
            %         C = LinOpCpx(sz);
    end
    
    %% ----------------------------------------- regularizer -------------------------------------------
    switch reg_type
        case 'TV'  % TV regularizer
            G = LinOpGrad(sz);       % Operator Gradient
            R_N12 = CostMixNorm21([sz,2],3);   % Mixed Norm 2-1
            
        case 'HS'  % Hessian-Shatten
            G = LinOpHess(sz);                 % Hessian Operator
            R_N12 = CostMixNormSchatt1([sz,3],1); % Mixed Norm 1-Schatten (p = 1)
            
        case ''
        otherwise
            
    end
    R_POS = CostNonNeg(sz);           % Non-Negativity
    Id = LinOpIdentity(sz);
    
    %% ---------------------------------------- Optimization --------------------------------------------
    switch solv_type
        case 'CP'
            % ------------------------------- Chambolle-Pock  LS + TV ---------------------------------
            optSolve = OptiChambPock(lamb*R_N12,G,F);
            optSolve.OutOp=OutputOptiSNR(1,im,20);

            optSolve.tau = OptPara;  % 1, algorithm parameters, 15
            optSolve.sig = 1/(optSolve.tau*G.norm^2)*0.99;  % sig x tau x ?H?2 <= 1
            
        case 'ADMM'
            %% ----------------------------------- ADMM LS + TV ----------------------------------------
            if ~isNonNeg % ADMM LS + TV
                Fn = {lamb*R_N12}; % Functionals F_n constituting the cost
                Hn = {G}; % Associated operators H_n
                rho_n = [OptPara]; % Multipliers rho_n, [1e-1];

                optSolve = OptiADMM(F,Fn,Hn,rho_n); % Declare optimizer
                optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10),[1 2]);
                
            else
                if  strcmp(cost_type, 'LS') % ADMM LS + TV + NonNeg
                    Fn = {lamb*R_N12, R_POS}; % Functionals F_n constituting the cost
                    Hn = {G, LinOpIdentity(size(im))}; % Associated operators H_n
                    rho_n = [OptPara, OptPara]; % Multipliers rho_n

                    optSolve = OptiADMM(F,Fn,Hn,rho_n); % Declare optimizer
                    optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10),[1 2]);
                    
                elseif strcmp(cost_type, 'KL') % ADMM LS + TV + NonNeg
                    Fn={KL,lamb*R_N12,R_POS};
                    Hn={H,G,LinOpDiag(sz)};
                    rho_n=[OptPara,OptPara,OptPara];

                    optSolve=OptiADMM([],Fn,Hn,rho_n);
                    optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10),[1 2]);
                end
            end
            
        case 'FISTA' % Forward-Backward Splitting optimization
            if strcmp(cost_type, 'KL')
                optSolve = OptiFBS(KL*H,R_POS);
                optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10));                
                
                optSolve.gam = 5;     % descent step
                optSolve.momRestart  = false; % true if the moment restart strategy is used
            elseif  strcmp(cost_type, 'LS')
                optSolve= OptiFBS(F,R_POS);
                optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10));
                optSolve.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4), 'StepRelative',1e-4);              

            end
            optSolve.fista = true;   % true if the accelerated version FISTA is used
            
        case 'RL' % Richardson-Lucy algorithm
            if  strcmp(cost_type, 'LS') % ADMM LS + TV + NonNeg
                optSolve = OptiRichLucy(F);
            elseif  strcmp(cost_type, 'KL')    % Richardson-Lucy KL + NonNeg (implicit)
                if isempty(reg_type)
                    optSolve = OptiRichLucy(KL*H);
                    optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10));

                else                    
                    optSolve = OptiRichLucy(KL*H,1,lamb);
                    optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10));

                end
            end
            
        case 'DR' % Richardson-Lucy algorithm

            optSolve = OptiDouglasRachford(F,R_POS,[],10,1.5);
            optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10));
            optSolve.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4), 'StepRelative',1e-4);
            
        case 'PD' % PrimalDual Condat KL
%             lamb = 1e-3;                  % Hyperparameter
            if ~isNonNeg
               
            else  % PD + LS + TV + NonNeg
                if  strcmp(cost_type, 'LS') % ADMM LS + TV + NonNeg
                    Fn = {lamb*R_N12};
                    Hn = {G};
                    optSolve = OptiPrimalDualCondat(F,R_POS,Fn,Hn);
                    optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10),[1 3]);
                    optSolve.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4,[1 3]), 'StepRelative',1e-4);

                    optSolve.tau = 1;          % set algorithm parameters
                    optSolve.sig = 1e-2;    %
                    optSolve.rho = 1.7;
                elseif strcmp(cost_type, 'KL') % ADMM LS + TV + NonNeg
                    Fn = {lamb*R_N12, KL};
                    Hn = {G,H};
                    optSolve = OptiPrimalDualCondat([],R_POS,Fn,Hn);
                    optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10),[2 3]);
                    optSolve.tau = 100;          % set algorithm parameters
                    optSolve.sig = 1e-2;    %
                    optSolve.rho = 1.2;
                end                
            end
             
        case 'VMLMB' % optSolve LS
            lamb = 1e-3;                  % Hyperparameter
            if ~isNonNeg
                hyperB = CostHyperBolic(G.sizeout, 1e-7, 3)*G*C;
                C1 = F + lamb*hyperB;
                C1.memoizeOpts.apply=true;
                optSolve = OptiVMLMB(C,0,[]);
                optSolve.m = 3;  % number of memorized step in hessian approximation (one step is enough for quadratic function)
            else
                hyperB = CostHyperBolic(G.sizeout, 1e-7, 3)*G;
                C1 = F + lamb*hyperB;
                C1.memoizeOpts.apply=true;
                optSolve = OptiVMLMB(C1,0.,[]);
                optSolve.m = 3;
            end
            
        case 'CG'  % ConjGrad LS
            A = H.makeHtH();
            b = H'*y;
            optSolve = OptiConjGrad(A,b);
            optSolve.OutOp=OutputOptiConjGrad(1,dot(y(:),y(:)),im,40);
            
        case 'GD'  %  Gradient Descent LS
            optSolve = OptiGradDsct(F);
            optSolve.OutOp = OutputOptiSNR(1,im,round(maxit/10));  % for simulation 
            
        case 'FCG'  % ConjGrad LS
            optSolve = OptiFGP(A,b);
    end
    
    optSolve.maxiter = maxit;                             % max number of iterations
%     optSolve.OutOp = OutputOptiSNR(1,im,round(maxit/10));  % for simulation
    optSolve.ItUpOut = round(maxit/10);         % call OutputOpti update every ItUpOut iterations
%     optSolve.CvOp = TestCvgCombine(TestCvgCostRelative(1e-10), 'StepRelative', 1e-10);
    optSolve.run(zeros(size(y)));             % run the algorithm
    
    save(['./output/', file_name, '.mat'], 'optSolve');
    
    temp = abs((optSolve.xopt)); if isGPU; temp = gather(temp); end; imwrite(uint8(255*mat2gray(temp)), ['./output/',file_name, '.png']);
    imdisp(optSolve.xopt,file_name,1);
    
    
    if isGPU; reset(gpuDevice(1)); end
    solve_lst = dir(['./output/',  '*.mat']);  run('PlotCmp.m');
    
end
