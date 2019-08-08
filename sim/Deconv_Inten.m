close all; clear; clc;
addpath(genpath('../function/'));


data_dir = './output/';
if any(size(dir([data_dir, '*.mat']),1))
    movefile([data_dir, '*.mat'], [data_dir, 'mat/'])
    % delete([data_dir, '*.mat']);
end


%% ------------------------ Parameters --------------------------------------------
isGPU = 0;
maxit = 200;       % Max iterations

%% fix the random seed (for reproductibility)
rng(1);
useGPU(isGPU);
case_type = 1;

if ismember(case_type,[7, 8, 9])
    [im,psf,y] = GenerateData('Poisson',100);
else
    [im,psf,y] = GenerateData('Gaussian',20);
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

if ~isempty(reg_type)
    file_name_head = [cost_type, '+', reg_type];
else
    file_name_head = [cost_type];
end

if isNonNeg
    out_name = [file_name_head, '(NonNeg)'];
else
    out_name = file_name_head;
end
    

for n = 1:length(solv_types)
    solv_type = solv_types{n};

    if isNonNeg
        file_name = [file_name_head, '(NonNeg)+', solv_type];
    else
        file_name = [file_name_head, '+', solv_type];
    end    
    
    %% --------------------------------------------- Cost ----------------------------------------------
    switch cost_type
        case 'LS'   % Least-Squares
            LS = CostL2([],y);
            Fwd = LS*H;
            
        case 'KL'  % Kullback-Leibler divergence
            KL = CostKullLeib([],y,1e-6);     % Kullback-Leibler divergence data term
            Fwd = KL*H;
    end
    Fwd.doPrecomputation = 1;
     %         C = LinOpCpx(sz);
    
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
            optSolve = OptiChambPock(lamb*R_N12,G,Fwd);
            optSolve.OutOp=OutputOptiSNR(1,im,20);

            optSolve.tau = tau_cp;  
            
            if  strcmp(reg_type, 'TV') 
                disp('CP + TV');
                optSolve.sig = 1/(optSolve.tau*G.norm^2)*0.99;  % sig x tau x ?H?2 <= 1, 1/(optSolve.tau*G.norm^2)*0.99
            elseif strcmp(reg_type, 'HS') 
                disp('CP + HS');
                optSolve.sig = sig_cp;  % sig x tau x ?H?2 <= 1, 1/(optSolve.tau*G.norm^2)*0.99
            end
            
            file_name = [file_name, '(L', num2str(lamb), ',T', num2str(optSolve.tau) ,')'];
            
        case 'ADMM'
            %% ----------------------------------- ADMM LS + TV ----------------------------------------
            if ~isNonNeg % ADMM LS + TV/HS
                disp('ADMM + LS + TV/HS');
                Fn = {lamb*R_N12}; % Functionals F_n constituting the cost
                Hn = {G}; % Associated operators H_n
                rho_n = [rho_admm]; % Multipliers rho_n, [1e-1];

                optSolve = OptiADMM(Fwd,Fn,Hn,rho_n); % Declare optimizer
                
            else
                if  strcmp(cost_type, 'LS') % ADMM LS + TV + NonNeg
                    disp('ADMM + LS + NonNeg');
                    Fn = {lamb*R_N12, R_POS}; % Functionals F_n constituting the cost
                    Hn = {G, LinOpIdentity(size(im))}; % Associated operators H_n
                    rho_n = [rho_admm, rho_admm]; % Multipliers rho_n

                    optSolve = OptiADMM(Fwd,Fn,Hn,rho_n); % Declare optimizer
                    
                elseif strcmp(cost_type, 'KL') % ADMM KL + TV + NonNeg
                    disp('ADMM + KL + NonNeg');
                    Fn={KL,lamb*R_N12,R_POS};
                    Hn={H,G,LinOpDiag(sz)};
                    rho_n=[rho_admm,rho_admm,rho_admm];

                    optSolve=OptiADMM([],Fn,Hn,rho_n);
                end
            end
            optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10),[1 2]);
            
            file_name = [file_name, '(L', num2str(lamb), ',R', num2str(rho_n(1)) ,')'];
            
        case 'FISTA' % Forward-Backward Splitting optimization
            if strcmp(cost_type, 'LS')
                disp('FISTA + LS');
                optSolve= OptiFBS(Fwd,R_POS);
                
            elseif strcmp(cost_type, 'KL')
                disp('FISTA + KL');
                optSolve = OptiFBS(Fwd,R_POS);
                
                optSolve.gam = gam_fista;     % descent step
                optSolve.momRestart  = false; % true if the moment restart strategy is used
                file_name = [file_name, '(G', num2str(optSolve.gam),')'];
            end
            optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10));
            optSolve.fista = true;   % true if the accelerated version FISTA is used            
            
        case 'RL' % Richardson-Lucy algorithm
            if  strcmp(cost_type, 'LS') 
                disp('RL + LS');
                optSolve = OptiRichLucy(Fwd);
            elseif  strcmp(cost_type, 'KL')    
                if isempty(reg_type)
                    disp('RL + KL');
                    optSolve = OptiRichLucy(Fwd);
                else
                    disp('RL + KL +TV');
                    optSolve = OptiRichLucy(Fwd,1,lamb);
                    file_name = [file_name, '(L', num2str(lamb),')'];
                end                
            end
            optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10));
            
        case 'DR' % Douglas-Rachford
            disp('DR + LS + NonNeg');
            optSolve = OptiDouglasRachford(Fwd,R_POS,[],10,1.5);
            optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/10));
            
        case 'PD' % PrimalDual Condat KL
            if ~isNonNeg
               
            else  % PD + LS + TV + NonNeg
                if  strcmp(cost_type, 'LS') % ADMM LS + TV + NonNeg
                    disp('PD + LS + NonNeg');
                    Fn = {lamb*R_N12};
                    Hn = {G};
                    optSolve = OptiPrimalDualCondat(Fwd,R_POS,Fn,Hn);
                    optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/5),[1 3]);
                    optSolve.CvOp=TestCvgCombine(TestCvgCostRelative(1e-8,[1 3]), 'StepRelative', 1e-8);

                    optSolve.tau = tau_pd;          % set algorithm parameters
                    if  strcmp(reg_type, 'TV')
                        optSolve.sig = (1/optSolve.tau-Fwd.lip/2)/G.norm^2*0.9;   % sig x tau x ?H?2 <= 1, 1/(optSolve.tau*G.norm^2)*0.99
                    elseif strcmp(reg_type, 'HS')
                        optSolve.sig = sig_pd;  % sig x tau x ?H?2 <= 1, 1/(optSolve.tau*G.norm^2)*0.99
                    end
            
                    optSolve.rho = rho_pd;
                elseif strcmp(cost_type, 'KL') % ADMM LS + TV + NonNeg
                    disp('PD + KL + NonNeg');
                    Fn = {lamb*R_N12, KL};
                    Hn = {G,H};
                    optSolve = OptiPrimalDualCondat([],R_POS,Fn,Hn);
                    optSolve.OutOp=OutputOptiSNR(1,im,round(maxit/5),[2 3]);
                    optSolve.tau = tau_pd;          % set algorithm parameters
                    optSolve.sig = sig_pd;    %
                    optSolve.rho = rho_pd;
                end                
            end
            file_name = [file_name, '(L', num2str(lamb), ',T', num2str(optSolve.tau) ,',S', num2str(optSolve.sig) ,',R', num2str(optSolve.rho) ,')'];
             
        case 'VMLMB' % optSolve LS
            if ~isNonNeg  % VMLMB + NonNeg 
               disp('VMLMB');
                H.memoizeOpts.applyHtH=true;
                optSolve=OptiVMLMB(Fwd,[],[]);
%                 optSolve.OutOp=OutputOptiSNR(1,im,10);
                optSolve.m = 2;  % number of memorized step in hessian approximation (one step is enough for quadratic function)
            else
                if strcmp(reg_type, 'TV') && strcmp(cost_type, 'LS') % VMLMB LS + TV NonNeg 
                    disp('VMLMB + KL + TV + NonNeg');
                    hyperB = CostHyperBolic(G.sizeout, 1e-7, 3)*G;
                    C = Fwd + lamb*hyperB;
                    C.memoizeOpts.apply=true;
                    optSolve=OptiVMLMB(C,0.,[]);
%                     optSolve.OutOp=OutputOptiSNR(1,im,10);
                    optSolve.m=3;                                     % number of memorized step in hessian approximation
                    file_name = [file_name, '(L', num2str(lamb),')'];
                    
                elseif strcmp(reg_type, 'TV') && strcmp(cost_type, 'KL')  % VMLMB KL + TV NonNeg 
                    disp('VMLMB+ KL + TV + NonNeg');
                    H.memoizeOpts.applyHtH = true;
                    hyperB = CostHyperBolic(G.sizeout, 1e-7, 3)*G;
                    C = Fwd + lamb*hyperB; 
                    C.memoizeOpts.apply=true;
                    
                    optSolve=OptiVMLMB(C,0.,[]);                   
                    optSolve.m=3;                                     % number of memorized step in hessian approximation
                    file_name = [file_name, '(L', num2str(lamb),')'];
                else % VMLMB LS +  NonNeg 
                    disp('VMLMB + LS + NonNeg');
                    H.memoizeOpts.applyHtH=true;
                    optSolve=OptiVMLMB(Fwd,0.,[]);              
                    optSolve.m=3;                                     % number of memorized step in hessian approximation
                end                 
            end
            optSolve.OutOp=OutputOptiSNR(1,im,10);   
            
        case 'CG'  % ConjGrad LS
            disp('CG');
            A = H.makeHtH();
            b = H'*y;
            optSolve = OptiConjGrad(A,b);
            optSolve.OutOp=OutputOptiConjGrad(1,dot(y(:),y(:)),im,10);
            
        case 'GD'  %  Gradient Descent LS
            disp('GD');
            optSolve = OptiGradDsct(Fwd);
            optSolve.OutOp = OutputOptiSNR(1,im,round(maxit/10));  % for simulation 
            
    end
    
    optSolve.maxiter = maxit;                             % max number of iterations
%     optSolve.OutOp = OutputOptiSNR(1,im,round(maxit/10));  % for simulation
    optSolve.ItUpOut = 10;         % call OutputOpti update every ItUpOut iterations
%     optSolve.CvOp = TestCvgCombine(TestCvgCostRelative(1e-10), 'StepRelative', 1e-10);
%     optSolve.run(zeros(size(y)));             % run the algorithm
    optSolve.run(y);             % run the algorithm
    
    save([data_dir, file_name, '.mat'], 'optSolve');
    
    temp = abs((optSolve.xopt)); if isGPU; temp = gather(temp); end; imwrite(uint8(255*mat2gray(temp)), [data_dir,file_name, '.png']);
    imdisp(optSolve.xopt,file_name,1);
    
    
    if isGPU; reset(gpuDevice(1)); end
    solve_lst = dir(['./output/',  '*.mat']);  run('PlotCmp.m');
    
end
