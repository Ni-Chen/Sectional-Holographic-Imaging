%% --------------------------------------------- Cost ----------------------------------------------
switch cost_type
    case 'LS'   % Least-Squares
        LS = CostL2([],y);
        F = LS*H;
        
        F.doPrecomputation = 1;
        C = LinOpCpx(sz);
        
    case 'KL'  % Kullback-Leibler divergence 
        KL = CostKullLeib([],y,1e-6);     % Kullback-Leibler divergence data term
        F = KL*H;
        
        F.doPrecomputation = 1;
        C = LinOpCpx(sz);
end

%% ----------------------------------------- regularizer -------------------------------------------
switch reg_type
    case 'TV'  % TV regularizer
        G = LinOpGrad(C.sizeout,[1,2]);       % Operator Gradient
        R_N12 = CostMixNorm21(G.sizeout,4);   % Mixed Norm 2-1
        
        R_POS = CostNonNeg(sz);           % Non-Negativity
        Id = LinOpIdentity(sz);
    case '3DTV'  % 3D TV regularizer
        G = LinOpGrad(C.sizeout, [1,2,3]);    % TV regularizer: Operator gradient
        R_N12 = CostMixNorm21(G.sizeout,4);   % TV regularizer: Mixed norm 2-1, check
        
        R_POS = CostNonNeg(sz);           % Non-Negativity
        Id = LinOpIdentity(sz);               % Identity Operator
        
    case 'HS'  % Hessian-Shatten
        G = LinOpHess(C.sizeout);                 % Hessian Operator
        R_N12 = CostMixNormSchatt1([sz, 3],1); % Mixed Norm 1-Schatten (p = 1)
        
        R_POS = CostNonNeg(sz);           % Non-Negativity
        Id = LinOpIdentity(sz);
end

%% ---------------------------------------- Optimization --------------------------------------------
switch solv_type   
    case 'CP'
        % ------------------------------- Chambolle-Pock  LS + TV ---------------------------------
%         lamb = 1e-3;  % circhelix: 1e-3       
        optSolve = OptiChambPock(lamb*R_N12,G*C,F);
        optSolve.tau = OptPara;  % 1, algorithm parameters, 15
        optSolve.sig = 1/(optSolve.tau*G.norm^2)*0.99;  % sig x tau x ?H?2 <= 1
%         optSolve.gam = 1;  % 1 or 2
%         optSolve.var = 1;  % 1: ; 2:
        
    case 'ADMM'
        %% ----------------------------------- ADMM LS + TV ----------------------------------------
%         lamb = 1e-3; % Hyperparameter, lamb = 1e-2;
        
        if ~isNonNeg % ADMM LS + TV 
            Fn = {lamb*R_N12}; % Functionals F_n constituting the cost
            Hn = {G*C}; % Associated operators H_n
            rho_n = [OptPara]; % Multipliers rho_n, [1e-1];
            
        else  % ADMM LS + TV + NonNeg 
            Fn = {lamb*R_N12, R_POS}; % Functionals F_n constituting the cost
            Hn = {G*C, Id}; % Associated operators H_n
            rho_n = [OptPara, OptPara]; % Multipliers rho_n           
        end
        
        % Here no solver needed in ADMM since the operator H'*H + alpha*G'*G is invertible
        optSolve = OptiADMM(F,Fn,Hn,rho_n); % Declare optimizer
        
    case 'FISTA' % Forward-Backward Splitting optimization  
        lamb = 1e-3;  % circhelix: 1e-3
        optSolve = OptiFBS(F,R_POS);
%         optSolve = OptiFBS(F, lamb*R_N12);
        optSolve.fista = true;   % true if the accelerated version FISTA is used
        optSolve.gam = 5;     % descent step
        optSolve.momRestart  = false; % true if the moment restart strategy is used
        
    case 'RL' % Richardson-Lucy algorithm
        lamb = 1e-2; % Hyperparameter for TV
        optSolve = OptiRichLucy(F, 1, lamb);
        optSolve.epsl = 1e-6; % smoothing parameter to make TV differentiable at 0 
        
    case 'PD' % PrimalDual Condat KL
        lamb = 1e-3;                  % Hyperparameter
%         if ~isNonNeg
%             Fn = {lamb*R_N12, KL};
%             Hn = {Hess,H};
%             optSolve = OptiPrimalDualCondat([],R_POS,Fn,Hn);
%         else  % PD + LS + TV + NonNeg
            Fn = {lamb*R_N12};
            Hn = {G*C};
            optSolve = OptiPrimalDualCondat(F,R_POS,Fn,Hn);
%         end
        optSolve.OutOp = OutputOptiSNR(1, im, round(maxit/10), [2 3]);
        optSolve.tau = 1;          % set algorithm parameters
        optSolve.sig = (1/optSolve.tau-F.lip/2)/G.norm^2*0.9;    %
        optSolve.rho = 1.95;          %
        
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
        optSolve.OutOp = OutputOptiConjGrad(1,dot(y(:),y(:)),im,40);  
        
    case 'FCG'  % ConjGrad LS 
        optSolve = OptiFGP(A,b);  
end

optSolve.maxiter = maxit;                             % max number of iterations
optSolve.OutOp = OutputOptiSNR(1,im,round(maxit/10));
optSolve.ItUpOut = round(maxit/10);         % call OutputOpti update every ItUpOut iterations
optSolve.CvOp = TestCvgCombine(TestCvgCostRelative(1e-8), 'StepRelative', 1e-8);
% optSolve.OutOpti = OutputOpti(true,round(maxit/10),costIndex)
optSolve.run(zeros(size(otf)));             % run the algorithm

save(['./output/', file_name, '.mat'], 'optSolve');
clear otf im y optSolve
