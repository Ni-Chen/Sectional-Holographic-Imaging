%% --------------------------------------------- Cost ----------------------------------------------
switch cost_type
    case 'LS'   % Least-Squares
        LS = CostL2([],y);
        Fwd = LS*H;
        
        Fwd.doPrecomputation = 1;
        C = LinOpCpx(sz);
        
    case 'KL'  % Kullback-Leibler divergence 
        KL = CostKullLeib([],y,1e-6);     % Kullback-Leibler divergence data term
        Fwd = KL*H;
        
        Fwd.doPrecomputation = 1;
        C = LinOpCpx(sz);
end

%% ----------------------------------------- regularizer -------------------------------------------
switch reg_type
    case 'TV'  % TV regularizer
        G = LinOpGrad(C.sizeout,[1,2]);       % Operator Gradient
        R_N12 = CostMixNorm21(G.sizeout,4);   % Mixed Norm 2-1
%         R_N12 = CostMixNorm21NonNeg(G.sizeout,4);   % Mixed Norm 2-1
        
    case '3DTV'  % 3D TV regularizer
        G = LinOpGrad(C.sizeout, [1,2,3]);    % TV regularizer: Operator gradient
        R_N12 = CostMixNorm21(G.sizeout,4);   % TV regularizer: Mixed norm 2-1, check
%         R_N12 = CostMixNorm21NonNeg(G.sizeout,4);   % Mixed Norm 2-1

    case 'HS'  % Hessian-Shatten
        G = LinOpHess(C.sizeout, 'circular', [1,2]);                 % Hessian Operator
        R_N12 = CostMixNormSchatt1(G.sizeout,1); % Mixed Norm 1-Schatten (p = 1)
end

R_POS = CostNonNeg(C.sizeout);           % Non-Negativity
Id = LinOpIdentity(C.sizeout);

% R_POS = CostNonNeg(sz);           % Non-Negativity
% Id = LinOpIdentity(sz);
        
%% ---------------------------------------- Optimization --------------------------------------------
switch solv_type   
    case 'CP'
        % ------------------------------- Chambolle-Pock  LS + TV --------------------------------- 
        optSolve = OptiChambPock(lamb*R_N12,G*C,Fwd);
        optSolve.tau = OptPara;  % 1, algorithm parameters, 15
%         optSolve.sig = 1/(optSolve.tau*G.norm^2)*0.99;  % sig x tau x ?H?2 <= 1
        optSolve.sig = 0.02; 
        
    case 'ADMM'
        %% ----------------------------------- ADMM LS + TV ----------------------------------------        
        if ~isNonNeg % ADMM LS + TV 
            Fn = {lamb*R_N12}; % Functionals F_n constituting the cost
            Hn = {G*C}; % Associated operators H_n
            rho_n = [OptPara]; % Multipliers rho_n, [1e-1];
            
        else  % ADMM LS + TV + NonNeg
            %             Fn = {lamb*R_N12, R_POS}; % Functionals F_n constituting the cost
            %             Hn = {G*C, Id*C}; % Associated operators H_n
            %             rho_n = [OptPara, OptPara]; % Multipliers rho_n
            if  strcmp(cost_type, 'LS') % ADMM LS + TV + NonNeg
                disp('ADMM + LS + NonNeg');
                Fn = {lamb*R_N12, R_POS}; % Functionals F_n constituting the cost
                Hn = {G*C, Id*C}; % Associated operators H_n
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
        
        optSolve = OptiADMM(Fwd,Fn,Hn,rho_n); % Declare optimizer
        
    case 'FISTA' % Forward-Backward Splitting optimization  
        optSolve = OptiFBS(Fwd,R_POS);
        optSolve.fista = true;   % true if the accelerated version FISTA is used
        optSolve.gam = 1;     % descent step
        optSolve.momRestart  = false; % true if the moment restart strategy is used
        
    case 'RL' % Richardson-Lucy algorithm
        optSolve = OptiRichLucy(Fwd, 1, lamb);
        optSolve.epsl = 1e-6; % smoothing parameter to make TV differentiable at 0 
        
    case 'PD' % PrimalDual Condat KL
%         lamb = 1e-3;                  % Hyperparameter
%         if ~isNonNeg
%             Fn = {lamb*R_N12, KL};
%             Hn = {Hess,H};
%             optSolve = OptiPrimalDualCondat([],R_POS,Fn,Hn);
%         else  % PD + LS + TV + NonNeg
            Fn = {lamb*R_N12};
            Hn = {G*C};
            optSolve = OptiPrimalDualCondat(Fwd,R_POS,Fn,Hn);
%         end
%         optSolve.OutOp = OutputOptiSNR(1, im, round(maxit/10), [2 3]);
        optSolve.tau = 1;          % set algorithm parameters
        optSolve.sig = (1/optSolve.tau-Fwd.lip/2)/G.norm^2*0.9;    %
        optSolve.rho = OptPara;          %
        
    case 'CG'  % ConjGrad LS 
        A = H;
        b = y;
        optSolve = OptiConjGrad(A,b*C);  
%         optSolve.OutOp = OutputOptiConjGrad(1,dot(y(:),y(:)),im,40);  
        
    case 'FCG'  % ConjGrad LS 
        optSolve = OptiFGP(A,b);  
end

optSolve.maxiter = maxit;                             % max number of iterations
% optSolve.OutOp = OutputOptiSNR(1,im,round(maxit/10));  % for simulation
optSolve.OutOp = OutputOpti(1,round(maxit/10));
% optSolve.OutOp = OutputOptiConjGrad(1,round(maxit/10));
optSolve.ItUpOut = round(maxit/10);         % call OutputOpti update every ItUpOut iterations
optSolve.CvOp = TestCvgCombine(TestCvgCostRelative(1e-8), 'StepRelative', 1e-8);
optSolve.run(zeros(size(otf)));             % run the algorithm


