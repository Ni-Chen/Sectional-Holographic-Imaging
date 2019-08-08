
%% --------------------------------------------- Cost ----------------------------------------------
switch cost_type
    case 'LS'   % Least-Squares
        LS = CostL2([], y);
        F = LS*H;        
        
    case 'KL'  % Kullback-Leibler divergence
        KL = CostKullLeib([], y, 1e-6);     % Kullback-Leibler divergence data term
        F = KL*H;
end
F.doPrecomputation = 1;
Cpx = LinOpCpx(sz);

%% ----------------------------------------- regularizer -------------------------------------------
switch reg_type
    case 'TV'  % TV regularizer
        G = LinOpGrad(Cpx.sizeout, [1, 2]);       % Operator Gradient
        R_N12 = CostMixNorm21(G.sizeout, 4);   % Mixed Norm 2-1
        
        G_Cpx = G*Cpx;
        
    case 'HS'  % Hessian-Shatten
        %             G = LinOpHess(sz);                 % Hessian Operator
        G = LinOpHess(Cpx.sizeout);                 % Hessian Operator
        R_N12 = CostMixNormSchatt1([G.sizeout], 1); % Mixed Norm 1-Schatten (p = 1)
        %             R_N12 = CostMixNormSchatt1([sz, 3], 1); % Mixed Norm 1-Schatten (p = 1)
        
        G_Cpx = G*Cpx;
end
if isNonNeg
    R_POS = CostNonNeg(Cpx.sizeout);           % Non-Negativity, Not work for complex value
    Id = LinOpIdentity(Cpx.sizeout);
    Id_Cpx = Id*Cpx;
end

%% ---------------------------------------- Optimization --------------------------------------------
switch solv_type
    case 'CP'
        optSolve = OptiChambPock(lamb*R_N12, G_Cpx, F);
        optSolve.OutOp=OutputOptiSNR(1, im, round(maxit/10));
        
        optSolve.tau = tau_cp;
        
        if  strcmp(reg_type, 'TV')
            disp('CP + TV');
            optSolve.sig = 1/(optSolve.tau*G.norm^2)*0.99;  % sig x tau x ?H?2 <= 1, 1/(optSolve.tau*G.norm^2)*0.99
        elseif strcmp(reg_type, 'HS')
            disp('CP + HS');
            optSolve.sig = sig_cp;  % sig x tau x H2 <= 1, 1/(optSolve.tau*G.norm^2)*0.99
        end
        
        method_name = [method_name, '(L',  num2str(lamb), ',T', num2str(optSolve.tau) , ',S', num2str(optSolve.sig) ,')'];
        
    case 'ADMM'
        %% ----------------------------------- ADMM LS + TV ----------------------------------------
        if ~isNonNeg % ADMM LS + TV/HS
            if  strcmp(reg_type, 'TV')
                disp('ADMM + LS + TV');
                Fn = {lamb*R_N12}; % Functionals F_n constituting the cost
                Hn = {G_Cpx}; % Associated operators H_n
                rho_n = [rho_admm]; % Multipliers rho_n, [1e-1];
                
                optSolve = OptiADMM(F, Fn, Hn, rho_n); % Declare optimizer
            elseif strcmp(reg_type, 'HS')   % Not work
                disp('ADMM + LS + HS');
                Fn = {lamb*R_N12}; % Functionals F_n constituting the cost
                Hn = {G_Cpx}; % Associated operators H_n
                rho_n = [rho_admm]; % Multipliers rho_n, [1e-1];
                
                optSolve = OptiADMM(F, Fn, Hn, rho_n); % Declare optimizer
            end
        else
            if  strcmp(cost_type, 'LS') % ADMM LS + TV + NonNeg
                disp('ADMM + LS + NonNeg');
                Fn = {lamb*R_N12, R_POS}; % Functionals F_n constituting the cost
                Hn = {G_Cpx, Id_Cpx}; % Associated operators H_n
                rho_n = [rho_admm, rho_admm]; % Multipliers rho_n
                
                optSolve = OptiADMM(F, Fn, Hn, rho_n); % Declare optimizer
                
            elseif strcmp(cost_type, 'KL') % ADMM KL + TV + NonNeg
                disp('ADMM + KL + NonNeg');
                Fn={KL, lamb*R_N12, R_POS};
                Hn={H, G_Cpx, LinOpDiag(sz)};
                rho_n=[rho_admm, rho_admm, rho_admm];
                
                optSolve=OptiADMM([], Fn, Hn, rho_n);
            end
        end
        optSolve.OutOp=OutputOptiSNR(1, im, round(maxit/10), [1 2]);
        
        method_name = [method_name, '(L', num2str(lamb), ',R', num2str(rho_n(1)) ,')'];
        
    case 'FISTA' % Forward-Backward Splitting optimization
        if strcmp(cost_type, 'LS')
            disp('FISTA + LS');
            optSolve= OptiFBS(F,R_POS);
            
            optSolve.gam = gam_fista;     % descent step
            optSolve.momRestart  = false; % true if the moment restart strategy is used
            method_name = [method_name, '(G', num2str(optSolve.gam), ')'];
            
        elseif strcmp(cost_type, 'KL')
            disp('FISTA + KL');
            optSolve = OptiFBS(F,R_POS);
            
            optSolve.gam = gam_fista;     % descent step
            optSolve.momRestart  = false; % true if the moment restart strategy is used
            method_name = [method_name, '(G', num2str(optSolve.gam),')'];
        end
        optSolve.OutOp=OutputOptiSNR(1, im, round(maxit/10));
        optSolve.fista = true;   % true if the accelerated version FISTA is used
        
        
    case 'PD' % PrimalDual Condat KL
        if ~isNonNeg
            
        else  % PD + LS + TV + NonNeg
            if  strcmp(cost_type, 'LS') % ADMM LS + TV + NonNeg
                disp('PD + LS + NonNeg');
                Fn = {lamb*R_N12};
                Hn = {G_Cpx};
                optSolve = OptiPrimalDualCondat(F, R_POS, Fn, Hn);
                optSolve.OutOp=OutputOptiSNR(1, im, round(maxit/5), [1 3]);
                optSolve.CvOp=TestCvgCombine(TestCvgCostRelative(1e-8, [1 3]), 'StepRelative', 1e-8);
                
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
                Hn = {G*Cpx, H};
                optSolve = OptiPrimalDualCondat([], R_POS*Cpx, Fn, Hn);
                optSolve.OutOp=OutputOptiSNR(1, im, round(maxit/5), [2 3]);
                optSolve.tau = tau_pd;          % set algorithm parameters
                optSolve.sig = sig_pd;    %
                optSolve.rho = rho_pd;
            end
        end
        method_name = [method_name, '(L', num2str(lamb), ',T', num2str(optSolve.tau) , ',S', num2str(optSolve.sig) , ',R', num2str(optSolve.rho) , ')'];
        
    case 'VMLMB' % optSolve LS, x must be real
        if ~isNonNeg  % VMLMB + NonNeg
            disp('VMLMB');
            H.memoizeOpts.applyHtH=true;
            optSolve=OptiVMLMB(F,[], []);
            %                 optSolve.OutOp=OutputOptiSNR(1,im,10);
            optSolve.m = 2;  % number of memorized step in hessian approximation (one step is enough for quadratic function)
        else
            if strcmp(reg_type, 'TV') && strcmp(cost_type, 'LS') % VMLMB LS + TV NonNeg
                disp('VMLMB + KL + TV + NonNeg');
                hyperB = CostHyperBolic(G.sizeout, 1e-7, 3)*G;
                Cpx = F + lamb*hyperB;
                Cpx.memoizeOpts.apply=true;
                optSolve=OptiVMLMB(Cpx, 0., []);
                %                     optSolve.OutOp=OutputOptiSNR(1, im, 10);
                optSolve.m=3;                                     % number of memorized step in hessian approximation
                method_name = [method_name, '(L', num2str(lamb), ')'];
                
            elseif strcmp(reg_type, 'TV') && strcmp(cost_type, 'KL')  % VMLMB KL + TV NonNeg
                disp('VMLMB+ KL + TV + NonNeg');
                H.memoizeOpts.applyHtH = true;
                hyperB = CostHyperBolic(G.sizeout, 1e-7, 3)*G*Cpx;
                Cpx = F + lamb*hyperB;
                Cpx.memoizeOpts.apply=true;
                
                optSolve=OptiVMLMB(Cpx, 0., []);
                optSolve.m=3;                                     % number of memorized step in hessian approximation
                method_name = [method_name, '(L', num2str(lamb), ')'];
            else % VMLMB LS +  NonNeg
                disp('VMLMB + LS + NonNeg');
                H.memoizeOpts.applyHtH=true;
                optSolve=OptiVMLMB(F, 0., []);
                optSolve.m=3;                                     % number of memorized step in hessian approximation
            end
        end
        optSolve.OutOp=OutputOptiSNR(1, im, 10);
        
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
                optSolve = OptiRichLucy(F, 1, lamb);
                method_name = [method_name, '(L', num2str(lamb),')'];
            end
        end
        optSolve.OutOp=OutputOptiSNR(1, im, round(maxit/10));
        
    case 'DR' % Douglas-Rachford
        disp('DR + LS + NonNeg');
        gama = 1; % [0 Inf]
        lamb = 1;  % [0 2]
        optSolve = OptiDouglasRachford(F, R_POS, [], gama, lamb);
        optSolve.OutOp=OutputOptiSNR(1, im, round(maxit/10));
        
    case 'CG'  % ConjGrad LS
        disp('CG');
        A = H.makeHtH();
        b = H'*y;
        %             optSolve = OptiConjGrad(A, b*Cpx);
        optSolve = OptiConjGrad(A, b);
        optSolve.OutOp=OutputOptiConjGrad(1, dot(y(:), y(:)), im, 10);
        
    case 'GD'  %  Gradient Descent LS
        disp('GD');
        optSolve = OptiGradDsct(F);
        optSolve.OutOp = OutputOptiSNR(1, im, round(maxit/10));  % for simulation
end

optSolve.maxiter = maxit;                             % max number of iterations
%     optSolve.OutOp = OutputOptiSNR(1, im, round(maxit/10));  % for simulation
optSolve.ItUpOut = 5;         % call OutputOpti update every ItUpOut iterations
%     optSolve.CvOp = TestCvgCombine(TestCvgCostRelative(1e-10),  'StepRelative',  1e-10);

% run the algorithm
if ~is3D  % 2D
    optSolve.run(y);
    %     optSolve.run(zeros(size(y)));
else  % 3D
    optSolve.run(zeros(size(otf)));
end



