close all; clear; clc;
%% Reading data
filename = 'CElegans-CY3';
imgDataPath = ['./', filename, '/'];
imgDataDir  = dir(imgDataPath);             % ??????
for i = 1:length(imgDataDir)
%     if(isequal(imgDataDir(i).name,'.')||... % ?????????????
%        isequal(imgDataDir(i).name,'..')||...
%        ~imgDataDir(i).isdir)                % ???????????
%            continue;
%     end
    imgDir = dir([imgDataPath imgDataDir(i).name '/*.tif']); 
    for j =1:length(imgDir)                 % ??????
        y(:,:,j) = imread([imgDataPath imgDataDir(i).name '/' imgDir(j).name]);
    end
end

imgDataPath = ['./', filename, '-PSF/'];
imgDataDir  = dir(imgDataPath);             % ??????
for i = 1:length(imgDataDir)
%     if(isequal(imgDataDir(i).name,'.')||... % ?????????????
%        isequal(imgDataDir(i).name,'..')||...
%        ~imgDataDir(i).isdir)                % ???????????
%            continue;
%     end
    imgDir = dir([imgDataPath imgDataDir(i).name '/*.tif']); 
    for j =1:length(imgDir)                 % ??????
        cpsf(:,:,j) = imread([imgDataPath imgDataDir(i).name '/' imgDir(j).name]);
    end
end

% cpsf=double(loadtiff(psfname));
cpsf = double(cpsf);
cpsf = cpsf/sum(cpsf(:));
% y=double(loadtiff(dataname));
y=double(y);
maxy=max(y(:));
y=y/maxy;
sz=size(y);

padsz=[0,0,52];
sznew=fft_best_dim(sz+2*padsz);
halfPad=(sznew-sz)/2;
cpsf=padarray(cpsf,halfPad,0,'both');


%% Least-squares data-fidelity term
H=LinOpConv(fftn(fftshift(cpsf)));                      % Convolution operator
H.memoizeOpts.apply=true;
S=LinOpSelectorPatch(sznew,halfPad+1,sznew-halfPad);   % Selector operator
L2=CostL2(S.sizeout,y);                                % L2 cost function
LS=L2*S;                                               % Least-squares data term
%% TV regularization
Freg=CostMixNorm21([sznew,3],4);      % TV regularizer: Mixed norm 2-1
Opreg=LinOpGrad(sznew);               % TV regularizer: Operator gradient
Opreg.memoizeOpts.apply=true;
lamb=2e-6;                            % Regularization parameter
%% Nonnegativity constraint
pos=CostNonNeg(sznew);                % Nonnegativity: Indicator function
Id=LinOpIdentity(sznew);              % Identity operator

dispIt=30;                     % Verbose every 30 iterations
maxIt=300;                     % Maximal number of iterations
Fn={LS,lamb*Freg,pos};         % Functionals F_n constituting the cost
Hn={H,Opreg,Id};               % Associated operators H_n
rho_n=[1e-3,1e-3,1e-3];        % Multipliers rho_n
ADMM=OptiADMM([],Fn,Hn,rho_n); % Declare optimizer
ADMM.OutOp=OutputOpti(1, '', round(maxIt/10),[1 2]); % build the output object
ADMM.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4), 'StepRelative',1e-4); % Set the convergence tests
ADMM.ItUpOut=dispIt;
ADMM.maxiter=maxIt;
ADMM.run(xopt);