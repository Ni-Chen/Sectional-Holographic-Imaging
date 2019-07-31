%-----------------------------------------------------------
% Deconv_LS_TV script: Deconvolution by minimizing the 
% Least-Squares function plus the TV regularizer:
%     0.5 ||Hx - y||^2  + i_{>0}(x) + lamb*TV(x)
% using 
%      - Chambolle-Pock
%      - ADMM 
%
% See LinOp, LinOpConv, LinOpGrad, Cost, CostL2,   
% CostMixNorm12, Opti, OptiChambPock, OptiADMM, OutpuOpti
%------------------------------------------------------------
close all;
% help Deconv_LS_TV
%--------------------------------------------------------------
%  Copyright (C) 2017 E. Soubies emmanuel.soubies@epfl.ch
%                     F. Soulez ferreol.soulez@univ-lyon1.fr
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%---------------------------------------------------------------
useGPU(1);

% -- fix the random seed (for reproductibility)
rng(1);

% -- Input image and psf
[im,otf,y] = GenerateData2D('Gaussian',30);
% y = y./max(y(:));

imdisp(abs(im),'Input Image (GT): mag', 1);
imdisp(angle(im),'Input Image (GT): phase', 1);

imdisp(abs(y),'Wave field: mag', 1);
imdisp(angle(y),'Wave field: phase', 1);

%% Conversion CPU/GPU is necessary
otf=gpuCpuConverter(otf);
im=gpuCpuConverter(im);
y=gpuCpuConverter(y);

sz = size(y);

% -- Convolution Operator definition
% H = LinOpConv(otf,false);
H = LinOpWaveProp(otf, false);

%% Back-propagation reconstruction
im_bp = LinOpAdjoint(H)*y;
imdisp(abs(im_bp),'BP Image (GT): mag', 1);
imdisp(angle(im_bp),'BP Image (GT): phase', 1);

% -- Functions definition
LS = CostL2([],y);            % Least-Squares data term
F = LS*H;
F.doPrecomputation = 1;
C = LinOpCpx(sz);

% -TV
G = LinOpGrad(C.sizeout,[1,2]);       % Operator Gradient
R_N12 = CostMixNorm21(G.sizeout,4);  % Mixed Norm 2-1
lamb = 2e-2;                  % Hyperparameter
% -Hessian Schatten norm
% G=LinOpHess(C.sizeout,'circular',[1,2]);
% R_N12=CostMixNormSchatt1(G.sizeout,1);
% lamb = 1e-2;                  % Hyperparameter

% -- Chambolle-Pock  LS + TV
maxit=200;
temp = G*C;
CP = OptiChambPock(lamb*R_N12,G*C,F);
CP.OutOp = OutputOptiSNR(1,im,round(maxit/10));
CP.CvOp = TestCvgCombine(TestCvgCostRelative(1e-10), 'StepRelative', 1e-10);  
CP.tau=0.1;               % algorithm parameters
CP.sig=0.02;   
CP.ItUpOut = round(maxit/10);                        % call OutputOpti update every ItUpOut iterations
CP.maxiter = maxit;                       % max number of iterations
CP.run(zeros(size(y)));               % run the algorithm 

% -- ADMM LS + TV
Fn = {lamb*R_N12}; % Functionals F_n constituting the cost
Hn = {G*C}; % Associated operators H_n
rho_n = [1e-1]; % Multipliers rho_n
% Here no solver needed in ADMM since the operator H'*H + alpha*G'*G is invertible
ADMM = OptiADMM(F,Fn,Hn,rho_n); % Declare optimizer
ADMM.CvOp = TestCvgCombine(TestCvgCostRelative(1e-10), 'StepRelative', 1e-10);
ADMM.OutOp = OutputOptiSNR(1,im,10);
ADMM.ItUpOut = 2;            % call OutputOpti update every ItUpOut iterations
ADMM.maxiter = 200;           % max number of iterations
ADMM.run(zeros(size(y)));   % run the algorithm 

% -- Display
temp = (CP.xopt);
imdisp(abs(temp),'LS + TV (CP):mag',1);
imdisp(angle(temp),'LS + TV (CP):phase',1);
temp = (ADMM.xopt);
imdisp(abs(temp),'LS + TV (ADMM):mag',1);
imdisp(angle(temp),'LS + TV (ADMM):phase',1);
figure;plot(CP.OutOp.iternum,CP.OutOp.evolcost,'LineWidth',1.5);grid; set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
hold all;plot(ADMM.OutOp.iternum,ADMM.OutOp.evolcost,'LineWidth',1.5); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
legend('CP','ADMM');title('Cost evolution');

figure;subplot(1,2,1); grid; hold all; title('Evolution SNR');set(gca,'FontSize',12);
semilogy(CP.OutOp.iternum,CP.OutOp.evolsnr,'LineWidth',1.5); 
semilogy(ADMM.OutOp.iternum,ADMM.OutOp.evolsnr,'LineWidth',1.5);
legend('LS+TV (CP)','LS+TV (ADMM)','Location','southeast');xlabel('Iterations');ylabel('SNR (dB)');
subplot(1,2,2);hold on; grid; title('Runing Time');set(gca,'FontSize',12);
orderCol = get(gca,'ColorOrder');
bar(1,[CP.time],'FaceColor',orderCol(1,:),'EdgeColor','k');
bar(2,[ADMM.time],'FaceColor',orderCol(2,:),'EdgeColor','k');
set(gca,'xtick',[1 2]);ylabel('Time (s)');
set(gca,'xticklabels',{'LS+TV (CP)','LS+TV (ADMM)'});set(gca,'XTickLabelRotation',45)  
