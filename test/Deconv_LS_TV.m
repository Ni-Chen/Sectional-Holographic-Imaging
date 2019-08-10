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
help Deconv_LS_TV
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
close all; clear; clc
addpath(genpath('../function/'));

% -- fix the random seed (for reproductibility)
rng(1);

N = 256;
lambda = 532e-9;   % Illumination wavelength
pps = 10e-6;       % pixel pitch of detector
dz = 30e-3;        % propagation distance

% -- Input image and psf
[im,otf,y] = GenerateData2D('Gaussian',30);
% otf = ifftshift(otf);
imdisp(abs(im),'Input Image (GT): mag', 1);
imdisp(angle(im),'Input Image (GT): phase', 1);

imdisp(abs(y),'Wave field: mag', 1);
imdisp(angle(y),'Wave field: phase', 1);

sz = size(y);

% -- Convolution Operator definition
H = LinOpConv(otf,false);


% -- Functions definition
LS = CostL2([],y);            % Least-Squares data term
F = LS*H;
F.doPrecomputation = 1;
C=LinOpCpx(sz);
% -TV
% Op = LinOpGrad(C.sizeout,[1,2]);       % Operator Gradient
% R = CostMixNorm21(Op.sizeout,4);  % Mixed Norm 2-1
% lamb = 2e-2;                  % Hyperparameter
% -Hessian Schatten norm
Op=LinOpHess(C.sizeout,'circular',[1,2]);
R=CostMixNormSchatt1(Op.sizeout,1);
lamb = 1e-2;                  % Hyperparameter

% -- Chambolle-Pock  LS + TV
maxit=100;
CP = OptiChambPock(lamb*R,Op*C,F);
CP.OutOp = OutputOptiSNR(1,im,round(maxit/10));
CP.CvOp = TestCvgCombine(TestCvgCostRelative(1e-15), 'StepRelative', 1e-15);  
CP.tau=1;               % algorithm parameters
CP.sig=0.02;   
CP.ItUpOut = round(maxit/10);                        % call OutputOpti update every ItUpOut iterations
CP.maxiter = maxit;                       % max number of iterations
CP.run(zeros(size(y)));               % run the algorithm 

% -- Display
imdisp(abs(CP.xopt),'LS + TV (CP):mag',1);
imdisp(angle(CP.xopt),'LS + TV (CP):phase',1);