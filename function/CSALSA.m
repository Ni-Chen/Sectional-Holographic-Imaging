function [x, numA, numAt, objective, distance1, distance2, criterion, times, mses] = CSALSA(y, A, mu1, mu2, sigma, varargin)

% Usage:
% [x_estimate, numA, numAt, objective, distance1, distance2, criterion, ...
% times, mses, ISNR] = csalsa(y, A, mu1, mu2, sigma, varargin)
%
% This function solves the constrained optimization problem 
%
%     arg min   phi( P^T x )
%          x
%       s.t.    || A x - y ||_2 <= epsilon
%
% where A is a generic linear operator, P^T is the analysis operator, and phi(.) is a regularizarion 
% function  such that the solution of the denoising problem 
%
%     Psi_tau(y) = arg min_x = 0.5*|| y - x ||_2^2 + tau \phi( x ), 
%
% is known, and the inverse of (alpha I + A^T A) can be computed easily.
% 
% -----------------------------------------------------------------------
%
% For further details about the SALSA and C-SALSA algorithms, see the papers: 
%
% [1] M. Afonso, J. Bioucas-Dias, and M. Figueiredo, "Fast image recovery
% using variable splitting and constrained optimization," IEEE Transactions 
% on Image Processing, vol. 19, no. 9, pp. 2345-2356, September, 2010.
%
% [2] M. Afonso, J. Bioucas-Dias, and M. Figueiredo, "An Augmented
% Lagrangian based Method for the Constrained Formulation of Imaging Inverse 
% Problems", IEEE Transactions on Image Processing, Vol. 20, no. 3, 
% pp 681 - 695, March, 2011.
%
% -----------------------------------------------------------------------
% Authors: Manya Afonso, Jose Bioucas-Dias and Mario Figueiredo
% Date: October, 2009.
% 
% For any technical queries/comments/bug reports, please contact:
% manya (dot) afonso (at) gmail (dot) com
% -----------------------------------------------------------------------

% Copyright (2009): Manya Afonso, Jose Bioucas-Dias, and Mario Figueiredo
% 
% SALSA is distributed under the terms of 
% the GNU General Public License 2.0.
% 
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
% 
%  ===== Required inputs =============
%
%  y: 1D vector or 2D array (image) of observations
%     
%  A: if y and x are both 1D vectors, A can be a 
%     k*n (where k is the size of y and n the size of x)
%     matrix or a handle to a function that computes
%     products of the form A*v, for some vector v.
%     In any other case (if y and/or x are 2D arrays), 
%     A can be passed as a handle to a function which computes 
%     products of the form A*x; another handle to a function 
%     AT which computes products of the form A'*y is also required 
%     in this case. A can also be the convolution kernel, or mask.
%     The size of x is determined as the size
%     of the result of applying AT.
%
%  mu1, mu2 : constraint weights. must be greater than zero.
%       min     Phi(u) + i_epsilon(v)
%       u,v,x
%           s.t.    u = x
%                   v = Ax - y
%           mu1 -> ||x-u||
%           mu2 -> ||Ax-y-v||
%
%   sigma : additive (Gaussian) noise variance
%
%  ===== Optional inputs =============
% 'CONTINUATIONFACTOR' = continuation factor on mu1, mu2.
%           must be greater than 1. default = 1 (no continuation)
%
% 'EPSILON' = epsilon parameter in the inequality ||Ax-y|| < epsilon .
%           default = sqrt(N+sqrt(8*N))*sigma; N = number of elements in y.
%  
%  'Psi' = denoising function handle; handle to denoising function
%          Default = soft threshold.
%
%  'Phi' = function handle to regularizer needed to compute the objective
%          function.
%          Default = ||x||_1
%
%  'AT'  = function handle for the function that implements
%          the multiplication by the conjugate of A, when A
%          is a function handle. 
%          If A is an array, AT is ignored.
%
% 'TVINITIALIZATION' = must be set to 1 if using TV with initialization of
% the dual variables. Phi and Psi are ignored, if specified.
%
% 'TViters' = number of iterations of Chambolle's algorithm, if using TV with initialization of
% the dual variables. Default = 5. Ignored if 'TVINITIALIZATION' is 0 or unspecified.
%
%  'AT'  = function handle for the function that implements
%          the multiplication by the conjugate of A, when A
%          is a function handle. 
%          If A is an array, AT is ignored.
%
% 'PT' = Transform in analysis prior formulation (inverse wavelet transform, in case of wavelet frames); 
%           default = identity
% 
% 'P' = inverse transform in analysis prior formulation (forward transform); default = identity
%
% 'LS' = function handle for multiplication with the term (mu I + A^T A)^(-1)
%           Can be omitted if A is a matrix.
%
%  'TOLERANCEA' = stopping threshold, typically < 1, but a positive integer in case 
%               stopping criterion is 3; Default = 0.001
% 
%  'MAXITERA' = maximum number of iterations allowed in the
%               main phase of the algorithm.
%               Default = 10000
%
%  'INITIALIZATION' must be one of {0,1,2,array}
%               0 -> Initialization at zero. 
%               1 -> Random initialization.
%               2 -> initialization with A'*y.
%               array -> initialization provided by the user.
%               Default = 0;
%
%  'TRUE_X' = if the true underlying x is passed in 
%                this argument, MSE evolution is computed
%
%  'VERBOSE'  = work silently (0) or verbosely (1), default = 1
%
%
% ===================================================  
% ============ Outputs ==============================
%   x = solution of the main algorithm
%
%   numA, numAt = number of calls to A and At
%           NOTE: if A,AT are function handles, the wrapper callcounter()
%           must be used to increment the global variable 'calls'
%
%   objective = sequence of values of the objective function phi( x )
%
%   distance1, distance2 = sequence of values of the constraint function
%
%   criterion = sequence of values of ||A x - y ||
%
%   times = CPU time after each iteration
%
%   mses = sequence of MSE values, with respect to True_x,
%          if it was given; if it was not given, mses is empty,
%          mses = [].
%
% ========================================================

%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
if (nargin-length(varargin)) ~= 5
     error('Wrong number of required parameters');
end

%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------
stopCriterion = 3;
compute_mse = 0;
maxiter = 10000; % outer loop iterations
init = 0;
AT = 0;


psi_ok = 0;
tolA = 0.001;
verbose = 1;

isTVinitialization = 0;
TViters = 5;


definedP = 0;
definedPT = 0;

isinvLS = 0;
numA = 0;
numAt = 0;

delta = 1;
epsilon = 0;

%--------------------------------------------------------------
% Read the optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    switch upper(varargin{i})
     case 'P'
         definedP = 1;
       P = varargin{i+1};   
     case 'PT'
         definedPT = 1;
       PT = varargin{i+1};  
     case 'PSI'
       psi = varargin{i+1};
     case 'PHI'
       phi = varargin{i+1};
     case 'TVINITIALIZATION'
        isTVinitialization = varargin{i+1};
     case 'TVITERS'
        TViters = varargin{i+1};
	 case 'STOPCRITERION'
       stopCriterion = varargin{i+1};
     case 'TOLERANCEA'       
       tolA = varargin{i+1};
     case 'MAXITERA'
       maxiter = varargin{i+1};
     case 'INITIALIZATION'
       if numel(varargin{i+1}) > 1   % we have an initial x
	 init = 33333;    % some flag to be used below
	 x = varargin{i+1};
       else 
	 init = varargin{i+1};
       end
     case 'TRUE_X'
       compute_mse = 1;
       true = varargin{i+1};
     case 'AT'
         isAT = 1;
       AT = varargin{i+1};
     case 'LS'
         invLS = varargin{i+1};
         isinvLS = 1;
     case 'VERBOSE'
         verbose = varargin{i+1};
     case 'CONTINUATIONFACTOR'
            delta = varargin{i+1};
     case 'EPSILON'
            epsilon = varargin{i+1};
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized option: ''' varargin{i} '''']);
    end
  end
end
%%%%%%%%%%%%%%

if (sum(stopCriterion == [1 2 3])==0)
   error(['Unknown stopping criterion']);
end

%%%% if P was given, make sure PT is also given, and that their
%%%% dimensions are compatible.
if xor( definedP,definedPT )
    error(['If you give P you must also give PT, and vice versa.']);
end

if (~definedP)
    P = @(x) x;
    PT = @(x) x;
end

% if A is a function handle, we have to check presence of AT,
if isa(A, 'function_handle') && ~isa(AT,'function_handle')
   error(['The function handle for transpose of A is missing']);
end 

% if A is a matrix, we find out dimensions of y and x,
% and create function handles for multiplication by A and A',
% so that the code below doesn't have to distinguish between
% the handle/not-handle cases
if ~( isa(A, 'function_handle') )
    
    if ~(size(y,1) == size(A,1) && size(y,2) == 1 )
        error('For a MxN observation matrix A, the measurement vector y must be a length M vector.');
    end
    
    if compute_mse && ~(size(true,1) == size(A,2) && size(true,2) == 1 )
        fprintf('\nWARNING: For a MxN observation matrix A, the vector x_true must be a length N vector.\n');
        compute_mse = 0;
    end
    matA = A;
    matAT = A';
   AT = @(x) matAT*x;
   A = @(x) matA*x;
end

% from this point down, A and AT are always function handles.

% Precompute A'*y since it'll be used a lot
ATy = AT(y);
numAt = numAt + 1;

% if A is a function handle, we have to check presence of invLS,
if isa(A, 'function_handle')
    
    if ( ~isinvLS || ~isa(invLS, 'function_handle'))
        error(['(A^T A + \mu I)^(-1) must be specified as a function handle.\n']);
    else
        dummy = invLS(ATy, mu1);
        if ~((size(dummy,1)==size(ATy,1))&&(size(dummy,2)==size(ATy,2)))
            error('Specified function handle for solving the LS step does not seem compatible with the specified A and AT.\n')
        end
    end
else %% A is a matrix
    [M,N] = size(matA);
    if ( M > N )
        ATA = transpose(matA)*matA;
        inverse_term = inv(ATA+(mu1+mu2)*eye(N,N));
        %inverse_term = @(x) invATA*x;
    else
        AAT = matA*transpose(matA);
        inverse_term = AT*inv((mu1+mu2)*eye(M,M)+AAT)*A;
    end
    invLS = @(x) inverse_term*x;        
end


% if psi was given, check to see if it is a handle and that it accepts two arguments
if exist('psi','var') 
    if isTVinitialization
        fprintf('Warning: user specified Phi and Psi will not be used as TV with initialization flag has been set to 1.\n')
    else
       if isa(psi,'function_handle')
          try  % check if phi can be used, using Aty, which we know has 
               % same size as x
               dummy = psi(PT(ATy),1/mu1); 
               psi_ok = 1;
          catch
             error(['Something is wrong with function handle for psi1'])
          end
       else
          error(['Psi1 does not seem to be a valid function handle']);
       end
    end
else %if nothing was given, use soft thresholding
   psi = @(x,tau) soft(x,tau);
end

% if psi exists, phi must also exist
if (psi_ok == 1)
   if exist('phi','var')
      if isa(phi,'function_handle')
         try  % check if phi can be used, using Aty, which we know has 
              % same size as x
              dummy = phi( PT(ATy) ); 
         catch
           error(['Something is wrong with function handle for phi'])
         end
      else
        error(['Phi does not seem to be a valid function handle']);
      end
   else
      error(['If you give Psi you must also give Phi']); 
   end
else  % if no psi and phi were given, simply use the l1 norm.
   if ~isTVinitialization
       phi = @(x) sum(abs(x(:))); 
   else
       phi = @(x) TVnorm(x);
   end
end

%--------------------------------------------------------------
% Initialization
%--------------------------------------------------------------
switch init
    case 0   % initialize at zero, using AT to find the size of x
       x = AT(zeros(size(y)));
    case 1   % initialize randomly, using AT to find the size of x
       x = randn(size(AT(zeros(size(y)))));
    case 2   % initialize x0 = A'*y
       x = ATy; 
       numAt = numAt + 1;
    case 33333
       % initial x was given as a function argument; just check size
       if size(A(x)) ~= size(y)
          error(['Size of initial x is not compatible with A']); 
       end
    otherwise
       error(['Unknown ''Initialization'' option']);
end

% if the true x was given, check its size
if compute_mse & (size(true) ~= size(x))  
   error(['Initial x has incompatible size']); 
end


% initializing
PTx = PT(x);

u = 0*PTx;
bu = 0*u;

v = 0*y;
bv = 0*y;

tau = mu1/mu2;

if (~epsilon)
    epsilon = sqrt(numel(y)+8*sqrt(numel(y)))*sigma;
end

Ax = A(x);
criterion(1) = norm(Ax(:)-y(:),2);

% Compute and store initial value of the objective function
resid =  y-Ax;
numA = numA + 1;

objective(1) = phi(x);

if verbose
    fprintf('Initial value of objective function = %3.3g\n',objective(1))
end

% start the clock
t0 = cputime;


times(1) = 0;

if compute_mse
   mses(1) = sum(sum((x-true).^2))/numel(x);
end

if isTVinitialization
    pux = 0*u;
    puy = 0*u;
end

Ax = A(x);
numA = numA + 1;
criterion(1) = norm(Ax(:)-y(:),2);
distance1(1) = norm(Ax(:)-y(:)-v(:),2);
distance2(1) = norm(PTx(:)-u(:),2);


if (verbose)
      fprintf(1,'iter = %d, restriction 1 = %g\t', 1, distance1(1) )
      fprintf(1,'restriction 2 = %g\t', distance2(1) )
      fprintf(1,'time = %g seconds\t||Ax-y|| = %g\t', times(1), criterion(1) )
      if compute_mse
          fprintf(1,'MSE = %g', mses(1) )
      end
      fprintf('\n')
end
    
for outer = 2:maxiter

    mu1inv = 1/mu1;
    mu2inv = 1/mu2;
    xprev = x;
    
    r = mu1*P(u+bu)+mu2*AT(y+v+bv);
    
    numAt = numAt + 1;
    
    x = invLS(r,mu1);
    
    PTx = PT(x);
    
%     if isTVinitialization
%         [u,pux,puy] = chambolle_prox_TV_stop(real(PTx-bu), 'lambda', mu1inv, 'maxiter', TViters, 'dualvars',[pux puy]);
%     else
        u = psi(real(PTx-bu), mu1inv);
%     end

    Ax = A(x);
    numA = numA + 1;
    ve = Ax-y-bv;
    n_ve = norm(ve(:));
    if n_ve <= epsilon;
        v = ve;
    else
        v =  ve/n_ve*epsilon;
    end
    
    bv = bv-(Ax-y-v);
    bu = bu-(PTx-u);
    
    criterion(outer) = norm(Ax(:)-y(:),2);
    distance1(outer) = norm(Ax(:)-y(:)-v(:),2);
    
    distance2(outer) = norm(PTx(:)-u(:),2);
    objective(outer) = phi(x);

    if compute_mse        
        mses(outer) = norm(x(:)-true(:),2)^2/numel(true);
    end
    
    times(outer) = cputime - t0;
     
    if (verbose)
      fprintf(1,'iter = %d, restriction 1 = %g\t', outer, distance1(outer) )
      fprintf(1,'restriction 2 = %g\t', distance2(outer) )
      fprintf(1,'objective 2 = %g\t', objective(outer) )
      fprintf(1,'time = %g seconds\t||Ax-y|| = %g\t', times(outer), criterion(outer) )
      if compute_mse
          fprintf(1,'MSE = %g', mses(outer) )
      end
      fprintf('\n')
    end

    mu1 = mu1*delta;
    mu2 = mu2*delta;
  

    if (outer>1)
        % take no less than miniter and no more than maxiter iterations
        switch stopCriterion
            case 1
                % compute the stopping criterion based on the relative
                % variation of the objective function.
                stopcriterion(outer) = abs(objective(outer)-objective(outer-1))/objective(outer);

                stop_flag = (stopcriterion(outer) < tolA)&&(criterion(outer)<= epsilon);

            case 2
                % compute the stopping criterion based on the relative
                % variation of the estimate.
                stopcriterion(outer) = abs(norm(x(:)-xprev(:))/norm(x(:)));
                stop_flag = (stopcriterion(outer) < tolA)&&(criterion(outer)<= epsilon);
            case 3
                % compute the stopping criterion based on the relative
                % variation of the estimate.
                stopcriterion(outer) = abs(criterion(outer)-criterion(outer-1))/criterion(outer);
                %stopcriterion(outer) = (criterion(outer-1)-criterion(outer))/criterion(outer-1);
                stop_flag = (stopcriterion(outer) < tolA)&&(criterion(outer)<= epsilon);
            case 4
                % minimum number of iterations

                stop_flag = (outer >= tolA)&&(criterion(outer)<= epsilon);
            otherwise
                error(['Unknown stopping criterion']);
        end

        if stop_flag
            if verbose
                fprintf('Stop criterion satisfied.\n');
            end
            break;
        end

    end
    
    times(outer) = cputime - t0;
    
end
end

%--------------------------------------------------------------
% soft for both real and complex numbers
%--------------------------------------------------------------
% function y = soft(x,T)
% if sum(abs(T(:)))==0
%    y = x;
% else
%    y = max(abs(x) - T, 0);
%    y = y./(y+T) .* x;
% end
% end

function y = soft(x,T)
    %y = sign(x).*max(abs(x)-tau,0);
    y = max(abs(x) - T, 0);
    y = y./(y+T).*x;
end


