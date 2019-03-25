function [x, numA, numAt, objective, distance, times, mses] = SALSA(y, A, tau, varargin)
%{
This function solves the regularization problem
    arg min_x = 0.5*||y - A x||_2^2 + tau phi_function(P^T x),
where A is a generic linear operator and phi_function(.) is a convex regularizarion function such that the 
solution of the denoising problem
    Psi_tau(y) = arg min_x = 0.5*||y - x||_2^2 + tau \phi(x),
is known.
----------------------------------------------------------------------------------------------------
For further details about the SALSA and C-SALSA algorithms, see the papers:
[1] M. Afonso, J. Bioucas-Dias, M. Figueiredo, "Fast image recovery using variable splitting and
constrained optimization," IEEE Transactions on Image Processing 19(9):2345-2356, 2010.
[2] M. Afonso, J. Bioucas-Dias, M. Figueiredo, "An Augmented Lagrangian based Method for the Constrained
Formulation of Imaging Inverse Problems", IEEE Transactions on Image Processing 20(3):681-695, 2011.
----------------------------------------------------------------------------------------------------
Authors: Manya Afonso, Jose Bioucas-Dias and Mario Figueiredo
Date: October, 2009.
Revised: August, 2011

For any technical queries/comments/bug reports, please contact: manya(dot)afonso(at)gmail(dot)com
----------------------------------------------------------------------------------------------------
Copyright (2009): Manya Afonso, Jose Bioucas-Dias, and Mario Figueiredo

SALSA is distributed under the terms of the GNU General Public License 2.0.

Permission to use, copy, modify, and distribute this software for any purpose without fee is hereby
granted, provided that this entire notice is included in all copies of any software which is or 
includes a copy or modification of this software and in all copies of the supporting documentation
for such software.

This software is being provided "as is", without any express or implied warranty. In particular, the
authors do not make any representation or warranty of any kind concerning the merchantability
of this software or its fitness for any particular purpose."
---------------------------------------------------------------------------------------------------

========================================= Required inputs ==========================================

 y: 1D vector or 2D array (image) of observations

 A: if y and x are both 1D vectors, A can be a k*n (where k is the size of y and n the size of x)
    matrix or a handle to a function that computes products of the form A*v, for some vector v.
    In any other case (if y and/or x are 2D arrays), A can be passed as a handle to a function which
    computes products of the form A*x; another handle to a function AT which computes products of the
    form A'*y is also required in this case. A can also be the convolution kernel, or mask.
    The size of x is determined as the size of the result of applying AT.

 tau: regularization parameter, usually a non-negative real parameter of the objective function(see above).
 
========================================== Optional inputs =========================================

 'MU' = parameter mu (weight for constraint function as a result of splitting); default = 1e-3;

 'Psi' = denoising function handle; handle to denoising function Default = soft threshold.

 'Phi' = function handle to regularizer needed to compute the objective function.
         Default = ||x||_1

 'TVINITIALIZATION' = must be set to 1 if using TV with initialization of the dual variables. Phi and 
   Psi are ignored, if specified.

 'TViters' = number of iterations of Chambolle's algorithm, if using TV with initialization of the
   dual variables. Default = 5. Ignored if 'TVINITIALIZATION' is 0 or unspecified.

 'AT'  = function handle for the function that implements the multiplication by the conjugate of A, 
         when A is a function handle. If A is an array, AT is ignored.

 'PT' = Transform in analysis prior formulation (inverse wavelet transform, in case of wavelet frames);
        default = identity

 'P' = inverse transform in analysis prior formulation (forward transform); default = identity

 'LS' = function handle for multiplication with the term (mu I + A^T A)^(-1), Can be omitted if A is a matrix.

 'STOPCRITERION' = type of stopping criterion to use
                   1 = stop when the relative change in the objective function falls below 'ToleranceA'
                   2 = stop when the relative norm of the difference between two consecutive estimates
                       falls below toleranceA
                   3 = stop when the objective function becomes equal or less than toleranceA.
                   Default = 1.

 'TOLERANCEA' = stopping threshold; Default = 0.001

 'MAXITERA' = maximum number of iterations allowed in the
              main phase of the algorithm.
              Default = 10000

 'INITIALIZATION' must be one of {0,1,2,array}
              0 -> Initialization at zero.
              1 -> Random initialization.
              2 -> initialization with A'*y.
              array -> initialization provided by the user.
              Default = 0;

 'TRUE_X' = if the true underlying x is passed in this argument, MSE evolution is computed

 'VERBOSE'  = work silently (0) or verbosely (1), default = 1

============================================== Outputs =============================================
  x_estimate = solution of the main algorithm

  numA, numAt = number of calls to A and At (if A is a matrix)
          NOTE: if A,AT are function handles, the wrapper callcounter()  must be used to increment 
          the global variable 'calls'
  
  objective = sequence of values of the objective function
              0.5*|| y - A x ||_2^2 + tau phi( x )

  distance = sequence of values of the constraint function

  times = CPU time after each iteration

  mses = sequence of MSE values, with respect to True_x, if it was given; if it was not given, mses 
         is empty, mses = [].
====================================================================================================
%}
    addpath('./utils/');  
    %-----------------------------------------------------------------------------------------------
    % test for number of required parametres
    %-----------------------------------------------------------------------------------------------
    if (nargin-length(varargin)) ~= 3
        error('Wrong number of required parameters');
    end

    %-----------------------------------------------------------------------------------------------
    % Set the defaults for the optional parameters
    %-----------------------------------------------------------------------------------------------
    stopCriterion = 3;
    compute_mse = 0;
    maxiter = 10000; % outer loop iterations
    init = 2;
    AT = 0;

    mu = 1e-3;

    psi_ok = 0;
    tolA = 0.001;
    isTVinitialization = 0;
    TViters = 5;
    verbose = 1;

    definedP = 0;
    definedPT = 0;

    isinvLS = 0;
    numA = 0;
    numAt = 0;

    %-----------------------------------------------------------------------------------------------
    % Read the optional parameters
    %-----------------------------------------------------------------------------------------------
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
                    psi_function = varargin{i+1};
                case 'PHI'
                    phi_function = varargin{i+1};
                case 'TVINITIALIZATION'
                    isTVinitialization = varargin{i+1};
                case 'TVITERS'
                    TViters = varargin{i+1};
                case 'MU'
                    mu = varargin{i+1};
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
                    AT = varargin{i+1};
                case 'VERBOSE'
                    verbose = varargin{i+1};
                case 'LS'
                    invLS = varargin{i+1};
                    isinvLS = 1;
                otherwise
                    % Hmmm, something wrong with the parameter string
                    error(['Unrecognized option: ''' varargin{i} '''']);
            end
        end
    end

    if (sum(stopCriterion == [1 2 3])==0)
        error('Unknown stopping criterion');
    end

    %%%% if P was given, make sure PT is also given, and that their dimensions are compatible.
    if xor( definedP,definedPT )
        error('If you give P you must also give PT, and vice versa.');
    end

    if (~definedP)
        P = @(x) x;
        PT = @(x) x;
    end

    % if A is a function handle, we have to check presence of AT,
    if isa(A, 'function_handle') && ~isa(AT,'function_handle')
        error('The function handle for transpose of A is missing');
    end

    % if A is a matrix, we find out dimensions of y and x, and create function handles for multiplication
    % by A and A', so that the code below doesn't have to distinguish between the handle/not-handle cases
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
            error('(A^T A + \mu I)^(-1) must be specified as a function handle.\n');
        else
            dummy = invLS(ATy);
            if ~((size(dummy,1)==size(ATy,1))&&(size(dummy,2)==size(ATy,2)))
                error('Specified function handle for solving the LS step does not seem compatible with the specified A and AT.\n')
            end
        end
    else %% A is a matrix
        [M,N] = size(matA);
        if ( M > N )
            ATA = transpose(matA)*matA;
            inverse_term = inv(ATA+mu*eye(N,N));
            %inverse_term = @(x) invATA*x;
        else
            AAT = matA*transpose(matA);
            inverse_term = AT*inv(mu*eye(M,M)+AAT)*A;
        end
        invLS = @(x) inverse_term*x;
    end

    % if psi was given, check to see if it is a handle and that it accepts two arguments
    if exist('psi_function','var') %& ~isTVinitialization
        if isTVinitialization
            fprintf('Warning: user specified Phi and Psi will not be used as TV with initialization flag has been set to 1.\n')
        else
            if isa(psi_function,'function_handle')
                try  % check if phi can be used, using Aty, which we know has same size as x
                    dummy = psi_function(PT(ATy), tau);
                    psi_ok = 1;
                catch
                    error('Something is wrong with function handle for psi');
                end
            else
                error('Psi1 does not seem to be a valid function handle');
            end
        end
    else %if nothing was given, use soft thresholding
        psi_function = @(x,tau) soft(x,tau);
    end
    
    % if psi exists, phi must also exist
    if (psi_ok == 1)
        if exist('phi_function','var')
            if isa(phi_function,'function_handle')
                try  % check if phi can be used, using Aty, which we know has
                    % same size as x
                    dummy = phi_function(PT(ATy));
                catch
                    error('Something is wrong with function handle for phi');
                end
            else
                error('Phi does not seem to be a valid function handle');
            end
        else
            error('If you give Psi you must also give Phi');
        end
    else  % if no psi and phi were given, simply use the l1 norm.
        if ~isTVinitialization
            phi_function = @(x) sum(abs(x(:)));
        else
            phi_function = @(x) TVnorm(x);
        end
    end

    %-----------------------------------------------------------------------------------------------
    % Initialization
    %-----------------------------------------------------------------------------------------------
    switch init
        case 0   % initialize at zero, using AT to find the size of x
            x = AT(zeros(size(y)));
        case 1   % initialize randomly, using AT to find the size of x
            x = randn(size(AT(zeros(size(y)))));
        case 2   % initialize x0 = A'*y
            x = ATy;
        case 33333
            % initial x was given as a function argument; just check size
            if size(A(x)) ~= size(y)
                error('Size of initial x is not compatible with A');
            end
        otherwise
            error('Unknown ''Initialization'' option');
    end

    % if the true x was given, check its size
%     if compute_mse && (size(true) ~= size(x))
    if compute_mse && (~isequal(size(true),size(x)))  % Modified by Ni Chen
        error('Initial x has incompatible size');
    end

    PTx = PT(x);

    % initializing
    u = PTx;
    bu = 0*u;
    threshold = tau/mu;

    criterion(1) = 1;

    % Compute and store initial value of the objective function
    resid = y-A(x);
    numA = numA + 1;
    prev_f = 0.5*(resid(:)'*resid(:)) + tau*phi_function(u);

    if verbose
        fprintf('Initial value of objective function = %3.3g\n',prev_f)
    end

    % start the clock
    t0 = cputime;

    times(1) = 0;
    objective(1) = prev_f;

    if compute_mse
        mses(1) = sum(sum((x-true).^2))/numel(x);
    end

%     if isTVinitialization
%         pux = 0*u;
%         puy = 0*u;
%         puz = 0*u;
%     end

    for outer = 1:maxiter
        xprev = x;
%         if isTVinitialization
%             [u, pux, puy, puz] = chambolle_prox_TV_stop(real(PTx-bu), 'lambda', threshold, 'tau', tau, 'maxiter', TViters, 'dualvars', [pux puy puz]);
            u = psi_function(real(PTx-bu), threshold);
%         else
%             u = psi_function(PTx-bu, tau);
%         end

        r = ATy + mu*P(u + bu);
        x = invLS(r);
        PTx = PT(x);
        bu = bu + (u - PTx);

        resid = y-A(x);
        numA = numA + 1;
        objective(outer+1) = 0.5*(resid(:)'*resid(:)) + tau*phi_function(u);

        if compute_mse
            err = x-true;
            mses(outer+1) = (err(:)'*err(:))/numel(x);
        end

        distance(outer) = norm(PTx(:)-u(:),2)/sqrt(norm(PTx(:),2)^2 + norm(u(:),2)^2);

        if (outer>1)
            % take no more than maxiter iterations
            switch stopCriterion
                case 1
                % compute the stopping criterion based on the relative variation of the objective function.
                    criterion(outer) = abs(objective(outer+1)-objective(outer))/objective(outer);
                case 2
                % compute the stopping criterion based on the relative variation of the estimate.
                    criterion(outer) = abs(norm(x(:)-xprev(:))/norm(x(:)));
                case 3
                % continue if not yet reached target value tolA
                    criterion(outer) = objective(outer+1);
                otherwise
                    error('Unknown stopping criterion');
            end

            if(criterion(outer) < tolA )
                if verbose
%                     fprintf('\niter = %d, obj = %3.3g, criterion = %3.3g, ( target = %3.3g )\t', outer, objective(outer+1), criterion(outer), tolA)
                    fprintf('\niter = %d, obj = %3.3g, criterion = %3.3g \t', outer, objective(outer+1), criterion(outer))
                    if compute_mse
                        fprintf('MSE = %g', mses(outer+1))
                    end

                    fprintf('\nConvergence reached.\n')
                end
                times(outer+1) = cputime - t0;
                break;
            end
        end

        if verbose
%             fprintf('\niter = %d, obj = %3.3g, stop criterion = %3.3g, ( target = %3.3g )\t', outer, objective(outer+1), criterion(outer), tolA)
            fprintf('\niter = %d, obj = %3.3g, criterion = %3.3g \t', outer, objective(outer+1), criterion(outer))
            if compute_mse
                fprintf('MSE = %g', mses(outer+1))
            end
        end        
       
        times(outer+1) = cputime - t0;
    end
    if verbose
        fprintf(1,'\nFinished the debiasing phase!\nResults:\n')
        fprintf(1,'||A x - y ||_2 = %10.3e\n',resid(:)'*resid(:))
        fprintf(1,'||x||_1 = %10.3e\n',sum(abs(x(:))))
        fprintf(1,'Objective function = %10.3e\n',objective(outer+1));
        fprintf(1,'CPU time so far = %10.3e\n', times(outer));
        fprintf(1,'\n');
    end
end

%--------------------------------------------------------------
% soft for both real and complex numbers
%--------------------------------------------------------------
function y = soft(x,T)
    %y = sign(x).*max(abs(x)-tau,0);
    y = max(abs(x) - T, 0);
    y = y./(y+T).*x;
end


