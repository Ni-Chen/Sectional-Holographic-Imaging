%   Solve the L1 penalized least-squares problem
%           min  mu|x|+.5||Ax-b||^2
%   using the solver FASTA.  
%
%  Inputs:
%    A   : A matrix or function handle
%    At  : The adjoint/transpose of A
%    b   : A column vector of measurements
%    mu  : Scalar regularization parameter
%    x0  : Initial guess of solution, often just a vector of zeros
%    opts: Optional inputs to FASTA
%
%   For this code to run, the solver "fasta.m" must be in your path.
%
%   For more details, see the FASTA user guide, or the paper "A field guide
%   to forward-backward splitting with a FASTA implementation."
%
%   Copyright: Tom Goldstein, 2014.



function [ solution, outs ] = fasta_sparseLeastSquares( A,At,b,mu,x0,opts )

%%  Check whether we have function handles or matrices
if ~isnumeric(A)
    assert(~isnumeric(At),'If A is a function handle, then At must be a handle as well.')
end
%  If we have matrices, create handles just to keep things uniform below
if isnumeric(A)
    At = @(x)A'*x;
    A = @(x) A*x;
end

%  Check for 'opts'  struct
if ~exist('opts','var') % if user didn't pass this arg, then create it
    opts = [];
end


%%  Define ingredients for FASTA
%  Note: fasta solves min f(Ax)+g(x).
%  f(z) = .5 ||z - b||^2
f    = @(z) .5*norm(z-b,'fro')^2;
grad = @(z) z-b;
% f    = @(z) .5*norm(b-z,'fro')^2;
% grad = @(z) b-z;
% g(z) = mu*|z|
% g = @(x) norm([real(x(:)); imag(x(:))],1)*mu;
g = @(x) norm(x(:),1)*mu;
% proxg(z,t) = argmin t*mu*|x|+.5||x-z||^2
prox = @(x,t) shrink(x,t*mu);

%% Call solver
[solution, outs] = norm_outputing_fasta(A,At,f,grad,g,prox,x0,opts);

end



%%  The vector shrink operator
function [ x ] = shrink( x,tau )
% re = real(x);
% im = imag(x);
% re = sign(re).*max(abs(re) - tau,0);
% im = sign(im).*max(abs(im) - tau,0);
% % x = sign(x).*max(abs(x) - tau,0);
% x = complex(re, im);

mag = abs(x);%x.*conj(x);
u = real(x);
v = imag(x);
u2 = u - tau*u./mag;
v2 = v - tau*v./mag;
% u2 = u - tau*abs(u)./mag;
% v2 = v - tau*abs(v)./mag;
idx = mag > tau;
x(idx) = complex(u2(idx), v2(idx));
x(~idx) = 0;
end

