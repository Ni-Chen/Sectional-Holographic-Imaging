function y=MyTVpsi(x,th,tau,iter,Nvx,Nvy,Nvz,weights)

if (nargin < 8), weights = 1;  end

X=reshape(x,Nvx,Nvy,Nvz);

Y=X-MyProjectionTV(X,tau,th*0.5,iter,weights);

y=reshape(Y,Nvx*Nvy*Nvz,1);