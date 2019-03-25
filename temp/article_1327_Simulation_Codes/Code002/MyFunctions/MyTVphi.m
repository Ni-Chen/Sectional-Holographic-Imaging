function y=MyTVphi(x,Nvx,Nvy,Nvz,weights)

if (nargin < 5), weights = 1;  end

X=reshape(x,Nvx,Nvy,Nvz);

[y,dif]=MyTVnorm(X,weights);
