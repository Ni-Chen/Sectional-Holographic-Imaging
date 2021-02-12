function y=MyTVphi(x,Nvx,Nvy,Nvz)

X=reshape(x,Nvx,Nvy,Nvz);

[y,dif]=MyTVnorm(X);
