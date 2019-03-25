function [y,dif]=MyTVnorm(x,weights)

if (nargin < 2), weights = 1;  end

TV=MyTV3D_conv(x,weights);

dif=sqrt(sum(abs(TV).^2,4));

y=sum(dif(:));

