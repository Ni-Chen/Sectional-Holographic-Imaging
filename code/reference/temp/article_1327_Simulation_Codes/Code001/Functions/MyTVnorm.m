function [y,dif]=MyTVnorm(x)

TV=MyTV3D_conv(x);

dif=sqrt(sum(TV.^2,4));

y=sum(dif(:));

