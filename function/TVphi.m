%{
TV based image restoration
    xe = arg min_x 0.5 ||Ax-y||^2 + tau TV(x)
the solution is given the Chambolle algorithm:

- A. Chambolle, "An algorithm for total variation minimization and applications", Journal of 
Mathematical Imaging and Vision 20:89-97, 2004.

%}
function y = TVphi(x, Nvy, Nvx, Nvz)
	X = reshape(x, Nvy, Nvx, Nvz);
	[y, dif] = TVnorm(X);
end

function [y, dif] = TVnorm(x)
	TV = gradient3D(x);
	dif = sqrt(sum(TV.^2, 4));
	y = sum(dif(:));
end
