function [y, dif] = TVnorm(x)

	TV = TV3D_conv(x);

	dif = sqrt(sum(TV.^2, 4));

	y = sum(dif(:));
end
