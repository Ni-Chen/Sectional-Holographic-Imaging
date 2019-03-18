function y = TVphi(x, Nvx, Nvy, Nvz)

	X = reshape(x, Nvx, Nvy, Nvz);

	[y, dif] = TVnorm(X);

end