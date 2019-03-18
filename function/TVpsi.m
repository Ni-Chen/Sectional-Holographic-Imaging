function y = TVpsi(x, th, tau, iter, Nvx, Nvy, Nvz)

	X = reshape(x, Nvx, Nvy, Nvz);

	Y = X - ProjectionTV(X, tau, th*0.5, iter);

	y = reshape(Y, Nvx*Nvy*Nvz, 1);

end