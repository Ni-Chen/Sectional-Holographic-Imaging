%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3d centered inverse Fourier transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Citation for this code and algorithm:
% Tatiana Latychevskaia and Hans-Werner Fink
% "Practical algorithms for simulation and reconstruction of digital in-line holograms",
% Appl. Optics 54, 2424 - 2434 (2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The code is written by Tatiana Latychevskaia, 2002
% The version of Matlab for this code is R2010b

function [out] = IFT3Dc(in)

[Ny, Nx, Nz] = size(in);

% f1 = zeros(Ny,Nx,Nz);
[x, y, z] = meshgrid((1:Nx), 1:Ny, 1:Nz);
% z = 1:Nz;
% for z = 1:Nz
%     f1(:,:,z) = exp(-1i*pi*(x + y + z - 3));
% end

f1 = exp(-1i*pi*(x + y + z - 3));

FT = ifftn(f1.*in);

out = f1.*FT;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%