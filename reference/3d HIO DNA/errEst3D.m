%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D HYBRID INPUT OUTPUT ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Citation for this code/algorithm or any of its parts:
% Tatiana Latychevskaia and Hans-Werner Fink
% "Three-dimensional double helical DNA structure directly revealed 
% from its X-ray fiber diffraction pattern by iterative phase retrieval",
% Optics Express 26(23), 30991 - 31017 (2018) 
% https://doi.org/10.1364/OE.26.030991
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program calculates error between experimental and fitted diffraction
% pattern according to RMS error as explained in
% J. Miao et al "Phase retrieval from the magnitude of the Fourier 
% transforms of nonperiodic objects" 
% J. Opt. Soc. Am. A/ Vol. 15, No. 6/June 1998
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [iterative_error] = errEst3D(object, support)

% [Ny, Nx, Nz] = size(support);

a = 0;
b = 0;
% a = zeros(Ny,Nx,Nz); 
% b = zeros(Ny,Nx,Nz); 

temp = ((support> 0) & (object > 0));
a = a + temp.*object.^2;
b = b + (1-temp).*object.^2;

a = sum(a(:));
b = sum(b(:));

% for ii = 1:Nx
%     for jj = 1:Ny
%         for kk = 1:Nz
%             if ((support(ii,jj,kk)>0) && (object(ii,jj,kk))>0)
%                 a = a + object(ii,jj,kk)*object(ii,jj,kk);
%             else
%                 b = b + object(ii,jj,kk)*object(ii,jj,kk);
%             end
%         end
%     end
% end

iterative_error = sqrt(b./a);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%