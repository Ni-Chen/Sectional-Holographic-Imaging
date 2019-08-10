function [im,otf3d,y]=GenerateData2D(noise,level)
%----------------------------------------------
% function [im,psh,y]=GenerateData(noise,level)
%
% Generate data for 2D deconvolution examples
%
% Input:   noise -> type of noise
%                    - 'Gaussian'  level is SNR
%                    - 'Poisson'   level is average number of photons per pixel
%
% Outputs: im    -> ground truth (star like object)
%          psf   -> psf
%          y     -> blurred and noisy data
%
%  Copyright (C) 2018 E. Soubies emmanuel.soubies@epfl.ch
%                     F. Soulez ferreol.soulez@univ-lyon1.fr
%----------------------------------------------

%% Ground truth
N=256;
im1=StarLikeSample(2,N,6,20,3,0.7);   % Star-like object (help StarLikeSample to see the details of parameters)
im2=StarLikeSample(2,N,4,10,3,0.7);   % Star-like object (help StarLikeSample to see the details of parameters)

offset=0.1*max(im1(:));
im = (im1 + offset).*exp(1i*1*pi*im2);   % Complex object field
%% PSF
% lamb=561;                % Illumination wavelength
% res=50;                  % Resolution (nm)
% Na=1.8;                  % Objective numerical aperture
% fc=2*Na/lamb*res;        % cut-off frequency
% ll=linspace(-0.5,0,sz(1)/2+1);
% lr=linspace(0,0.5,sz(1)/2);
% [X,Y]=meshgrid([ll,lr(2:end)],[ll,lr(2:end)]);
% [th,rho]=cart2pol(X,Y);
% OTF=fftshift(1/pi*(2*acos(abs(rho)/fc)-sin(2*acos(abs(rho)/fc))).*(rho<fc));
% psf=real(ifft2(OTF));

lambda = 532e-9;                % Illumination wavelength
pps = 10e-6;      % pixel pitch of CCD camera
dz = 30e-3;     % propagation distance

%% Data
% H=LinOpConv(otf3d);
H = LinOpWavePropKernel(lambda, N, N, pps, dz);
y_noNoise=H*im;
if strcmp(noise,'Gaussian')
    y = y_noNoise +  10^(-level/20).*random('Normal',zeros(size(y_noNoise)),ones(size(y_noNoise))) + ...
        1i* 10^(-level/20).*random('Normal',zeros(size(y_noNoise)),ones(size(y_noNoise)))  ;
%    y = y_noNoise;
elseif strcmp(noise,'Poisson')
    factor = level./mean(y_noNoise(:)) ;
    y_noNoise = y_noNoise.* factor;
    im = im.*factor;
    y = random('Poisson',y_noNoise);
else
    error('Wrong type of noise');
end
disp(['SNR data : ',num2str(20*log10(norm(y_noNoise(:))/norm(y_noNoise(:)-y(:))))]);

otf3d = H.mtf;
end
