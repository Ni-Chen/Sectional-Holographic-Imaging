function [im,otf,y]=GenerateDataCpx(noise, level)
%----------------------------------------------
%
% Input:   noise -> type of noise
%                    - 'Gaussian'  level is PSNR
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
%% Ground truth
N=128;
% im1=StarLikeSample(2,N,6,20,3,0.7);   % Star-like object (help StarLikeSample to see the details of parameters)
% im2=StarLikeSample(2,N,4,10,3,0.7);   % Star-like object (help StarLikeSample to see the details of parameters)

im1=StarLikeSample(2,N,4,20,3,0.8);   % Star-like object (help StarLikeSample to see the details of parameters)
im2=StarLikeSample(2,N,2,20,3,0.6);   % Star-like object (help StarLikeSample to see the details of parameters)


im1=1-mat2gray(im1);
im2=1-mat2gray(im2);

offset=0.1*max(im1(:));
% im = im1.*exp(1i*1*pi*im2);   % Complex object field
% im = (im1 + offset).*exp(1i*1*pi*im2);   % Complex object field
im = exp(log(0.5)*(im1)).*exp(-1i*0.99*pi*(im2));   % Complex object field

lambda = 532e-9;                % Illumination wavelength
pps = 3.45e-6;      % pixel pitch of CCD camera
dz = 50e-3;     % propagation distance

%% Data
% H=LinOpConv(otf3d);
H = LinOpWavePropKernel(lambda, N, N, pps, dz);
otf = H.mtf;
y_noNoise=H*im;

if strcmp(noise,'Gaussian')
    noise = 10^(-level/20).*random('Normal',zeros(size(y_noNoise)),ones(size(y_noNoise))) + ...
        1i* 10^(-level/20).*random('Normal',zeros(size(y_noNoise)),ones(size(y_noNoise)))  ;
    
    y = y_noNoise + noise;
    y = y./max(abs(y(:)));
    
elseif strcmp(noise,'Poisson')
    factor = level./mean((y_noNoise(:))) ;
    y_noNoise = (y_noNoise).* factor;
    
%     factor = level./mean(real(y_noNoise(:))) ;
%     y_noNoise_re = real(y_noNoise).* factor;
%     
%     factor = level./mean(imag(y_noNoise(:))) ;
%     y_noNoise_im = imag(y_noNoise).* factor;
    
    im = im.*factor;
    y = random('Poisson', real(y_noNoise)) + 1i*random('Poisson', imag(y_noNoise));
    
%     y = random('Poisson', y_noNoise_re) + 1i*imag(y_noNoise);
    
    y = y./max(abs(y(:)));
else
    error('Wrong type of noise');
end

disp(['SNR data : ',num2str(20*log10(norm(y_noNoise(:))/norm(y_noNoise(:)-y(:))))]);

end
