function [im,otf,y]= setData(obj_name, noise, level, isCpx, is3D)
%----------------------------------------------
% function [im,psh,y]=GenerateData(noise,level)
%
% Generate data for 2D deconvolution examples
%
% Input:   noise -> type of noise
%                    - 'Gaussian'  level is PSNR
%                    - 'Poisson'   level is average number of photons per pixel
%
% Outputs: im    -> ground truth (star like object)
%          otf   -> otf
%          y     -> blurred and noisy data
%
%  Copyright (C) 2018 E. Soubies emmanuel.soubies@epfl.ch
%                     F. Soulez ferreol.soulez@univ-lyon1.fr
%----------------------------------------------
    global isSim;
    if ~any(strcmp(obj_name, {'random', 'conhelix', 'circhelix', 'star'}))
        % Experimental data
        
        isSim = 0;        
        
        run(['./data/', obj_name, '_param.m']);  % Parameters of the object and hologram
        y = holo;

        H = LinOpWavePropKernel(lambda, Ny, Nx, pps, z);
        otf = H.mtf;

        im = [];
        
        
    else
        % Simulation data
        isSim = 1;
        
        Nxy= 128;
        Nz = Nxy;

        if isCpx
            lambda = 532e-9;                % Illumination wavelength
            pps = 3.45e-6;      % pixel pitch of CCD camera
            z0 = 50e-3;     % propagation distance

            if is3D
                switch obj_name
                    case 'random'
                        im = randomScatter(Nxy, Nz, 1);
                        dz = 500e-6;
                    case 'conhelix'
                        im = conHelix(Nxy, Nz, 0.8, 8);
                        dz = 500e-6;
                    case 'circhelix'
                        im = circHelix(Nxy, Nz, 0.4, 6);
                        dz = 500e-6;
                    case 'star' % Star-like object (help StarLikeSample to see the details of parameters)
                        im = StarLikeSample(3, Nxy, 2, 20, 20, 0.7);
                        dz = 1000e-6;
                    otherwise  % experimetnal data
                end
                z = z0 + ((1:Nz)-Nz/2)*dz;
            else
                z = z0;

                im1=StarLikeSample(2,Nxy,4,20,3,0.8);   % Star-like object (help StarLikeSample to see the details of parameters)
                im2=StarLikeSample(2,Nxy,2,20,3,0.6);   % Star-like object (help StarLikeSample to see the details of parameters)

%                 im1=1-mat2gray(im1);
%                 im2=1-mat2gray(im2);

                offset=0.1*max(im1(:));
                im = (im1 + offset).*exp(1i*1*pi*im2);   % Complex object field
%                 im = exp(log(0.5)*(im1 + offset)).*exp(-1i*0.99*pi*(im2));   % Complex object field

            end
            H = LinOpWavePropKernel(lambda, Nxy, Nxy, pps, z);
            otf = H.mtf;

            y_noNoise=H*im;

        else
            % real data
            if is3D
                sz=[Nxy Nxy Nxy];
                im=StarLikeSample(3,Nxy,4,20,5,0.9);   % Star-like object (help StarLikeSample to see the details of parameters)
            else
                sz=[Nxy Nxy];                       % Image size
                im=StarLikeSample(2,Nxy,4,20,5,0.9);   % Star-like object (help StarLikeSample to see the details of parameters)
            end

            %% PSF
            lamb=561;                % Illumination wavelength
            res=30;                  % Resolution (nm)
            Na=1.4;                  % Objective numerical aperture
            fc=2*Na/lamb*res;        % cut-off frequency
            ll=linspace(-0.5,0,sz(1)/2+1);
            lr=linspace(0,0.5,sz(1)/2);
            [X,Y]=meshgrid([ll,lr(2:end)],[ll,lr(2:end)]);
            [th,rho]=cart2pol(X,Y);
            MTF=fftshift(1/pi*(2*acos(abs(rho)/fc)-sin(2*acos(abs(rho)/fc))).*(rho<fc));
            otf=real(ifft2(MTF));

            %% Data
            H=LinOpConv(MTF);
            y_noNoise=H*im;

        end

        if strcmp(noise,'Gaussian')
            if isCpx
                ynoise = 10^(-level/20).*random('Normal',zeros(size(y_noNoise)),ones(size(y_noNoise))) + ...
                    1i* 10^(-level/20).*random('Normal',zeros(size(y_noNoise)),ones(size(y_noNoise)))  ;
            else
                ynoise = max(y_noNoise(:)) .* 10^(-level/20).*random('Normal',zeros(size(y_noNoise)),ones(size(y_noNoise)));
            end
            y = y_noNoise + ynoise;
        elseif strcmp(noise,'Poisson')
            factor = level./mean(y_noNoise(:)) ;
            y_noNoise = y_noNoise.* factor;
            im = im.*factor;
            y = random('Poisson',y_noNoise);

            %          factor = level./mean((y_noNoise(:))) ;
            %     y_noNoise = (y_noNoise).* factor;
            %
            % %     factor = level./mean(real(y_noNoise(:))) ;
            % %     y_noNoise_re = real(y_noNoise).* factor;
            % %
            % %     factor = level./mean(imag(y_noNoise(:))) ;
            % %     y_noNoise_im = imag(y_noNoise).* factor;
            %
            %     im = im.*factor;
            %     y = random('Poisson', real(y_noNoise)) + 1i*random('Poisson', imag(y_noNoise));
            %
            % %     y = random('Poisson', y_noNoise_re) + 1i*imag(y_noNoise);
            %
            %     y = y./max(abs(y(:)));
        else
            error('Wrong type of noise');
        end
        disp(['SNR data : ',num2str(20*log10(norm(y_noNoise(:))/norm(y_noNoise(:)-y(:))))]);
    end
   
end
