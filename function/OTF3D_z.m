function [otf3d, psf3d, pupil3d] = OTF3D_z(Ny, Nx, lambda, ppxy, z)

    Nz = length(z);
    otf3d = zeros(Ny, Nx, Nz);
    pupil3d = zeros(Ny, Nx, Nz);
    psf3d = zeros(Ny, Nx, Nz);  

    Lx = Nx*ppxy;
    rd = Lx/2;
    NA = rd./sqrt(z.^2 + rd.^2);
    for iz = 1:Nz
        [otf3d(:,:,iz), psf3d(:,:,iz), pupil3d(:,:,iz)] = RSOTF2D(Ny, Nx, z(iz), lambda, ppxy);
    end
    
    psf3d = psf3d/max(max(max(abs(psf3d))));
end

% OTF and PSF of "rayleigh sommerfeld" diffraction (angular spectrum approach)
function [otf2d, psf2d, pupil2d] = RSOTF2D(Ny, Nx, z, lambda, ppxy)
    iFT2 = @(x) ifftshift(ifft2(fftshift(x)));
    % Sampling at the spectrum plane
%     KX = (ceil(-Nx/2):1:ceil(Nx/2-1))'.*(1/(Nx*ppxy));
%     KY = (ceil(-Ny/2):1:ceil(Ny/2-1)).*(1/(Ny*ppxy));
%     kx = repmat(KX, 1, Ny);
%     ky = repmat(KY, Nx, 1);
        
%     KY = (ceil(-Nx/2):1:ceil(Nx/2-1))'.*(1/(Nx*ppxy));
%     KX = (ceil(-Ny/2):1:ceil(Ny/2-1)).*(1/(Ny*ppxy));
%     kx = repmat(KX, Ny, 1);
%     ky = repmat(KY, 1, Nx);

    Y = ((1:Ny)-round(Ny/2))*ppxy;
    X = ((1:Nx)-round(Nx/2))*ppxy;      
    [x, y] = meshgrid(X, Y);
          
    KY = ((1:Ny)-round(Ny/2))*(1/(Ny*ppxy));
    KX = ((1:Nx)-round(Nx/2))*(1/(Nx*ppxy));      
    [kx, ky] = meshgrid(KX, KY);
    
    % ========================================== OTF =============================================== 
    k = 1/lambda;
    term = k.^2 - (kx.^2 + ky.^2);
    term(term<0) = 0;

    otf2d = exp(1i*2*pi*z*sqrt(term));  % Transfer function of angular spectrum method
      
    % ========================================== PSF =============================================== 
%     psf2d = iFT2(otf2d);    % analytical will be better (todo)
    r = sqrt(x.^2 + y.^2 + z^2);
    psf2d = exp(1i*2*pi/lambda*r)./r;
    
    % ========================================= Circular Pupil ============================================== 
    rs = Nx*ppxy/2;  % radius of aperture in spatial domain
    pupil2d = double(sqrt(x.^2 + y.^2)<=rs);
    
    rf  = rs/lambda/z;  % radius of aperture in Frequency domain
    pupil2d_fft = double(sqrt(kx.^2 + ky.^2)<=rf).*exp(1i*pi*lambda*z*(kx.^2 + ky.^2));  % Gaussian plane wave
    
    % ========================================= Rectangle Pupil ============================================== 
    rs = Nx*ppxy/2;  % radius of aperture in spatial domain
    pupil2d = double(abs(x)<=rs & abs(y)<=rs);
    
    rf  = rs/lambda/z;  % radius of aperture in Frequency domain
    pupil2d_fft = double(abs(kx)<=rf & abs(ky)<=rf).*exp(1i*pi*lambda*z*(kx.^2 + ky.^2)); % Gaussian plane wave
    pupil2d_fft = double(abs(kx)<=rf & abs(ky)<=rf); % rectangle aperture
    pupil2d = pupil2d_fft;
    
%     pupil2d = ones(Ny, Nx);    % No aperture case, Enough for pure simulation cases
end

% The Fresnel scaled "rayleigh sommerfeld" diffraction (angular spectrum approach)
function [otf2d, psf2d, pupil2d] = FSRSOTF2D(Nx, Ny, z, lambda, ppxy, NA)
      % Sampling at the spectrum plane
    KX = (ceil(-Nx/2):1:ceil(Nx/2-1))'.*(1/(Nx*ppxy));
    KY = (ceil(-Ny/2):1:ceil(Ny/2-1)).*(1/(Ny*ppxy));

    kx = repmat(KX, 1, Ny);
    ky = repmat(KY, Nx, 1);
    
    k = 1/lambda;
    term = k.^2 - (kx.^2 + ky.^2);
    term(term<0) = 0;

    otf2d = exp(1i*2*pi*z*sqrt(term));  % Transfer function of angular spectrum method
    pupil2d = ones(Ny, Nx);
    
    % ====================================== Modifed by Ni Chen ====================================
%     psf2d = ifftshift(ifft2(fftshift(otf2d)));    
    psf2d = (ifft2(fftshift(otf2d)));    
    
end

% Fresnel diffraction (angular spectrum approach)
function [otf2d, psf2d, pupil2d] = FROTF2D(Ny, Nx, z, lambda, ppxy, NA)
    
    % Sampling at the spectrum plane
    KX = (ceil(-Nx/2):1:ceil(Nx/2-1))'.*(1/(Nx*ppxy));
    KY = (ceil(-Ny/2):1:ceil(Ny/2-1)).*(1/(Ny*ppxy));

    kx = repmat(KX, 1, Ny);
    ky = repmat(KY, Nx, 1);
    
    k = 1/lambda;
    term = k.^2 - (kx.^2 + ky.^2);
    term(term<0) = 0;

    otf2d = exp(1i*2*pi*z*sqrt(term));  % Transfer function of angular spectrum method
    pupil2d = ones(Ny, Nx);
    
    % ====================================== Modifed by Ni Chen ====================================
    psf2d = ifftshift(ifft2(fftshift(otf2d)));    
end
