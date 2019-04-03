function [otf3d, psf3d, pupil3d] = OTF3D_z(Ny, Nx, lambda, ppxy, z)

    % The object center is located at the origin of the coordinate system, z0 is the detect plane
%     z = zd - ((1:Nz)- round(Nz/2)).*dz;  % Sampling of the z space, z = zd + (0:Nz-1).*dz;  
%     z = zd + (0:Nz-1).*dz;   % inline twin image eliminating
    
    Nz = length(z);
    otf3d = zeros(Ny, Nx, Nz);
    pupil3d = zeros(Ny, Nx, Nz);
    psf3d = zeros(Ny, Nx, Nz);  % Ni Chen

    Lx = Nx*ppxy;
    rd = Lx/2;
    NA = rd./sqrt(z.^2 + rd.^2);
    for iz = 1:Nz
        [otf3d(:, :, iz), psf3d(:, :, iz), pupil3d(:, :, iz)] ...
         = RSOTF2D(Ny, Nx, z(iz), lambda, ppxy, NA(iz));
    end
    
    psf3d = psf3d/max(max(max(abs(psf3d))));  % Ni Chen
end

% OTF and PSF of "rayleigh sommerfeld" diffraction (angular spectrum approach)
function [otf2d, psf2d, pupil2d] = RSOTF2D(Ny, Nx, z, lambda, ppxy, NA)
    
    % Sampling at the spectrum plane
%     KX = (ceil(-Nx/2):1:ceil(Nx/2-1))'.*(1/(Nx*ppxy));
%     KY = (ceil(-Ny/2):1:ceil(Ny/2-1)).*(1/(Ny*ppxy));
%     kx = repmat(KX, 1, Ny);
%     ky = repmat(KY, Nx, 1);
        
%     KY = (ceil(-Nx/2):1:ceil(Nx/2-1))'.*(1/(Nx*ppxy));
%     KX = (ceil(-Ny/2):1:ceil(Ny/2-1)).*(1/(Ny*ppxy));
%     kx = repmat(KX, Ny, 1);
%     ky = repmat(KY, 1, Nx);
        
    KY = ((1:Ny)-round(Ny/2))*(1/(Ny*ppxy));
    KX = ((1:Nx)-round(Nx/2))*(1/(Nx*ppxy));   
   
    [kx, ky] = meshgrid(KX, KY);
    
    k = 1/lambda;
    term = k.^2 - (kx.^2 + ky.^2);
    term(term<0) = 0;

    otf2d = exp(1i*2*pi*z*sqrt(term));  % Transfer function of angular spectrum method
    pupil2d = ones(Ny, Nx);
    
    % ====================================== Modifed by Ni Chen ====================================
    psf2d = ifftshift(ifft2(fftshift(otf2d)));    
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
    psf2d = ifftshift(ifft2(fftshift(otf2d)));    
    
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
