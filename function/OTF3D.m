function [otf3d, psf3d, pupil3d] = OTF3D(Nx, Ny, Nz, lambda, dx, dy, dz, zd, sensor_size)

    % The object center is located at the origin of the coordinate system, z0 is the detect plane
    z = zd - ((1:Nz)- round(Nz/2)).*dz;  % Sampling of the z space
    
    otf3d = zeros(Ny, Nx, Nz);
    pupil3d = zeros(Ny, Nx, Nz);
    psf3d = zeros(Ny, Nx, Nz);  % Ni Chen

    rd = sensor_size/2;
    NA = rd./sqrt(z.^2 + rd.^2);
    for iz = 1:Nz
        [otf3d(:, :, iz), psf3d(:, :, iz), pupil3d(:, :, iz)] ...
         = RSOTF2D(Nx, Ny, z(iz), lambda, dx, dy, NA(iz));
    end
    
    psf3d = psf3d/max(max(max(abs(psf3d))));  % Ni Chen
end

% OTF and PSF of "rayleigh sommerfeld" diffraction (angular spectrum approach)
function [otf2d, psf2d, pupil2d] = RSOTF2D(Nx, Ny, z, lambda, dx, dy, NA)
    
    % Sampling at the spectrum plane
    KX = (ceil(-Nx/2):1:ceil(Nx/2-1))'.*(1/(Nx*dx));
    KY = (ceil(-Ny/2):1:ceil(Ny/2-1)).*(1/(Ny*dy));

    kx = repmat(KX, 1, Ny);
    ky = repmat(KY, Nx, 1);
    
    k = 1/lambda;
    term = k.^2 - (kx.^2 + ky.^2);
    term(term<0) = 0;

    otf2d = exp(1i*2*pi*z*sqrt(term));  % Transfer function of angular spectrum method
    pupil2d = ones(Nx, Ny);
    
    % ====================================== Modifed by Ni Chen ====================================
    psf2d = ifftshift(ifft2(fftshift(otf2d))); 
    
end

% % The Fresnel propagation
% function [otf2d, psf2d, pupil2d] = FROTF2D(Nx, Ny, z, lambda, dx, dy, NA)
% 
%     k = 1/lambda;
% 
%     KX = (ceil(-Nx/2):1:ceil(Nx/2-1))'.*(1/(Nx*dx));
%     KY = (ceil(-Ny/2):1:ceil(Ny/2-1)).*(1/(Ny*dy));
% 
%     kx = repmat(KX, 1, Ny);
%     ky = repmat(KY, Nx, 1);
%         
%     term = k.^2 - (kx.^2 + ky.^2);
%     term(term<0) = 0;
% 
%     otf2d = exp(1i*2*pi*z*sqrt(term));  % Transfer function of angular spectrum method
%     pupil2d = ones(Nx, Ny);
%     
%     % ====================================== Modifed by Ni Chen ====================================
%     psf2d = ifftshift(ifft2(fftshift(otf2d))); 
%     
% end