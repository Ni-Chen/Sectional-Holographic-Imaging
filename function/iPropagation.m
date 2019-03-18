function vol_field = iPropagation(plane_field, otf3d, pupil)
    [Ny, Nx, ~] = size(otf3d);
    
    plane_field = reshape(V2C(plane_field), Nx, Ny);

    vol_field = iPropagation3D(plane_field, otf3d, pupil);

    vol_field = C2V(vol_field(:));
end

% function vol_field = iPropagation3D(plane_field, otf3d, pupil)
%     [Ny, Nx, ~] = size(otf3d);
%    
%     cEsp = conj(ifftshift(ifft2(fftshift(conj(plane_field)))));  % Modified by Ni Chen, delete real
%     cEs = conj(otf3d).*conj(pupil).*cEsp;
% 
%     E0 = ones(Nx, Ny);  % Plane wave
%     pinhole = ifftshift(fft2(fftshift(E0)));  % point object
%     
%     cE = pupil.*conj(otf3d).*pinhole;  % cE = pinhole.*cOTF3D(:, :, i).*pupil(:, :, i);
%      
%     E = ifftshift(ifft2(fftshift(cE)));
%     vol_field = ifftshift(fft2(ifftshift(conj(cEs))));
% 
%     vol_field = conj(E).*conj(vol_field);
% end