% Free-space wave propagation
function plane_field = Propagation(obj3d, otf3d, pupil)
    [Ny, Nx, Nz] = size(otf3d);

    obj3d = reshape(V2C(obj3d), Nx, Ny, Nz);

    plane_field = Propagation3D(obj3d, otf3d, pupil);

    plane_field = C2V(plane_field(:));
end

% %{
%     Propagate 3D to one single plane gabor hologram	
% %}
% function plane_field = Propagation3D(obj3d, otf3d, pupil)
%     
%     [Ny, Nx, ~] = size(otf3d);
% 
%     E0 = ones(Nx, Ny); 
%     pinhole = ifftshift(fft2(fftshift(E0))); % point object
% 
%     % cE = pinhole.*cOTF3D(:, :, i).*pupil(:, :, i);
%     cE = pupil.*conj(otf3d).*pinhole;
%     
%     temp = ifftshift(ifft2(fftshift(cE)));
%     E_vol = ifftshift(fft2(fftshift(obj3d.*temp)));  % g = iFT(O*OTF)
%   
%     plane_field = sum(E_vol.*otf3d.*pupil, 3);  % intergration along z axis
%     plane_field = ifftshift(ifft2(fftshift(plane_field)));
%     
% end