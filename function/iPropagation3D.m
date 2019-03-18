%{
    Propagate a single plane field to the 3D object space
%}
function vol_field = iPropagation3D(plane_field, otf3d, pupil3d)
    [Ny, Nx, ~] = size(otf3d);
   
    plane_field_ft = conj(ifftshift(ifft2(fftshift(conj(plane_field)))));  
%     plane_field_ft = ifftshift(ifft2(fftshift(plane_field)));  
    vol_field_ft = conj(otf3d).*conj(pupil3d).*plane_field_ft;  % Illuminate with plane wave, and back propagation

    pinhole = ifftshift(fft2(fftshift(ones(Ny, Nx))));   % delta = FT(1)    
    point_field = pupil3d.*conj(otf3d).*pinhole;       
    point_field_ft = ifftshift(ifft2(fftshift(point_field)));
    
    volume_field = ifftshift(fft2(ifftshift(conj(vol_field_ft))));

    vol_field = conj(point_field_ft).*conj(volume_field);
end