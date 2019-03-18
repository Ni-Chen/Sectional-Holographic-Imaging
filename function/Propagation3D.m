%{
    Propagate 3D volume field to a single plane
%}
function plane_field = Propagation3D(vol_field, otf3d, pupil3d)
    
    [Ny, Nx, ~] = size(otf3d);

    % Field of plane wave, corresponds to Fourier transform of point light source
    pinhole = ifftshift(fft2(fftshift(ones(Ny, Nx))));  % delta = FT(1)
    point_field = pupil3d.*conj(otf3d).*pinhole;
    plane_wave = ifftshift(ifft2(fftshift(point_field)));  
    
    vol_field_ft = ifftshift(fft2(fftshift(vol_field.*plane_wave)));  % Illuminate object with plane wave
  
    plane_field_ft = sum(vol_field_ft.*otf3d.*pupil3d, 3);  % intergration along z axis
    plane_field = ifftshift(ifft2(fftshift(plane_field_ft)));
end