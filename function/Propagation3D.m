%{
    Propagate a 3D volume field to a single plane
%}
function plane_field = Propagation3D(vol_field, otf3d, pupil3d, holo_type)
    
    [Ny, Nx, ~] = size(otf3d);

    % Field of plane wave, corresponds to Fourier transform of point light source
    pinhole = ifftshift(fft2(fftshift(ones(Ny, Nx))));  % delta = FT(1)
    point_field = pupil3d.*conj(otf3d).*pinhole;
    plane_wave = ifftshift(ifft2(fftshift(point_field)));  
    
    vol_field_ft = ifftshift(fft2(fftshift(vol_field.*plane_wave)));  % Illuminate object with plane wave
  
    plane_field_ft = sum(vol_field_ft.*otf3d.*pupil3d, 3);  % intergration along z axis
    plane_field = ifftshift(ifft2(fftshift(plane_field_ft)));

    if nargin>3 && strcmp(holo_type, 'inline')
        %{
          Gabor hologram: I = |R|^2 + |O|^2 + R*O + RO*
          Substract background: I_prime = |O|^2 + R*O + RO* = 2 real(O) + |O|^2 = 2 real(O) + err
        %}
        plane_field = 2*real(plane_field);
    end

%     % Filter
%     sigma = 1;
%     gausFilter = fspecial('gaussian', [3,3], sigma);
%     plane_field = imfilter(plane_field, gausFilter, 'replicate');
end