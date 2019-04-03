%{
    Propagate a single complex field to the 3D object space
%}
function field3d = iMatProp3D(field2d, otf3d, pupil3d, holo_type)
    [Ny, Nx, ~] = size(otf3d);

    if  nargin>3 && strcmp(holo_type, 'inline')
        %{
          Gabor hologram: I = |R|^2 + |O|^2 + R*O + RO*
          Substract background: I_prime = |O|^2 + R*O + RO* = 2 real(O) + |O|^2 = 2 real(O) + err
        %}
        field2d = 2*real(field2d);
    end

    plane_field_ft = conj(ifftshift(ifft2(fftshift(conj(field2d)))));  
    vol_field_ft = conj(otf3d).*conj(pupil3d).*plane_field_ft;  % Illuminate with plane wave, and back propagation

    pinhole = ifftshift(fft2(fftshift(ones(Ny, Nx))));  % delta = FT(1)    
    point_field = pupil3d.*conj(otf3d).*pinhole;
    point_field_ft = ifftshift(ifft2(fftshift(point_field)));
    
    volume_field = ifftshift(fft2(ifftshift(conj(vol_field_ft))));

    field3d = conj(point_field_ft).*conj(volume_field);
    
end