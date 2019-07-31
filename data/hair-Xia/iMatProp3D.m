%{
    Propagate a single complex field to the 3D object space
%}
function field3d = iMatProp3D(field2d, otf3d, pupil3d)
    FT2 = @(x) ifftshift(fft2(fftshift(x)));
    iFT2 = @(x) ifftshift(ifft2(fftshift(x)));
    
%     [Ny, Nx, ~] = size(otf3d);
% 
%     plane_field_ft = conj(iFT2(conj(field2d)));
%     vol_field_ft = conj(otf3d).*conj(pupil3d).*plane_field_ft;  % Illuminate with plane wave, and back propagation
% 
%     pinhole = FT2(ones(Ny, Nx));  % delta = FT(1)   
%     point_field = pupil3d.*conj(otf3d).*pinhole;
%     point_field_ft = iFT2(point_field);
%     
%     volume_field = FT2(conj(vol_field_ft));
% 
%     field3d = conj(point_field_ft).*conj(volume_field);  
%     
    field3d = iFT2(FT2(field2d).*conj(otf3d).*pupil3d);
%     field3d = FT2(iFT2(field2d).*conj(otf3d).*pupil3d);
%     field3d = FT2(conj(conj(otf3d).*conj(iFT2(field2d))));
    
end