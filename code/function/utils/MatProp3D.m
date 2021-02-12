%{
    Propagate a 3D volume field to a single plane
%}
function field2d = MatProp3D(field3d, otf3d, pupil3d)
    FT2 = @(x) ifftshift(fft2(fftshift(x)));
    iFT2 = @(x) ifftshift(ifft2(fftshift(x)));

%     [Ny, Nx, ~] = size(otf3d);
% 
%     % Field of plane wave, corresponds to Fourier transform of point light source
%     pinhole = FT2(ones(Ny, Nx));  % delta = FT(1)
%     point_field = pupil3d.*conj(otf3d).*pinhole;
%     plane_wave = iFT2(point_field);  
%     
%     vol_field_ft = FT2(field3d.*plane_wave);  % Illuminate object with plane wave
%   
%     plane_field_ft = sum(vol_field_ft.*otf3d.*pupil3d, 3);  % intergration along z axis
%     field2d = iFT2(plane_field_ft);
    
    field2d = iFT2(FT2(field3d).*otf3d.*pupil3d);
    field2d = sum(field2d, 3);
%     field2d = 2*real(field2d);
%     u = iFT2(FT2(eta.*E).*phase);
end