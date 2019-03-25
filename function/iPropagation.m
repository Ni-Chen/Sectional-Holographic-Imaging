function vol_field = iPropagation(plane_field, otf3d, pupil, holo_type)
    [Ny, Nx, ~] = size(otf3d);
    
    plane_field = reshape(V2C(plane_field), Ny, Nx);

    vol_field = iPropagation3D(plane_field, otf3d, pupil, holo_type);

    vol_field = C2V(vol_field(:));
end
