% Free-space wave propagation
function plane_field = Propagation(obj3d, otf3d, pupil, holo_type)
    [Ny, Nx, Nz] = size(otf3d);

    obj3d = reshape(V2C(obj3d), Ny, Nx, Nz);

    plane_field = Propagation3D(obj3d, otf3d, pupil, holo_type);

    plane_field = C2V(plane_field(:));
end

