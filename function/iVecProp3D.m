function field3d_vec = iVecProp3D(field2d_vec, otf3d, pupil, holo_type)
    [Ny, Nx, ~] = size(otf3d);
    
    field2d_vec_complex = V2C(field2d_vec);
    field2d_mat = reshape(field2d_vec_complex, Ny, Nx);

    field3d_mat = iMatProp3D(field2d_mat, otf3d, pupil, holo_type);

    field3d_vec = C2V(field3d_mat(:));
end
