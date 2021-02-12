% Free-space wave propagation, 3D to 2D
function field2d_vec = VecProp3D(field3d_vec, otf3d, pupil, holo_type)
    [Ny, Nx, Nz] = size(otf3d);

    field3d_vec_complex = V2C(field3d_vec);
    field3d_mat = reshape(field3d_vec_complex, Ny, Nx, Nz);

    field2d_mat = MatProp3D(field3d_mat, otf3d, pupil);  % 3D volume propagation

%     field2d_vec = C2V(field2d_mat);
    field2d_vec = C2V(field2d_mat(:));
end

