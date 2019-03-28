% Vector to complex
function x_complex = V2C(x_vec)
    [Ny, Nx] = size(x_vec);
    x_complex = zeros(Ny*Nx/2, 1);
    x_complex = x_vec(1:Ny*Nx/2) + 1i*x_vec(Ny*Nx/2+1:Ny*Nx);
end