function TV = Grad3D(x)
%{
    Discrete Total Variation is:
    1D: V(x) = sum{|x_{i+1}-x_i|}
    2D: V(x) = sum{|x_{i+1,j}-x_{i,j}| + |x_{i,j+1}-x_{i,j}|}
    3D: V(x) = 
%}
    [ny,  nx,  nz] = size(x);
    TV = zeros(ny,  nx,  nz,  3);
    
    TV(:, :, :, 1) = circshift(x, [-1 0 0]) - x;
    TV(ny, :, :, 1) = 0.0;

    TV(:, :, :, 2) = circshift(x, [0 -1 0]) - x;
    TV(:, nx, :, 2) = 0.0;

    TV(:, :, :, 3) = circshift(x, [0 0 -1]) - x;
    TV(:, :, nz, 3) = 0.0;
    
%     TV(:, :, :, 3) = TV(:, :, :, 3).*(1.0);    
end