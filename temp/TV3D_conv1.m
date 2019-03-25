function TV = TV3D_conv(x)
%{
    V(x) = sum{|x_{i+1}?x_i|}
    V(x) = sum{|x_{i+1,j}?x_{i,j}| + |x_{i,j+1}?x_{i,j}|}
%}
    [nx,  ny,  nz] = size(x);
    TV = zeros(nx,  ny,  nz,  3);
    
    TV(:, :, :, 1) = circshift(x, [-1 0 0]) - x;
    TV(nx, :, :, 1) = 0.0;

    TV(:, :, :, 2) = circshift(x, [0 -1 0]) - x;
    TV(:, ny, :, 2) = 0.0;

    TV(:, :, :, 3) = circshift(x, [0 0 -1]) - x;
    TV(:, :, nz,  3) = 0.0;
    
    TV(:, :, :, 3) = TV(:, :, :, 3).*(1.0);
    
end