% Vector to complex
function x_complex = V2C(x_vec)
[Nyy, Nxx] = size(x_vec);
x_complex = zeros(Nyy/2, Nxx);
x_complex = x_vec(1:Nyy/2, :) + 1i*x_vec(Nyy/2+1:Nyy, :);
end