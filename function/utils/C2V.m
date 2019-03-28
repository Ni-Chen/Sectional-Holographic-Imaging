% Complex to vector
% x_vec = @(x_complex) C2V([real(x_complex); imag(x_complex)]) 
function x_vec = C2V(x_complex)

	x_vec = [real(x_complex); imag(x_complex)];

end