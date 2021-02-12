function y = L1phi(x)

    x_complex = V2C(x(:));
    y = sum(sqrt(real(x_complex).^2 + imag(x_complex).^2));
    
end
