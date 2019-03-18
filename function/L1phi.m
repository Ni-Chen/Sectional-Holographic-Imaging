function y = L1phi(x)

    X= V2C(x(:));
    y=sum(sqrt(real(X).^2+imag(X).^2));
    
end
