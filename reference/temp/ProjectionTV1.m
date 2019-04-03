function p = ProjectionTV(g, tau, lam, iter)

    [nx, ny, nz] = size(g);
    pn = zeros(nx, ny, nz, 3);
    div_pn = zeros(nx, ny, nz);
    b = pn;

    for i = 1:iter

        a = TV3D_conv(div_pn-g./lam);

        b(:, :, :, 1) = sqrt(a(:, :, :, 1).^2 + a(:, :, :, 2).^2 + a(:, :, :, 3).^2);
        b(:, :, :, 2) = b(:, :, :, 1);
        b(:, :, :, 3) = b(:, :, :, 1);
        
        pn = (pn+tau.*a)./(1.0+tau.*b);

        div_pn = Div3D(pn);
    end

    p = lam.*Div3D(pn);
end

function y = Div3D(TV)

	n = size(TV);

	x_shift  =  circshift(TV(:, :, :, 1), [1 0 0 0]);
	yx = TV(:, :, :, 1) - x_shift;
	yx(1, :, :) = TV(1, :, :, 1);
	yx(n(1), :, :) = -x_shift(n(1), :, :);

	y_shift = circshift(TV(:, :, :, 2), [0 1 0 0]);
	yy = TV(:, :, :, 2) - y_shift;
	yy(:, 1, :) = TV(:, 1, :, 2);
	yy(:, n(2), :) = -y_shift(:, n(2), :);

	z_shift = circshift(TV(:, :, :, 3), [0 0 1 0]);
	yz = TV(:, :, :, 3) - z_shift;
	yz(:, :, 1) = TV(:, :, 1, 3);
	yz(:, :, n(3)) = -z_shift(:, :, n(3));

	y = yx + yy + yz;
end

