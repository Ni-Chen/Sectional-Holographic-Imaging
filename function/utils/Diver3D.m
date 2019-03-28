function divp = Diver3D(TV)
% 3D divergence
	n = size(TV);

	y_shift  =  circshift(TV(:, :, :, 1), [1 0 0 0]);
	yy = TV(:, :, :, 1) - y_shift;
	yy(1, :, :) = TV(1, :, :, 1);
	yy(n(1), :, :) = -y_shift(n(1), :, :);

	x_shift = circshift(TV(:, :, :, 2), [0 1 0 0]);
	yx = TV(:, :, :, 2) - x_shift;
	yx(:, 1, :) = TV(:, 1, :, 2);
	yx(:, n(2), :) = -x_shift(:, n(2), :);

	z_shift = circshift(TV(:, :, :, 3), [0 0 1 0]);
	yz = TV(:, :, :, 3) - z_shift;
	yz(:, :, 1) = TV(:, :, 1, 3);
	yz(:, :, n(3)) = -z_shift(:, :, n(3));

	divp = yx + yy + yz;
end


