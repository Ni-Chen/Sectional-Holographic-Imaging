function [holo_shrink, pps_new] = holoResize(holo, pps, Ny_prime, pad_size)
    [Ny, Nx] = size(holo);

    shrink_ratio = Ny_prime/Ny;
    holo_shrink = imresize(holo, shrink_ratio);
    holo_shrink(1:pad_size, :) = 0;
    holo_shrink(Ny_prime-pad_size:Ny_prime, :) = 0;
    holo_shrink(:, 1:pad_size) = 0;
    holo_shrink(:, Ny_prime-pad_size:Ny_prime) = 0;

    pps_new = pps/shrink_ratio;
end
