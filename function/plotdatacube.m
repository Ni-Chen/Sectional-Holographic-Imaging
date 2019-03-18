function spatialcube = plotdatacube(data)
    data = padarray(data, [3 3], max(data(:)), 'both');
%     data = padarray(data, [3 3], min(data(:)), 'both');

    [Ny, Nx, Nz] = size(data);
    
    if Nz>=20
        cols = 10;
    else
        cols = 5;
    end
    
    rows = ceil(Nz/cols);

    spatialcube = zeros(rows*Ny, cols*Nx);

    figscount = 1;
    for r = 1:rows
        for c = 1:cols
            if figscount <= Nz
                spatialcube((r-1)*Ny+1:(r-1)*Ny+Ny, (c-1)*Nx+1:(c-1)*Nx+Nx) = squeeze(data(:, :, figscount)); %Nz-figscount+1));
            else
                spatialcube((r-1)*Ny+1:(r-1)*Ny+Ny, (c-1)*Nx+1:(c-1)*Nx+Nx) = zeros(Ny, Nx);
            end
            figscount = figscount+1;
        end
    end

end