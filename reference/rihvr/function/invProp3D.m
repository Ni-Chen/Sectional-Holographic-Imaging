% Reconstruction of hologram to form a 3D optical field
%
% Inputs:
%   holo    : ny-by-nx 2D real data (often the residual hologram)
%   params  : Input parameters structure. Must contain the following
%       nx, ny, nz  : Size of the volume (in voxels) in x, y, and z
%       z_list      : List of reconstruction planes (um)
%       resolution  : Pixel size (um) of the image
%       wavelength  : Illumination wavelength (um)
% 
% Outputs:
%   volume  : 3D estimated optical field, complex valued

function volume = invProp3D(holo, params)
    %holo = real(holo);

    % Abreviations
    nx = params.nx;
    ny = params.ny;
    nz = params.nz;
    reso = params.resolution;
    lambda = params.wavelength;
    z_list = params.z_list;


    % Constant frequencies
    % write it this way to avoid fftshifts
    x = 0:(nx-1);
    y = 0:(ny-1);
    [X,Y] = meshgrid(x,y);
    fx = (mod(X + nx/2, nx) - floor(nx/2)) / nx;
    fy = (mod(Y + ny/2, ny) - floor(ny/2)) / ny;
    f2 = fx.*fx + fy.*fy;

    sqrt_input = 1 - f2*(lambda/reso)^2;
    sqrt_input(sqrt_input < 0) = 0;

    H = -2*pi*1i*sqrt(sqrt_input)/lambda;

    Fholo = fft2(holo);
    volume = zeros(ny, nx, nz);

    for zid = 1:length(z_list)
        z = z_list(zid);
        phase = exp(2*pi*1i*z/lambda);

        Hz = exp(z*H);
        Fplane = Fholo .* Hz;
        plane = ifft2(Fplane) .* phase;
        volume(:,:,zid) = (plane);
    end
end

