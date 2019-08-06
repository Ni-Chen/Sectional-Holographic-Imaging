% Simulate forward formation of a hologram from a 3D optical field
%
% Inputs:
%   volume  : 3D estimated optical field, generally complex valued
%   params  : Input parameters structure. Must contain the following
%       nx, ny, nz  : Size of the volume (in voxels) in x, y, and z
%       z_list      : List of reconstruction planes (um)
%       resolution  : Pixel size (um) of the image
%       wavelength  : Illumination wavelength (um)
% 
% Outputs:
%   holo    : ny-by-nx 2D real hologram (estimated image)

function holo = Prop3D(volume, params)

    % Abreviations
    nx = params.nx;
    ny = params.ny;
    reso = params.resolution;
    lambda = params.wavelength;
    z_list = params.z_list;

    Fholo = zeros(ny, nx);

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

    for zid = 1:length(z_list)
        z = -z_list(zid);
        phase = exp(1i*2*pi*z/lambda);

        Fplane = fft2(volume(:,:,zid)*phase);
        Hz = exp(z*H);
        Fholo = Fholo + (Fplane .* Hz);
    end

    % Can sum in Fourier domain and only inverse after
    holo = ifft2(Fholo);
    holo = 1*real(holo);
    % holo = (holo).^2;
end