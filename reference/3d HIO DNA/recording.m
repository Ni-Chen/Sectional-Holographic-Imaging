%{
Simulates 3d angular-averaged diffraction pattern of discontinuous double helix. 

Citation for this code/algorithm or any of its parts:
T. Latychevskaia and H.-W. Fink, "Three-dimensional double helical DNA structure directly revealed
from its X-ray fiber diffraction pattern by iterative phase retrieval", OE 26(23):30991-31017(2018)

To reconstruct the 3D diffraction patter, start "b_HIO_3d_sequence.m" and
select the 2D slices simulated by "a_simulate_3d_dp_no_phase.m".
Obtain several reconstructions, 300 iterations for each reconstruction is sufficient.
%}

close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 128;                    % number of pixels
r = 10^(-9);                % DNA radius in meter
H = 3.4*10^(-9);            % helix large period in meter
h = 0.34*10^(-9);           % base-to-base period in meter
delta0 = 0.17*10^(-9);      % pixel size in meter
Dz = 1.275*10^(-9);         % z-shift between the helices  

deltak = 2*pi/(100*delta0); % pixels size in K-domain
kz_h = 2*pi/h;              % k_z position corresponding to h
kz_h_px = kz_h/deltak;      % k_z in pixels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

[kx, ky] = meshgrid(((1:N) - N/2 - 1)*deltak, ((1:N) - N/2 - 1)*deltak);
R = sqrt(kx.^2 + ky.^2);
a = R*r;

U1 = zeros(N, N, N); % term1
kz1 = ((1:N) - N/2 - 1)*deltak*H/(2*pi);

U2 = zeros(N, N, N); % term2, the first diffraction order because of atomic lattice
kz2 = ((1:N) - N/2 - 1 - kz_h_px)*deltak*H/(2*pi);

U3 = zeros(N, N, N); % term3, the minus first diffraction order because of atomic lattice
kz3 = ((1:N) - N/2 - 1 - kz_h_px)*deltak*H/(2*pi);

for kk = 1:N
    n1 = kz1(kk);
    nn1 = round(n1);
    n_nn1 = abs(n1 - round(nn1));
    if ((n1 == 0) || (n_nn1 < 10^(-8)))
        U1(:, :, kk) = besselj(nn1,a);
    end
    
    n2 = kz2(kk);
    nn2 = round(n2);
    n_nn2 = abs(n2 - round(nn2));
    if ((n2 == 0) || (n_nn2 < 10^(-8)))
        U2(:, :, kk) = besselj(nn2,a);
    end
    
    n3 = kz3(kk);
    nn3 = round(n3);
    n_nn3 = abs(n3 - round(nn3));
    if ((n3 == 0) || (n_nn3 < 10^(-8)))
        U3(:, :, kk) = besselj(nn3,a);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modulating function because of the second helix
M = zeros(N, N, N);

% x = N/2 - (1:N);
% y = N/2 - (1:N);
kz = ((1:N) - N/2 -1)*deltak;
temp = 2*(1 + cos(kz*Dz));

for kk = 1:N
      M(:,:,kk) = temp(kk);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I1 = abs(U1).^2 + abs(U2).^2 + abs(U3).^2;
dp = I1.*M;
Ishow(:,:)= dp(:,:,N/2+1);
imshow(rot90(Ishow), []);

toc

figure; show3d(abs(U1), 0.05); axis normal;
figure; show3d(abs(U2), 0.05); axis normal;
figure; show3d(abs(U3), 0.05); axis normal;

save('3d_dp.mat', 'dp');


