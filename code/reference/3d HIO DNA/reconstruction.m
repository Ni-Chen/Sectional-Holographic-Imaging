%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
3D ITERATIVE PHASE RETRIEVAL BY HYBRID INPUT OUTPUT ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Citation for this code/algorithm or any of its parts:
T. Latychevskaia and H.-W. Fink, "Three-dimensional double helical DNA structure directly revealed
from its X-ray fiber diffraction pattern by iterative phase retrieval", OE 26(23):30991-31017(2018)

To reconstruct the 3D diffraction patter, start "b_HIO_3d_sequence.m".
Obtain several reconstructions, 300 iterations for each reconstruction is sufficient.
%}

close all; clear; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.9;             % parameter in HIO algorithm
Iterations = 100;       % number of iterative loops, typically 200 iterations are enough
Reconstructions = 5;    % reconstuctions
N = 128;                % number of pixels
support_r = 9;          % radius of the support mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reading z-slices of 3d diffraction pattern

tic

load('3d_dp.mat');
figure; show3d(abs(dp), 0.05); axis normal;

Ishow(:,:)= dp(N/2,:,:);
imshow(rot90(Ishow), []);

amplitude_experimental = sqrt(dp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% creating support in object domain
[x, y] = meshgrid(N/2 - (1:N));
temp = (sqrt(x.^2 + y.^2) < support_r);
support = repmat(temp, [1 1 N]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reconstructions loop

for iRec = 1:Reconstructions    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % creating initial complex-valued field distribution at detector plane
    phase = (2*rand(N,N,N)-1)*pi;
    
    field_detector_0 = amplitude_experimental.*exp(1i*phase);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % getting initial object distribution
    
    object_0 = IFT3Dc(field_detector_0);
    gk = real(object_0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % iterative loop
    qq = 0;    
    for iLoop = 1:Iterations
        
        disp(['Reconstruction:', num2str(iRec), ', Iteration:', num2str(iLoop)]);
        
        field_detector_updated = FT3Dc(gk);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calculating error
        % amplitude_updated = abs(field_detector_updated);
        
        % replacing updated amplitude for measured amplitude
        field_detector_updated = amplitude_experimental.*exp(1i*angle(field_detector_updated));
        
        % getting updated object distribution
        gk_prime = real(IFT3Dc(field_detector_updated));
        Er(iLoop) = errEst3D(gk_prime, support);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gk_prime_show(:,:) = gk_prime(:, N/2, :);
        imshow(rot90(gk_prime_show), []);

        % HIO filter        
        temp = ((gk_prime > 0) & (support > 0));
        gk1= temp.*gk_prime + (1-temp).*(gk - beta.*gk_prime);
        
        gk = gk1;
        
    end % Iterations
    
    field_detector_updated = FT3Dc(gk);
    % amplitude_updated = abs(field_detector_updated);
    % Er = errEst3D(amplitude_experimental, amplitude_updated);
    
    save('rec_obj.mat', 'gk');
end 
toc

figure; plot(Er);
reconstruction_show(:,:)= gk(:,:,N/2);
figure;imshow(rot90(reconstruction_show), []);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D PLOT OF RECONSTRUCTED OBJECT
figure; show3d(abs(gk), 0.05); axis normal;
