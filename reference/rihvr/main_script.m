%   FASTA Copyright: Tom Goldstein, 2014.
%   For more details, see the paper
%   T. Goldstein, C. Studer, and R. Baraniuk, “A Field Guide to Forward-
%   Backward Splitting with a FASTA Implementation,” arXiv:1411.3406, 2014.
%

close all;
clearvars();

addpath(genpath('.\function\'));

rng(0);

if ~isdir('outputs')
    mkdir('outputs');
end

%% Define Parameters
% Optional parameters, set to empty because not included in all mat files
display_planes = [];    % List of plane indices to display at end
z_list = [];            % List of plane locations to reconstruct
opts = [];

% Choose regularization method
% method = 'L1';          % L1-norm enforces sparsity
% method = 'TV';          % Total Variation enforces smoothness (3D)
% method = 'TV2D';        % 2D variant of TV (XY plane only)
% method = 'FUSED_LASSO'; % Combination of L1 and TV (3D TV)
method = 'FUSED_LASSO2D'; % 2D variant of FUSED_LASSO (2D TV)

% The following     parameters are defined in the preceding mat files
% holo                 Hologram, may be scaled to better size for computation
% pixel_num         Number of pixels (on one side) of the original hologram
%                   prior to resizing
% pps     Size in microns of camera pixels (camera property)
% lambda            Wavelength (um) of the laser used to form the hologram
% offsetZ           Location of first Z plane (um)
% deltaZ            Spacing between Z planes (um)
% nz                Number of Z planes
% pad_size          Size of zero padding buffer to used (if any)
% TV_subproblem_its Number of iterations to use in TV proximal
% mu                L1 regularization parameter
% tv_mu             Total variation regularization parameter
% num_its           Number of iterations to perform
% force_L           Used to identify whether to force the fasta step size
% fista_L           Stepsize value to use when force_L is true
% display_planes    List of plane indices to display after processing
% z_list            List of plane locations to reconstruct (um). Will
%                   override deltaZ, offsetZ, and nz parameters

img = im2double(imread([ './data/holo_t001_8bit_crop975_1045_1024.tif']));
bg = im2double(imread([ './data/background_8bit_crop975_1045_1024.tif']));
[pixel_num,~] = size(img);
pps = 0.07;

img = imresize(img, 0.5);
bg = imresize(bg, 0.5);
holo = (img - bg) ./ sqrt(bg);
% holo = imresize(holo, 0.5);

clear img bg;
lambda = 0.447;
offsetZ = 0;
deltaZ = 0.14*2;
nz = 150;
pad_size = 0;
fista_L = 0;
TV_subproblem_its = 5;
force_L = false;

holo = imresize(holo, 0.5);

% num_its = 10;
% mu = 1;
% % mu = 1.5;
% tv_mu = 1;
% % tv_mu = 10000;
% tv_mu = 0.1;
% num_its = 50;

% % For L1 regularization
% mu = 1.5;
% num_its = 50;

% % For TV regularization
% tv_mu = 0.3;
% num_its = 500;

% % For FL regularization
% mu = 1;
% tv_mu = 100;
% num_its = 100;

% mu = 1.3;
% tv_mu = 100;
% num_its = 10;

% BG division data L1
% mu = 0.28;
% num_its = 10;

% BG division data FL
mu = 0.05;
tv_mu = 0.1;
num_its = 10;

% mu = 0.2;
% tv_mu = 0.20;
% num_its = 10;
% 
% mu = 0.15;
% tv_mu = 20;
% num_its = 10;
% mu = 0.2;

%%

% Optional boolean flags for viewing results/steps
opts.plot_steps = true;
view_standard_reconstruction = false;
any_plots = true;

%% Prepare for execution (derived parameters)

% Scaling is to match with older version
% mu = mu/sqrt(nz);
% tv_mu = tv_mu/sqrt(nz);

% Determine size of pixels in actual input image
shrinkage_factor=pixel_num/size(holo,1);
sensor_size=pixel_num*pps;
deltaX=pps*shrinkage_factor;
deltaY=pps*shrinkage_factor;

if any_plots
    figure;imagesc(holo);title('CapturedData');axis image;
end

% Zero pad image
% Note: zero padding helps with fft to prevent image data from wrapping
% around the edges. Can also enable detection of objects outside the
% original FOV
holo=padarray(holo,[pad_size pad_size]);
range=pad_size*2+pixel_num;

[nx, ny]=size(holo);

% recon_params are hologram reconstruction parameters
recon_params.nx = nx;
recon_params.ny = ny;
recon_params.nz = nz;
recon_params.dz = deltaZ;
recon_params.offsetZ = offsetZ;
recon_params.resolution = pps*shrinkage_factor;
recon_params.wavelength = lambda;

% If z_list is given by user, override defaults
if isempty(z_list)
    recon_params.z_list = offsetZ:deltaZ:(nz-1)*deltaZ;
else
    recon_params.z_list = z_list;
    nz = length(z_list);
    recon_params.nz = nz;
end

% Define forward and adjoint operators for fasta
% The forward operator goes from a volume to a plane (object field to
% hologram). The adjoint operator is the opposite (standard reconstruction
% of hologram to get volume)
A = @(volume) Prop3D(volume, recon_params);
AT = @(plane) invProp3D(plane, recon_params);

% opts are parameters required for the fasta implementation of fista
opts.tol = 1e-8;  % Use super strict tolerance
opts.recordObjective = true; %  Record the objective function so we can plot it
opts.verbose=2;
opts.stringHeader='';
opts.accelerate = true;
opts.adaptive = false;
opts.stopRule = 'iterations';
opts.maxIters = num_its;
opts.TV_subproblem_its = TV_subproblem_its;
% opts.backtrack = false;
if force_L
    opts.L = fista_L; % To elliminate randomness with GPU
end

% Initial guess is all zeros
x0 = zeros(ny, nx, length(recon_params.z_list));

%% Standard reconstruction for visualization
if view_standard_reconstruction
    f_reconstruct = AT(holo);
    
    intensity = abs(f_reconstruct);
    cmbxy = max(intensity, [], 3);
    cmbxz = max(intensity, [], 1);
    cmbxz = rot90(flipud(squeeze(cmbxz)),-1);
    
    figure; imagesc(cmbxy);
    title('XY'); xlabel('X'); ylabel('Y'); axis image; axis equal;
    figure; imagesc(cmbxz);
    title('XZ'); xlabel('X'); ylabel('Z'); axis image; axis equal;
    
    cmax = max(intensity(:));
    cmin = min(intensity(:));
    figure;
    for zid = 1:length(recon_params.z_list)
        imagesc(intensity(:,:,zid), [cmin, cmax]);
        colorbar();
        axis image; axis equal;
        title(sprintf('Plane %d, z = %f um', zid, recon_params.z_list(zid)));
        drawnow();
        pause(0.1);
    end
end

%% Run the inverse reconstruction
tic;
switch method
    case 'L1'
        [f_reconstruct, outs_accel] = fasta_sparseLeastSquares(A,AT,holo,mu,x0, opts);
        tv_mu = 0;
    case 'TV'
        [f_reconstruct, outs_accel] = fasta_fusedLasso(A,AT,holo,0,tv_mu,x0, opts);
        mu = 0;
    case 'TV2D'
        [f_reconstruct, outs_accel] = fasta_fusedLasso2D(A,AT,holo,0,tv_mu,x0, opts);
        mu = 0;
    case 'FUSED_LASSO'
        [f_reconstruct, outs_accel] = fasta_fusedLasso(A,AT,holo,mu,tv_mu,x0, opts);
    case 'FUSED_LASSO2D'
        [f_reconstruct, outs_accel] = fasta_fusedLasso2D(A,AT,holo,mu,tv_mu,x0, opts);
end
toc;

save Outputs/result_fasta.mat f_reconstruct outs_accel recon_params opts;

%% Display results
% These are sample outputs, best tuned to your own needs for your data

fprintf('backtracks = %d, L = %f, initialStepsize = %f\n',...
    outs_accel.backtracks, outs_accel.L, outs_accel.initialStepsize);
fprintf('Sparsity = %f\n', nnz(f_reconstruct(:)) / numel(f_reconstruct(:)));

if any_plots
    figure; semilogy(outs_accel.residuals); title('residuals');
    % figure; semilogy(outs_accel.stepsizes); title('stepsizes');
    % figure; semilogy(outs_accel.normalizedResiduals); title('normalizedResiduals');
    % figure; semilogy(outs_accel.objective); title('objective');
    % figure; semilogy(outs_accel.funcValues); title('funcValues');
    figure; semilogy(outs_accel.fVals); title('fvals (Objective function)');
    
    figure;
    hold all;
    semilogy(outs_accel.error_norm);
    semilogy(outs_accel.L1_norm);
    semilogy(outs_accel.TV_norm);
    title('Norms');
    legend('||Ax - b||_2^2', '||x||_1', 'TV(x)');
    set(gca, 'YScale', 'log');
    
    % Multiply norms by regularization
    L1_norm_reg = outs_accel.L1_norm*mu;
    TV_norm_reg = outs_accel.TV_norm*tv_mu;
    objective_norm = outs_accel.error_norm + L1_norm_reg + TV_norm_reg;
    figure;
    hold all;
    semilogy(objective_norm);
    semilogy(outs_accel.error_norm);
    if ~any(strcmp(method,{'TV','TV2D'}))
        semilogy(L1_norm_reg);
    end
    if ~strcmp(method, 'L1')
        semilogy(TV_norm_reg);
    end
    title('Scaled Norms');
    if any(strcmp(method, {'FUSED_LASSO','FUSED_LASSO2D'}))
        legend('Objective', '||Ax - b||_2^2', '\lambda_L_1*||x||_1', '\lambda_T_V*TV(x)');
    elseif any(strcmp(method,{'TV','TV2D'}))
        legend('Objective', '||Ax - b||_2^2', '\lambda_T_V*TV(x)');
    else
        legend('Objective', '||Ax - b||_2^2', '\lambda_L_1*||x||_1');
    end
    
    set(gca, 'YScale', 'log');
end

intensity = abs(f_reconstruct);
minI = min(intensity(:));
maxI = max(intensity(:));
intensity = (intensity - minI) / (maxI - minI);
cmbxy = max(intensity, [], 3);
cmbxz = max(intensity, [], 1);
cmbxz = rot90(flipud(squeeze(cmbxz)),-1);
imwrite(cmbxy, './Outputs/cmbxy.tif');
imwrite(cmbxz, './Outputs/cmbxz.tif');
% for z = 1:size(intensity, 3)
%     imwrite(intensity(:,:,z), sprintf('./Outputs/Planes/plane_%03d.tif', z));
% end
if any_plots
    plotVolumeImage(intensity);
end

props = regionprops(intensity > 0, 'Centroid', 'Area');
num_props = length(props);
fprintf('Number of objects segmented: %d\n', num_props);

if any_plots
%     % View projections
%     figure; imagesc(cmbxy); title('XY'); xlabel('X'); ylabel('Y');
%     axis image; colormap('gray'); axis equal;
%     figure; imagesc(cmbxz); title('XZ'); xlabel('X'); ylabel('Z');
%     axis image; colormap('gray'); axis equal;
    
%     % View thresholded projections
%     figure; imagesc(cmbxy>0); title('Binarized XY'); xlabel('X'); ylabel('Y');
%     axis image; colormap('gray'); axis equal;
%     figure; imagesc(cmbxz>0); title('Binarized XZ'); xlabel('X'); ylabel('Z');
%     axis image; colormap('gray'); axis equal;
    
    % Display some selected planes
    if ~isempty(display_planes)
        for i = 1:length(display_planes)
            figure; imagesc(intensity(:,:,display_planes(i)));
            title(sprintf('Plane %d', display_planes(i)));
            axis image; colormap('gray'); xlabel('X'); ylabel('Y');
        end
    end
    
    figure;
    xpix = 0:deltaX*(size(cmbxz,2)-1);
    ypix = recon_params.z_list;
    imagesc(xpix, ypix, cmbxz);
    axis image;
    axis equal;
    xlabel('X (um)');
    ylabel('Z (um)');
end

% Write csv file
% nonzero_list = find(f_reconstruct ~= 0);
% [y,x,z] = ind2sub([ny,nx,nz], nonzero_list);
% value = f_reconstruct(nonzero_list);
% csvwrite('Outputs/inverse_voxels_0001.csv', [x,y,z, value]);

if any_plots
    est = A(f_reconstruct);
    figure; imagesc(est); axis image; title('Estimated'); colorbar();
    resid = holo - est;
    figure; imagesc(resid); axis image; title('Residual'); colorbar();
end