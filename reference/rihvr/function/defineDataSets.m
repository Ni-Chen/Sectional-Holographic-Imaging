function defineDataSets(varargin)

if size(varargin) < 1
    main_pathn = './';
else
    main_pathn = varargin{1};
end

%% Define Brady's dataset
% Contains g, the image data
load([main_pathn, 'Data/dandelion.mat']);
g = imresize(g, 0.5);

% number of detector pixels
pixel_num=1024;

% size of detector pixels (um)
detector_size=5.2;

% wavelength (um)
lambda=0.633;

% distance between each axial plane (um)
deltaZ=8000;

% distance from detector to first reconstructed plane (um)
offsetZ=0;

% number of axial planes
nz=10;

% number of zeros to pad matrix by in each direction
% pad_size=100;
pad_size = 0;

mu = 0.001;
tv_mu = 0.05;
num_its = 100;

force_L = false;
L = 2.230210;

% test parameters
num_its = 10;
TV_subproblem_its = 5;
mu = 0.1;
tv_mu = 0.5;

mu = 0.10;
tv_mu = 0.01;
num_its = 50;

% L1 norm
mu = 0.15;
num_its = 50;

% TV, 100 planes
mu = 0.15;
tv_mu = 0.05;
num_its = 10;

% TV, 10 planes
tv_mu = 0.01;
num_its = 500;

% Fuse lasso, 10 planes
tv_mu = 0.02;
mu = 0.01;
num_its = 500;

ratio = 2;
mu = 0.05;
tv_mu = ratio * mu;

% nz = nz * 10;
% deltaZ = deltaZ / 10;
display_planes = nz*[2/10,7/10] + 1;

save([main_pathn, 'Data/Brady_dandelion_parameters.mat']);
clearvars -except main_pathn


%% Define Microchannel images (smallest size)
% img = double(imread('../Images/MC_128_Holo_00001.tif'));
% bg = double(imread('../Images/MC_128_BG_00001.tif'));
img = im2double(imread([main_pathn, 'Data/MC_128_Holo_00001.tif']));
bg = im2double(imread([main_pathn, 'Data/MC_128_BG_00001.tif']));
g = (img - bg) ./ sqrt(bg);
pixel_num = 256;
detector_size = 1.1;
lambda = 0.632;
offsetZ = 0;
deltaZ = 4;
nz = 256;
pad_size = 0;
mu = 0.2;
fista_L = 0;
num_its = 20;
TV_subproblem_its = 5;
force_L = false;
L = 11.405158;

% Trial parameters
% mu = 0.25;
% tv_mu = 0.5;
% tv_iters = 20;
% 
% mu = 0;
% tv_mu = 0.5;
% tv_iters = 10;
% num_its = 5;

save([main_pathn, 'Data/microchannel_params.mat']);
clearvars -except main_pathn

%% Nanowire in suspension

% img = im2double(imread([main_pathn, 'Data/nanowire_0001.tif']));
% bg = mean(img(:));

img = im2double(imread([main_pathn, 'Data/holo_t001_8bit_crop975_1045_1024.tif']));
bg = im2double(imread([main_pathn, 'Data/background_8bit_crop975_1045_1024.tif']));
pixel_num = 1024;
detector_size = 0.07;
img = imresize(img, 0.5);
bg = imresize(bg, 0.5);

g = (img - bg) ./ sqrt(bg);
% imwrite(imresize(img,0.5),[main_pathn, 'Data/holo_nanowire_shrunk.tif']);
% imwrite(imresize(bg,0.5),[main_pathn, 'Data/background_nanowire_shrunk.tif']);
clear img bg;
% pixel_num = 256;
% detector_size = 0.14;
lambda = 0.447;
offsetZ = 0;
deltaZ = 0.14;
nz = 300;
pad_size = 0;
mu = 0.2;
tv_mu = 0.1;
fista_L = 0;
num_its = 20;
TV_subproblem_its = 5;
force_L = false;

g = imresize(g, 0.5);
nz = nz/2;
deltaZ = deltaZ*2;

num_its = 10;
mu = 1;
% mu = 1.5;
tv_mu = 1;
% tv_mu = 10000;
tv_mu = 0.1;
num_its = 50;

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
mu = 0.28;
num_its = 10;

% BG division data FL
mu = 0.05;
tv_mu = 0.1;
num_its = 10;

mu = 0.2;
tv_mu = 0.20;
num_its = 10;

mu = 0.15;
tv_mu = 20;
num_its = 10;
mu = 0.2;

save([main_pathn, 'Data/nanowire_params.mat']);
clearvars -except main_pathn

%% Underwater hologram from Madison Lake June 2016
this_pathn = [main_pathn, 'Data/Underwater_Madison_6-16/'];
img = im2double(imread([this_pathn, 'image_1756.png']));
bg = im2double(imread([this_pathn, 'background.tif']));

img = mean(img, 3);

g = (img - bg) ./ sqrt(bg);
clear img bg;
pixel_num = 1024;
detector_size = 2.2;
lambda = 0.632;
offsetZ = 0;
deltaZ = 70;
nz = 300;
pad_size = 0;
fista_L = 0;
TV_subproblem_its = 5;
force_L = false;

mu = 0.45;
tv_mu = 5;
num_its = 10;

mu = 0.6;

g = g(513:1024, 1:512);
pixel_num = 512;
g = imresize(g, 0.25);

save([main_pathn, 'Data/underwater_madison_params.mat']);
clearvars -except main_pathn;

%% Synthetic image of letter alpha
img = im2double(imread('Data/holo_single_alpha_image.tif'));
bg = im2double(imread('Data/bg_single_alpha_image.tif'));

g = (img - bg) ./ sqrt(bg);
clear img bg;
pixel_num = 512;
detector_size = 1.1;
lambda = 0.632;
offsetZ = 0;
deltaZ = 5;
nz = 200;
pad_size = 0;
fista_L = 0;
TV_subproblem_its = 4;
force_L = false;

mu = 0;
tv_mu = 50;
num_its = 500;

% g = imresize(g, 0.25);

deltaZ = 250;
nz = 6;

display_planes = round(500/deltaZ) + 1;

save([main_pathn, 'Data/synthetic_alpha_params.mat']);
clearvars -except main_pathn;

%% Synthetic image of Endo phantom
img = im2double(imread('Data/endo_phantom_holo.tif'));
bg = im2double(imread('Data/endo_phantom_bg.tif'));

g = (img - bg) ./ sqrt(bg);

% load('Data/endo_phantom_complex_holo.mat');
% g = holo;

clear img bg;
pixel_num = 1024;
detector_size = 9;
lambda = 0.632;
offsetZ = 0;
deltaZ = 2000;
nz = 100;
pad_size = 0;
fista_L = 0;
TV_subproblem_its = 10;
force_L = false;

% From Endo paper
mu = 0.2;
tv_mu = 0.01;
num_its = 1000;

mu = 1;
tv_mu = 0.01;
num_its = 1000;

% g = imresize(g, 0.25);

display_planes = [100, 102, 104, 106]*1000;
display_planes = round(display_planes/deltaZ) + 1;

z_list = [0, 100, 102, 104, 106]*1000;
display_planes = [1,2,3,4,5];

save([main_pathn, 'Data/endo_phantom_large_params.mat']);
clearvars -except main_pathn;

%% Dunaliella
this_pathn = [main_pathn, 'Data/Dunaliella/'];
img = im2double(imread([this_pathn, 'CoreView_284_Camera_0001.tif']));
bg = im2double(imread([this_pathn, 'background.tif']));
bg_med = im2double(imread([this_pathn, 'background_MED.tif']));
img = img(1:1024, 1:1024);
bg = bg(1:1024, 1:1024);

g = (img - bg) ./ sqrt(bg);
pixel_num = 1024;
detector_size = 1.1;
lambda = 0.632;
offsetZ = 0;
deltaZ = 12;
nz = 250;
pad_size = 0;
fista_L = 0;
TV_subproblem_its = 5;
force_L = false;

mu = 0.3;
tv_mu = 0.05;
num_its = 30;

save([main_pathn, 'Data/dunaliella_params.mat']);
clearvars -except main_pathn;


%% Various unsaved

clearvars();
% Orcas Island
input_pathn = 'H:/Lab Members/Current Members/Graduate Students/Kevin Mallery/Tool Development/TestFiles/Orcas_Island/HighSpeed_LowMagnification/';
img = im2double(imread([input_pathn, 'w_7m_1000001.tif']));
bg = im2double(imread([input_pathn, 'background.tif']));
g = (img - bg)./sqrt(bg);
pixel_num = 1024;
detector_size = 1.71;
lambda = 0.632;
offsetZ = 0;
deltaZ = 400;
nz = 20;
pad_size = 0;
TV_subproblem_its = 5;
mu = 0.07;
tv_mu = 0.1;
num_its = 100;
force_L = false;
fista_L = 1000;
g = imresize(g, 0.25);

clearvars();