close all;
clear all;clc;
addpath('./Functions');


%% Paprameters (1)
nx=64;  % data size
ny=64;
nz=5;
lambda=0.633;  % wavelength (um)
detector_size=30;  % pixel pitch (um)
sensor_size=nx*detector_size;  % detector size (um)
deltaZ=20*1000;  % axial spacing (um)
offsetZ=0*1000;  % distance from detector to first reconstructed plane (um)
% offsetZ=0*1000;  % distance from detector to first reconstructed plane (um)
offsetZ=50*1000;  % distance from detector to first reconstructed plane (um)

deltaX=detector_size;
deltaY=detector_size;
Nx=nx;
Ny=ny*nz*2;
Nz=1;


%% Object generation (2)
load('D48.mat');load('I48.mat');load('S48.mat');load('P48.mat');
% D48 = (imread('S.tif'))./255;
% I48 = (imread('N.tif'))./255;
% S48 = (imread('U.tif'))./255;
% P48 = (imread('E.tif'))./255;

f=zeros(nx,ny,nz);
f(:,:,2)=0.7*D48;
f(:,:,3)=0.5*I48;
f(:,:,4)=1.0*S48;
f(:,:,5)=0.9*P48;

figure;imagesc(plotdatacube(abs(f)));title('3D object');axis image;drawnow;
axis off; colormap(hot); colorbar;


%% Propagation kernel (3)
[Phase3D Pupil]=MyMakingPhase3D(nx,ny,nz,lambda,deltaX,deltaY,deltaZ,offsetZ,sensor_size);
figure;imagesc(plotdatacube(angle(Phase3D)));title('Phase of kernel');axis image;drawnow;
axis off; colormap(hot); colorbar;
E0=ones(nx,ny);  % illumination light
E0=E0.*exp(i.*pi*0.0);
E=MyFieldsPropagation(E0,nx,ny,nz,Phase3D,Pupil);  % propagation of illumination light


%% Field measurement and backpropagation (4)
cEs=zeros(nx,ny,nz);
Es=f.*E;
for i=1:nz
    cEs(:,:,i)=fftshift(fft2(Es(:,:,i)));
end
cEsp=sum(cEs.*Phase3D.*Pupil,3);
S=(ifft2(ifftshift(cEsp)));
% squared field
% g=S.*conj(S);
% sf = ifftshift(fft2(fftshift(g)));
% rf = ifftshift(fft2(fftshift(ones(ny,nx))));
% g = ifftshift(ifft2(fftshift(sf-rf)));
% g = real(g);

% diffracted field

g=S+conj(S)+f(:,:,2).^2+f(:,:,3).^2+f(:,:,4).^2+f(:,:,5).^2;
% S = S./max(max(abs(S)));
% g= S + conj(S) + sum(f.^2,3)/Nz;

g= S+f(:,:,2).^2+f(:,:,3).^2+f(:,:,4).^2+f(:,:,5).^2;
% g= S+conj(S)+f(:,:,2).^2+f(:,:,3).^2+f(:,:,4).^2+f(:,:,5).^2;



figure;imagesc(abs(g));title('Diffracted field');axis image;
axis off; colormap(hot); colorbar;
g=MyC2V(g(:));
transf=MyAdjointOperatorPropagation(g,E,nx,ny,nz,Phase3D,Pupil);
transf=reshape(MyV2C(transf),nx,ny,nz);
figure;imagesc(plotdatacube(abs(transf)));title('Numerical backpropagation');axis image;drawnow;
axis off; colormap(hot); colorbar;


%% Propagation operator (5)
A = @(f_twist) MyForwardOperatorPropagation(f_twist,E,nx,ny,nz,Phase3D,Pupil);  % forward propagation operator
AT = @(g) MyAdjointOperatorPropagation(g,E,nx,ny,nz,Phase3D,Pupil);  % backward propagation operator


%% TwIST algorithm (6)
% twist parameters
tau = 0.01; 
piter = 4;
tolA = 1e-6;
iterations = 200;

Psi = @(f,th) MyTVpsi(f,th,0.05,piter,Nx,Ny,Nz);
Phi = @(f) MyTVphi(f,Nx,Ny,Nz);
vol3d_vec = C2V(f(:));
[f_reconstruct,dummy,obj_twist,...
    times_twist,dummy,mse_twist]= ...
    TwIST(g,A,tau,...
    'AT', AT, ...
    'Psi', Psi, ...
    'Phi',Phi, ...
    'Initialization',2,...
    'Monotone',1,...
    'StopCriterion',1,...
    'MaxIterA',iterations,...
    'MinIterA',iterations,...
    'ToleranceA',tolA,...
    'True_x', vol3d_vec, ...
    'Verbose', 1);



f_reconstruct=reshape(MyV2C(f_reconstruct),nx,ny,nz);
figure;imagesc(plotdatacube(abs(f_reconstruct)));title('Compressive reconstruction');axis image;drawnow;
axis off; colormap(hot); colorbar;

figure; semilogy(times_twist, obj_twist, 'LineWidth',2); ylabel('objective distance'); xlabel('CPU time (sec)');
figure; plot(times_twist, mse_twist, 'LineWidth',2); ylabel('MSE'); xlabel('CPU time (sec)');
