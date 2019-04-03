close all;
clear all;clc;
addpath('./Functions');


%% parameters (1)
% data size
nx=128; ny=128; nz=1;
% wavelength (um)
lambda=0.633;
% pixels pitch (um)
detector_size=100;
% detector size (um)
sensor_size=nx*detector_size;
% distance between each axial plane (um)
deltaZ=2000*1000;
% distance from detector to first reconstructed plane (um)
offsetZ=deltaZ;
deltaX=detector_size;
deltaY=detector_size;
Nx=nx;
Ny=ny*nz*2;
Nz=1;



%% Object generation (2)
obj=zeros(nx,ny,nz);
obj=phantom(length(obj));
figure;imagesc(plotdatacube(obj,nz));title('target object');axis image;drawnow;
axis off;colormap(hot);colorbar;



%% Sparse aperture generation (3)
%random sparse aperture
sparse=randn(ny,nx);
threshold=1.0;   %1.1
sparse(sparse<threshold)=0;
sparse(sparse>=threshold)=1;
openning=sum(sparse(:))/(ny*nx)*100;
figure(); imagesc(sparse); axis image;title(['sparse aperture - ' num2str(openning)]);
axis off;colormap(hot);colorbar;



%% Propagation operator (4)
% Forward propagation
A = @(f_twist) FresnelScalingAngularSpectrumFor(f_twist,deltaX,nx,ny,offsetZ,lambda,sparse);
% Backward propagation
AT = @(g) FresnelScalingAngularSpectrumBack(g,deltaX,nx,ny,offsetZ,lambda,sparse);



%% Holographic measurement, back-propagation, and averaging (5)
num_measure = 30;
G_avg=0;
if 1
    for num=1:num_measure
        f_r=sqrt(obj/2).*randn(nx,ny,nz);
        f_i=sqrt(obj/2).*randn(nx,ny,nz);
        f=f_r+i.*f_i;
        G=reshape(MyV2C(A(MyC2V(f(:)))),nx,ny);
        G_avg=G_avg+G;
        figure(99);imagesc(plotdatacube(abs(G_avg),nz));axis image;colorbar;
        title(['number of realization:',num2str(num)]);drawnow;
    end;
    G_avg=G_avg./num_measure;   
end;
g=MyC2V(G_avg(:));
transf=AT(g);
transf=reshape(MyV2C(transf),nx,ny,nz);
figure;imagesc(plotdatacube(abs(transf),nz));title('Numerical backpropagation');axis image;drawnow;
axis off;colormap(hot);colorbar;



%% TwIST algorithm (6)
tau = 0.01;
piter = 4;
tolA = 1e-6;
iterations = 500; %1000

Psi = @(f,th) MyTVpsi(f,th,0.05,piter,Nx,Ny,Nz);
Phi = @(f) MyTVphi(f,Nx,Ny,Nz);

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
    'Verbose', 1);

f_reconstruct=reshape(MyV2C(f_reconstruct),nx,ny,nz);
figure;imagesc(plotdatacube(abs(f_reconstruct),nz));title('Compressive reconstruction');axis image;drawnow;
axis off;colormap(hot);colorbar;

