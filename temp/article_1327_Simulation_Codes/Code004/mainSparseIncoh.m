close all;
clear all;clc;
addpath('./MyFunctions');
% break

rand('state',1);
randn('state',1);


%% Parameters (1)
nx=128; ny=128; nz=1;
% wavelength (um)
lambda=[0.633];
% pixels pitch (um)
detector_size=100;
% detector size (um)
sensor_size=nx*detector_size;
% distance from detector to first reconstructed plane (um)
Z=[2000*1000];
deltaX=detector_size;
deltaY=detector_size;
% number of measurements
num_measure=30;
% Parameter for Tikhonov preconditioning
alpha=0.00001;



%% Object generation (2)
obj=zeros(nx,ny,nz)+0.1;
obj=phantom(length(obj));
figure;imagesc(plotdatacube(obj,nz));axis image;
title(['Intensity of object']);drawnow;
axis off;colormap(hot);colorbar;



%% Sparse aperture generation (3)
%random sparse aperture
sparse=randn(nx,ny);
threshold=1.0;   %1.1 1.0
sparse(sparse<threshold)=0;
sparse(sparse>=threshold)=1;
openning=sum(sparse(:))/(ny*nx)*100;
figure(); imagesc(sparse); axis image;title(['sparse aperture - ' num2str(openning)]);
axis off;colormap(hot);colorbar;




%% Propagation operators (4)
% Forward propagation
A = @(f_twist) FresnelScalingAngularSpectrumFor(f_twist,deltaX,nx,ny,Z,lambda,sparse);
% Backward propagation
AT = @(g) FresnelScalingAngularSpectrumBack(g,deltaX,nx,ny,Z,lambda,sparse);



%% Holographic Mmeasurement, back-propagation, and averaging (5)
average=zeros(nx*ny*nz,1);
if 1
    for num=1:num_measure
        f_r=sqrt(obj/2).*randn(nx,ny,nz);
        f_i=sqrt(obj/2).*randn(nx,ny,nz);
        f=f_r+i.*f_i;
        G=reshape(MyV2C(A(MyC2V(f(:)))),nx,ny);
%         G=G+r*(conj(G)+abs(G).^2);
%         G=awgn(G,snr,'measured');
        average=average+(abs(MyV2C(AT(MyC2V(G(:)))))).^2;
        figure(99);imagesc(plotdatacube(reshape(average./num,nx,ny,nz),nz));axis image;colorbar;
        title(['number of realization:',num2str(num)]);drawnow;
    end;
    average=average./num_measure;
    f_average=reshape(average(:),nx,ny,nz);
end;
figure;imagesc(plotdatacube(f_average,nz));axis image;
title(['Average of back-propagations where number of realization is ', num2str(num)]);drawnow;
axis off;colormap(hot);colorbar;



%% Incoherent system matrix (6)
newH=@(delta) AT(A(delta));
Q=MyQ4DeblurringImageSpace(newH,nx,ny,nz);
[pinvQ combQ]=MyCombinedPinvQ4DeblurringImageSpace(Q,alpha);



%% TwIST algorithm (7)
B=@(f_twist) MyFunc4DeblurringImageSpace(f_twist,Q);
C=@(f_twist) MyFunc4DeblurringImageSpace(f_twist,pinvQ);
D=@(f_twist) MyFunc4DeblurringImageSpace(f_twist,combQ);                % forward projection
DT=@(g) MyFunc4DeblurringImageSpace(g,permute(conj(combQ),[1 2 4 3]));  % backward projection

if 1
    tau = 0.4;  %0.2
    piter = 4;
    tolA = 1e-6;
    iterations = 500;

    Psi = @(f,th) MyTVpsi(f,th,0.05,piter,nx,ny,nz,[1 1 0]);
    Phi = @(f) MyTVphi(f,nx,ny,nz,[1 1 0]);

    [f_reconstruct,dummy,obj_twist,...
        times_twist,dummy,mse_twist]= ...
        MyTwIST(C(f_average(:)),D,tau,nx,ny,nz,...
        'AT', DT, ...
        'Psi', Psi, ...
        'Phi',Phi, ...
        'Initialization',2,...
        'Monotone',1,...
        'StopCriterion',1,...
        'MaxIterA',iterations,...
        'MinIterA',iterations,...
        'ToleranceA',tolA,...
        'Verbose', 1);

    f_reconstruct=reshape(f_reconstruct,nx,ny,nz);
    figure;imagesc(plotdatacube(f_reconstruct,nz));axis image;
    title({'TwIST reconstruction', 'where number of realizations is ', num2str(num)});
    axis off;colormap(hot);colorbar;
    drawnow;
 
end;

