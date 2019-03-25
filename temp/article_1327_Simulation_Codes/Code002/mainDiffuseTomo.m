close all;
clear all;clc;
addpath('./MyFunctions');

rand('state',1);
randn('state',1);



%% Parameters (1)
% data size
nx=64;
ny=64*2;
nz=3;
% wavelength (um)
lambda=[0.633];
% detector pixel pitch (um)
detector_size=50;
% detector size (um)
sensor_size=nx*detector_size;
% distance from detector to first reconstructed plane (um)
Z=[-50*1000,0,+50*1000];
deltaX=detector_size;
deltaY=detector_size;
mag=1;
ref_index=1;
NA=[lambda/(2*deltaX/mag),lambda/(2*deltaY/mag)];
% number of measurements
num_measure=30;
% Parameter for Tikhonov preconditioning
alpha=0.00001;



%% Object generation (2)
obj=zeros(nx,ny,nz)+0.1;
obj(5:60,32:42,1)=1;
obj(25:35,5:123,2)=1;
obj(5:60,96:106,3)=1;
figure;imagesc(plotdatacube(obj,nz));axis image;
title(['Intensity of object']);drawnow;
axis off;colormap(hot);colorbar;



%% Propagation kernel (3)
nl=size(lambda,2);
E0=ones(nx,ny);
E0=E0.*exp(1i.*pi*1*ones(nx,ny));
[Phase3D Pupil]=MyMakingPhase3DMultiWaves(nx,ny,nz,nl,lambda,deltaX,deltaY,Z,NA,mag,ref_index);
figure;imagesc(plotdatacube(angle(Phase3D.*1)));axis image;
title(['real(propagation kernel)']);
E=MyFieldsPropagationMultiWaves(E0,nx,ny,nz,nl,Phase3D,Pupil);
axis off;colormap(hot);colorbar;


%% Propagation operators (4)
% Forward propagation
A = @(f_twist) MyForwardOperatorPropagationMultiWaves(f_twist,E,nx,ny,nz,nl,Phase3D,Pupil);
% Backward propagation
AT = @(g) MyAdjointOperatorPropagationMultiWaves(g,E,nx,ny,nz,nl,Phase3D,Pupil);



%% Holographic measurement, back-propagation, and averaging (5)
average=zeros(nx*ny*nz,1);
if 1
    for num=1:num_measure
        f_r=sqrt(obj/2).*randn(nx,ny,nz);
        f_i=sqrt(obj/2).*randn(nx,ny,nz);
        f=f_r+i.*f_i;
        G=reshape(MyV2C(A(MyC2V(f(:)))),nx,ny);
        average=average+(abs(MyV2C(AT(MyC2V(G(:)))))).^2;
        figure(99);imagesc(plotdatacube(reshape(average./num,nx,ny,nz),nz));axis image;colorbar;
        axis off;colormap(hot);colorbar;
        title(['number of realization:',num2str(num)]);drawnow;
    end;
    average=average./num_measure;
    f_average=reshape(average(:),nx,ny,nz);
end;
figure;imagesc(plotdatacube(f_average,nz));axis image;
title(['Average of back-propagations where number of realization is ', num2str(num)]);drawnow;
axis off;colormap(hot);colorbar;
% break


%% Incoherent system matrix (6)
newH=@(delta) AT(A(delta));
Q=MyQ4DeblurringImageSpace(newH,nx,ny,nz);
[pinvQ combQ]=MyCombinedPinvQ4DeblurringImageSpace(Q,alpha);



%% TwIST algorithm (7)
B=@(f_twist) MyFunc4DeblurringImageSpace(f_twist,Q);
C=@(f_twist) MyFunc4DeblurringImageSpace(f_twist,pinvQ);
D=@(f_twist) MyFunc4DeblurringImageSpace(f_twist,combQ);
DT=@(g) MyFunc4DeblurringImageSpace(g,conj(combQ));

if 1
    tau = 1.5;
    piter = 4;
    tolA = 1e-6;
    iterations = 500;

    Psi = @(f,th) MyTVpsi(f,th,0.05,piter,nx,ny,nz,[1 1 0]);
    Phi = @(f) MyTVphi(f,nx,ny,nz,[1 1 0]);
    
    invf=reshape(C(f_average(:)),nx,ny,nz);
    figure;imagesc(plotdatacube(invf,nz));axis image;colorbar;
    title({'Tikhonov inverse of averaged backpropagations', 'where number of realizations is ', num2str(num)});
    axis off;colormap(hot);colorbar;
    drawnow;

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