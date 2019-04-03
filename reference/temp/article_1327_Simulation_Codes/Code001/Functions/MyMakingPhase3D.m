function [Phase3D Pupil]=MyMakingPhase3D(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,deltaZ,offsetZ,sensor_size)

Z=[0:Nz-1].*deltaZ+offsetZ;

s=sensor_size/2;

Phase3D=zeros(Nx,Ny,Nz);
Pupil=zeros(Nx,Ny,Nz);

for i=1:Nz
    NA=s/sqrt(Z(i).^2+s.^2);
    [Phase3D(:,:,i) Pupil(:,:,i)]=...
        MyMakingPhase(Nx,Ny,Z(i),lambda,deltaX,deltaY,NA);
end;
