function [Phase3D Pupil]=MyMakingPhase3D(Nx,Ny,Nz,lambda,...
    deltaX,deltaY,Z,NA,mag,ref_index,offsetZ)

if (nargin < 8), NA = 1;  end
if (nargin < 9), mag = 1;  end
if (nargin < 10), ref_index = 1;  end
if (nargin < 11), offsetZ = 0;  end

if(double(size(Z)==[1 1])*double(Nz>1))
    Z=[0:Nz-1].*Z+offsetZ;
else
    Z=Z+offsetZ;
end;


Phase3D=zeros(Nx,Ny,Nz);
Pupil=zeros(Nx,Ny,Nz);

for i=1:Nz

    [Phase3D(:,:,i) Pupil(:,:,i)]=...
        MyMakingPhase(Nx,Ny,Z(i),lambda,deltaX,deltaY,NA,mag,ref_index);
end;
