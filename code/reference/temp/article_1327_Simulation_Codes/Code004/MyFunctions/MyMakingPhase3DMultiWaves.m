function [Phase3D Pupil]=MyMakingPhase3DMultiWaves(Nx,Ny,Nz,Nl,lambda,...
    deltaX,deltaY,Z,NA,mag,ref_index,offsetZ)

if (nargin < 9), NA = 1;  end
if (nargin < 10), mag = 1;  end
if (nargin < 11), ref_index = 1;  end
if (nargin < 12), offsetZ = 0;  end

Phase3D=zeros(Nx,Ny,Nz,Nl);
Pupil=zeros(Nx,Ny,Nz,Nl);

for i=1:Nl
    [Phase3D(:,:,:,i) Pupil(:,:,:,i)]=MyMakingPhase3D(Nx,Ny,Nz,lambda(i),...
    deltaX,deltaY,Z,NA,mag,ref_index,offsetZ);
end;
