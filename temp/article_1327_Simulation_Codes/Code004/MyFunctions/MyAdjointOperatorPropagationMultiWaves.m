function eta=MyAdjointOperatorPropagationMultiWaves(S,E,Nx,Ny,Nz,Nl,phase3D,pupil,aperture,propW)

if (nargin < 10), propW = sqrt(Nz);  end

eta=zeros(Nx*Ny*Nz,1);
S=reshape(MyV2C(S),Nx,Ny,Nl);

for i=1:Nl
    s=S(:,:,i);
    s=MyC2V(s(:));
    eta=eta+MyV2C(MyAdjointOperatorPropagation(s,E(:,:,:,i),Nx,Ny,Nz,...
        phase3D(:,:,:,i),pupil(:,:,:,i),aperture,propW))./sqrt(Nl);
end;

eta=MyC2V(eta(:));
