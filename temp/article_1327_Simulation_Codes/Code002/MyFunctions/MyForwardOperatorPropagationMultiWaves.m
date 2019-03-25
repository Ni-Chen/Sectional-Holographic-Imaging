function S=MyForwardOperatorPropagationMultiWaves(eta,E,Nx,Ny,Nz,Nl,phase3D,pupil,propW)

if (nargin < 9), propW = sqrt(Nz);  end

S=zeros(Nx*Ny,Nl);

for i=1:Nl
    S(:,i)=MyV2C(MyForwardOperatorPropagation(eta,E(:,:,:,i),Nx,Ny,Nz,...
        phase3D(:,:,:,i),pupil(:,:,:,i),propW))./sqrt(Nl);
end;

S=MyC2V(S(:));