function S=MyForwardOperatorPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil,propW)

if (nargin < 8), propW = sqrt(Nz);  end

eta=reshape(MyV2C(eta),Nx,Ny,Nz);

S=MyForwardPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil,propW);

S=MyC2V(S(:));
