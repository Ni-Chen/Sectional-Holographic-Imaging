function S=MyForwardOperatorPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil,propW,aperture)

if (nargin < 9), propW = sqrt(Nz);  end

eta=reshape(MyV2C(eta),Nx,Ny,Nz);

S=MyForwardPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil,propW,aperture);

S=MyC2V(S(:));
