function eta=MyAdjointOperatorPropagation(S,E,Nx,Ny,Nz,phase3D,pupil,aperture,propW)

if (nargin < 9), propW = sqrt(Nz);  end

S=reshape(MyV2C(S),Nx,Ny);

eta=MyAdjointPropagation(S,E,Nx,Ny,Nz,phase3D,pupil,aperture,propW);

eta=MyC2V(eta(:));
