function eta=MyAdjointOperatorPropagation(S,E,Nx,Ny,Nz,phase3D,pupil,propW)

if (nargin < 8), propW = sqrt(Nz);  end

S=reshape(MyV2C(S),Nx,Ny);

eta=MyAdjointPropagation(S,E,Nx,Ny,Nz,phase3D,pupil,propW);

eta=MyC2V(eta(:));
