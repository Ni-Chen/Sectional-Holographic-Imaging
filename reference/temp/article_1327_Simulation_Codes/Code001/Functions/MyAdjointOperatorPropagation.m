function eta=MyAdjointOperatorPropagation(S,E,Nx,Ny,Nz,phase3D,pupil)

S=reshape(MyV2C(S),Nx,Ny);

eta=MyAdjointPropagation(S,E,Nx,Ny,Nz,phase3D,pupil);

eta=MyC2V(eta(:));
