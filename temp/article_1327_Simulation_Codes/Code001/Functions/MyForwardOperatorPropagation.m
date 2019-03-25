function S=MyForwardOperatorPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil)

eta=reshape(MyV2C(eta),Nx,Ny,Nz);

S=MyForwardPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil);

S=MyC2V(S(:));
