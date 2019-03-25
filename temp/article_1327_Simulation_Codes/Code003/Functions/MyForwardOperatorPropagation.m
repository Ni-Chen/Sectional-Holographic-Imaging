function S=MyForwardOperatorPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil,Apert)

eta=reshape(MyV2C(eta),Nx,Ny,Nz);

S=MyForwardPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil,Apert);

S=MyC2V(S(:));
