function E=MyFieldsPropagationMultiWaves(E0,Nx,Ny,Nz,Nl,phase3D,pupil)

E=zeros(Nx,Ny,Nz,Nl);

for i=1:Nl
    E(:,:,:,i)=MyFieldsPropagation(E0(:,:,i),Nx,Ny,Nz,phase3D(:,:,:,i),pupil(:,:,:,i));
end;