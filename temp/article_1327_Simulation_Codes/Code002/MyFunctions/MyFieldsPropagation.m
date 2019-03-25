function E=MyFieldsPropagation(E0,Nx,Ny,Nz,phase3D,pupil)

E=zeros(Nx,Ny,Nz);

for i=1:Nz
    cE0=fftshift(fft2(E0));

    cE=cE0.*conj(phase3D(:,:,i)).*pupil(:,:,i);

    E(:,:,i)=ifft2(ifftshift(cE));
end