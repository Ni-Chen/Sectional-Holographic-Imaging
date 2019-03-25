function eta=MyAdjointPropagation(S,E,Nx,Ny,Nz,phase3D,pupil,Apert)


% cEsp=fftshift(conj(ifft2(conj(real(S)))));
cEsp=fftshift(conj(ifft2(conj((S.*Apert)))));

cEs=conj(phase3D).*conj(pupil).*repmat(cEsp,[1 1 Nz]);

eta=zeros(Nx,Ny,Nz);
for i=1:Nz
    eta(:,:,i)=conj(fft2(conj(ifftshift(cEs(:,:,i)))));
end

eta=conj(E).*eta;