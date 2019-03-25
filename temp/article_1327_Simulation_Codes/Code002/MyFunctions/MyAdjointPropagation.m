function eta=MyAdjointPropagation(S,E,Nx,Ny,Nz,phase3D,pupil,propW)

if (nargin < 8), propW = sqrt(Nz);  end

cEsp=fftshift(conj(ifft2(conj(S))));

cEs=conj(phase3D).*conj(pupil).*repmat(cEsp./propW,[1 1 Nz]);

eta=zeros(Nx,Ny,Nz);
for i=1:Nz
    eta(:,:,i)=conj(fft2(conj(ifftshift(cEs(:,:,i)))));
end

eta=conj(E).*eta;