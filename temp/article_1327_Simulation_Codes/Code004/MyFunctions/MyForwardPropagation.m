function S=MyForwardPropagation(eta,E,Nx,Ny,Nz,phase3D,pupil,propW,aperture)

if (nargin < 9), propW = sqrt(Nz);  end

cEs=zeros(Nx,Ny,Nz);
Es=eta.*E;

for i=1:Nz
    cEs(:,:,i)=fftshift(fft2(Es(:,:,i)));
end

cEsp=sum(cEs.*phase3D.*pupil,3)./propW;

S=(ifft2(ifftshift(cEsp))).*aperture;