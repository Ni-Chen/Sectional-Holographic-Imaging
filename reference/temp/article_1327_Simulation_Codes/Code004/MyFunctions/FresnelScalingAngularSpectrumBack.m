function eE=FresnelScalingAngularSpectrumBack(S,resol,Nx,Ny,zo,lam,mask)

S=reshape(MyV2C(S),Ny,Nx);
[X,Y]=meshgrid(-size(S,2)/2:1:size(S,2)/2-1, size(S,1)/2-1:-1:-size(S,1)/2+0);
fx=X./(Nx.*resol);
fy=Y./(Ny.*resol);
% clear X Y;
% fx=fx+1/(Nx.*dx)*px*ones(size(fx));
% fy=fy+1/(Ny.*dx)*py*ones(size(fy));

k0=2*pi/lam;
kx=2*pi.*fx;
ky=2*pi.*fy;

tR=exp(1i.*(sqrt(k0.^2-kx.^2-ky.^2)-k0+(kx.^2+ky.^2)./(2*k0)));
% figure(22);imagesc([], [],angle(tR()));colormap;axis image;title('Corrector');colorbar;

Quad=exp(1i.*k0*resol.^2/(2*zo).*((resol.*X).^2+(resol.*Y).^2));
% figure(23);imagesc([], [],angle(Quad()));colormap;axis
% image;title('Fresnel quad');colorbar;

Phase=exp(1i.*k0.*zo+1i.*k0./(2*zo).*(fx.^2+fy.^2));

eE=ifft2(ifftshift(fftshift(fft2(ifft2(ifftshift(S.*(mask).*conj(Phase).*-zo/(1i.*k0))).*conj(Quad))).*conj(tR)));
% figure(2041);imagesc([], [],abs(eE()));colormap;axis image;title('Backward propagation');colorbar;

eE=MyC2V(eE(:));