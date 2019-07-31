function S=FresnelScalingAngularSpectrumFor(obj,resol,Nx,Ny,zo,lam,mask)

obj=reshape(MyV2C(obj),Ny,Nx);

tobj=fftshift(fft2(obj));
% figure(21);imagesc([], [],abs(tobj()));colormap;axis image;title('FT of patch object');colorbar;

[X,Y]=meshgrid(-size(tobj,2)/2:1:size(tobj,2)/2-1, size(tobj,1)/2-1:-1:-size(tobj,1)/2+0);
fx=X./(Nx.*resol);
fy=Y./(Ny.*resol);
% clear X Y;

k0=2*pi/lam;
kx=2*pi.*fx;
ky=2*pi.*fy;

tR=exp(1i.*(sqrt(k0.^2-kx.^2-ky.^2)-k0+(kx.^2+ky.^2)./(2*k0)));
% figure(22);imagesc([], [],angle(tR()));colormap;axis image;title('Corrector');colorbar;

Quad=exp(1i.*k0*resol.^2/(2*zo).*((resol.*X).^2+(resol.*Y).^2));
% figure(23);imagesc([], [],angle(Quad()));colormap;axis image;title('Fresnel quad');colorbar;

Phase=exp(1i.*k0.*zo+1i.*k0./(2*zo).*(fx.^2+fy.^2));

S=mask.*Phase.*-1i.*k0/zo.*fftshift(fft2(ifft2(ifftshift(tobj.*tR)).*Quad));
% figure(2031);imagesc([], [],abs(S()));colormap;axis
% image;title('Propagated patch object');colorbar;

S=MyC2V(S(:));