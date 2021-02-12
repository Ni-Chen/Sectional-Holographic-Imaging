function Q=MyQ4DeblurringImageSpace(func,nx,ny,nz)

Q=zeros(nx,ny,nz,nz);


for in_layer=1:nz
    delta=zeros(nx,ny,nz);
    delta(1,1,in_layer)=1;
    
    psf=(abs(MyV2C(func(MyC2V(delta(:)))))).^2;
    psf=reshape(psf,nx,ny,nz);
    
    for out_layer=1:nz
%         Q(:,:,in_layer,out_layer)=(psf(:,:,out_layer));
        Q(:,:,in_layer,out_layer)=fft2(psf(:,:,out_layer));
    end;
end;