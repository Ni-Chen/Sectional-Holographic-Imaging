function g=MyFunc4DeblurringImageSpace(f,Q)

[nx,ny,nzi,nzo]=size(Q);

f=reshape(f,nx,ny,nzi);
g=zeros(nx,ny,nzo);


for layer=1:nzi
    f(:,:,layer)=fft2(f(:,:,layer));
end;


for out_layer=1:nzo
    %         layer=abs(in_layer-out_layer)+1;
    g(:,:,out_layer)=sum(f.*Q(:,:,:,out_layer),3);
end;

for layer=1:nzo
    g(:,:,layer)=ifft2(g(:,:,layer));
end;


g=g(:);