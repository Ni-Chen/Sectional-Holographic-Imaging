function g=MyFastFunc4DeblurringImageSpace(f,nx,ny,nz)


f=reshape(f,nx,ny,nz);

g=zeros(nx,ny,nz);

s=sum(squeeze(sum(f,1)),1)/(nx*ny*nz);


for z1=1:nz;
    for z2=1:nz;
        if(z1==z2)
            g(:,:,z1)=g(:,:,z1)+f(:,:,z2)-s(z2)*(nz-1);
        else
            g(:,:,z1)=g(:,:,z1)+s(z2);
        end;
    end;
end;

g=g(:);