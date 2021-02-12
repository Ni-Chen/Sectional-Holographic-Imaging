function g=MyDerivationTV2D(f)

[nx ny]=size(f);

f_d=circshift(f,[1 0]);
f_dl=circshift(f,[1 -1]);
f_u=circshift(f,[-1 0]);
f_l=circshift(f,[0 -1]);
f_r=circshift(f,[0 1]);
f_ur=circshift(f,[-1 1]);

A1=(f-f_d)./sqrt((abs(f-f_d)).^2+(abs(f_dl-f_d)).^2);

A2=(2*f-f_u-f_l)./sqrt((abs(f_u-f)).^2+(abs(f_l-f)).^2);

A3=(f-f_r)./sqrt((abs(f_ur-f_r)).^2+(abs(f-f_r)).^2);

if 0
    g=A1+A2+A3;
else
    g=zeros(nx,ny);
    g(1:nx-1,1:ny-1)=A1(1:nx-1,1:ny-1)+A2(1:nx-1,1:ny-1)+A3(1:nx-1,1:ny-1);
    g(nx,1:ny-1)=A1(nx,1:ny-1);
    g(1:nx-1,ny)=A3(1:nx-1,ny);

    if(((f(nx,ny)-f_d(nx,ny))*(f(nx,ny)-f_r(nx,ny)))>0)
        g(nx,ny)=2;
    end;

end;