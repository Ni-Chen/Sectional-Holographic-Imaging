function [Phase Pupil]=MyMakingPhase(Nx,Ny,z,lambda,deltaX,deltaY,NA,mag,ref_index)

if(size(NA,2)==1)
    NA(2)=NA;
end;

k=1/(lambda*mag)*ref_index;

X=[ceil(-Nx/2):1:ceil(Nx/2-1)]'.*(1/(Nx*deltaX));
Y=[ceil(-Ny/2):1:ceil(Ny/2-1)].*(1/(Ny*deltaY));

kx=repmat(X,1,Ny);
ky=repmat(Y,Nx,1);
kp=sqrt(kx.^2+ky.^2);

term=k.^2-kp.^2;
term(term<0)=0;


Phase=exp(j*2*pi*z*(mag.^1)*sqrt(term));


Pupil=zeros(Nx,Ny);
Pupil(kx.^2./(k*NA(1)/ref_index).^2+ky.^2./(k*NA(2)/ref_index).^2<=1)=1;
% Pupil(kx.^2./(k*NA(1)/(ref_index*mag)).^2+ky.^2./(k*NA(2)/(ref_index*mag)).^2<=1)=1;
% Pupil(kp<=k*NA/ref_index)=1;
% Pupil((abs(kx)<=k*NA/shrinkage_factor)&(abs(ky)<=k*NA/shrinkage_factor))=1;
% Pupil=ones(Nx,Ny);