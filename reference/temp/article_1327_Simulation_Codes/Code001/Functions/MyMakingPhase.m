function [Phase Pupil]=MyMakingPhase(Nx,Ny,z,lambda,deltaX,deltaY,NA)

k=1/lambda;

X=[ceil(-Nx/2):1:ceil(Nx/2-1)]'.*(1/(Nx*deltaX));
Y=[ceil(-Ny/2):1:ceil(Ny/2-1)].*(1/(Ny*deltaY));

kx=repmat(X,1,Ny);
ky=repmat(Y,Nx,1);
kp=sqrt(kx.^2+ky.^2);

term=k.^2-kp.^2;
term(term<0)=0;

Phase=exp(j*2*pi*z*sqrt(term));

Pupil=ones(Nx,Ny);