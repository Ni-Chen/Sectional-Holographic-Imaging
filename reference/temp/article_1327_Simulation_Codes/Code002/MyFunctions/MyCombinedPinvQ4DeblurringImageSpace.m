function [pinvQ combQ]=MyCombinedPinvQ4DeblurringImageSpace(Q,lambda)

if (nargin < 2), lambda=0; end;

[nx,ny,nzi,nzo]=size(Q);

pinvQ=MyPinvQ4DeblurringImageSpace(Q,lambda);


combQ=zeros(nx,ny,nzi,nzi);


for z1=1:nzi
    for z2=1:nzi
        combQ(:,:,z2,z1)=combQ(:,:,z2,z1)+...
            sum(squeeze(pinvQ(:,:,:,z1)).*squeeze(Q(:,:,z2,:)),3);    
    end;
end;