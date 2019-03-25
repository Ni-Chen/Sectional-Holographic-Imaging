function Q=MyQ4DeblurringImageSpaceMultiAngles(func,nx,ny,nz,deltaX,deltaY,deltaZ,angles,mag)

if (nargin < 9), mag = 1;  end

Q=MyQ4DeblurringImageSpace(func,nx,ny,nz);

na=size(angles,1);

% deltas spatial domain
delta=zeros(nx,ny);
delta(1,1)=1/na;
psfs=zeros(nx,ny,nz);

for z=1:nz
    shifts=deltaZ*mag*(z-1)*tan(-angles);
    shifts(:,1)=shifts(:,1)./(deltaX*mag);
    shifts(:,2)=shifts(:,2)./(deltaY*mag);
    shifts=round(shifts);
    for s=1:na
        psfs(:,:,z)=psfs(:,:,z)+circshift(delta,shifts(s,:));
    end;
end;

% figure;imagesc(plotdatacube(psfs,nz));axis image;

% frequency domain
for z=1:nz
    deltasQ(:,:,z)=fft2(psfs(:,:,z));
end;

for i=1:nz
    for o=1:nz
        Q(:,:,i,o)=Q(:,:,i,o).*deltasQ(:,:,o);
    end;
    deltasQ=circshift(deltasQ,[0 0 1]);
end;