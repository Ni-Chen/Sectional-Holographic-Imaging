function pinvQ=MyPinvQ4DeblurringImageSpace(Q,lambda)

if (nargin < 2), lambda=0; end;

[nx,ny,nzi,nzo]=size(Q);
pinvQ=zeros(nx,ny,nzo,nzi);

if (lambda == 0), func=@(a) pinv(a);  
else func=@(a) inv(a'*a + eye(nzi)*lambda)*a'; 
end


for y=1:ny
    for x=1:nx
        pinvQ(x,y,:,:)=func(squeeze(Q(x,y,:,:)).').'; 
    end;
    sprintf('%d, %d\n',x,y)
end;

% thresholding
% pinvQ(abs(pinvQ)>10^15)=0;