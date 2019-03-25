function TV=MyTV3D_conv(x,weights)

if (nargin < 2), weights = 1;  end

[nx,ny,nz]=size(x);
TV=zeros(nx,ny,nz,3);

TV(:,:,:,1)=circshift(x,[-1 0 0])-x;
TV(nx,:,:,1)=0.0;

TV(:,:,:,2)=circshift(x,[0 -1 0])-x;
TV(:,ny,:,2)=0.0;

TV(:,:,:,3)=circshift(x,[0 0 -1])-x;
TV(:,:,nz,3)=0.0;

[nx_w ny_w nz_w n_w]=size(weights);
if([nx_w ny_w nz_w n_w]==[1 3 1 1])
    TV(:,:,:,1)=TV(:,:,:,1).*weights(1);
    TV(:,:,:,2)=TV(:,:,:,2).*weights(2);
    TV(:,:,:,3)=TV(:,:,:,3).*weights(3);
else
    TV=TV.*weights;
end;
