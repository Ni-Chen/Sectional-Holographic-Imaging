function p=MyProjectionTV(g,tau,lam,iter,weights)

if (nargin < 5), weights = 1;  end

[nx,ny,nz]=size(g);
pn=zeros(nx,ny,nz,3);
div_pn=zeros(nx,ny,nz);
b=pn;


for i=1:iter

    a=MyTV3D_conv(div_pn-g./lam,weights);

    b(:,:,:,1)=sqrt(a(:,:,:,1).^2+a(:,:,:,2).^2+a(:,:,:,3).^2);
    b(:,:,:,2)=b(:,:,:,1);
    b(:,:,:,3)=b(:,:,:,1);
    pn=(pn+tau.*a)./(1.0+tau.*b);

%     b(:,:,:,1)=a(:,:,:,1).^2+a(:,:,:,2).^2+a(:,:,:,3).^2;
%     b(:,:,:,2)=b(:,:,:,1);
%     b(:,:,:,3)=b(:,:,:,1);    
%     pn=(pn+tau.*a)./sqrt(1.0+(tau.^2).*b);


    div_pn=MyDiv3D(pn);
end;

p=lam.*MyDiv3D(pn);

