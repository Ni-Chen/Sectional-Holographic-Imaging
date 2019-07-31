function obj3d = randomScatter(Nxy,Nz,sr)
%{
------------------------------------------------
Generates circluar helix object

Inputs: 
    sr -> radius of one single scatter

Example: 
    im = randomScatter(128,128,1);

Copyright (C) 2019, Ni Chen, nichen@snu.ac.kr
------------------------------------------------
%}

%     xy=linspace(-1, 1, Nxy);
%     z=linspace(-1, 1, Nz);
%     dx = 2/Nxy;
%     dz = 2/Nz;
%     
%     [X,Y,Z]=meshgrid(xy,xy,z);
%     
%     xy_idx = randperm(Nxy, round(sqrt(Ns)));
%     z_idx = randperm(Nz, round(sqrt(Ns)));
%     
%     X-X(xy_idx)
%     
%     [~,ph1,r1]=cart2sph(Y,Z,X);
%     [~,ph2,r2]=cart2sph(X,Y,Z);
%     [~,ph3,r3]=cart2sph(X,Z,Y);
    
%     im=(1+cos(4*w*ph1)) + (1+cos(4*w*ph2)) + (1+cos(4*w*ph3));
%     im = rand(r1);
%     obj3d = im./(1+exp(s*(sqrt(X.^2+Y.^2+Z.^2)-x0)));
    

    obj3d = zeros(Nxy, Nxy, Nz);
    for iz = 1:Nz
        ix = ceil(rand(1)*Nxy);
        iy = ceil(rand(1)*Nxy);

        ix_range = ix:ix+sr;
        iy_range = iy:iy+sr;

        ix_range(ix_range>Nxy) = Nxy;
        iy_range(iy_range>Nxy) = Nxy;

        obj3d(iy_range, ix_range, iz) = rand(1);     % rand(1)
    end
    
end

