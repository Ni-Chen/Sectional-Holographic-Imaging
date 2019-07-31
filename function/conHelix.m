function obj3d = conHelix(Nxy, Nz, ra, turns)
%{
------------------------------------------------
Generates conical helix object

Inputs : 
    ra -> ratio of the radius
    turns -> number of rings

Example : 
    im = conHelix(128,128,0.8,6);

Copyright (C) 2019, Ni Chen, nichen@snu.ac.kr
------------------------------------------------
%}

    t = 0:0.0001:1;
    r = Nxy/2*ra;  % r<Nxy/2 , 0<ra<1

    obj3d = zeros(Nxy,Nxy,Nz);

    Nz_off = round(Nz*0.1);   % off size of z
    c = 1/(Nz-Nz_off*2);

    x = round(r*t.*sin(t*2*pi*turns)) + Nxy/2;
    y = round(r*t.*cos(t*2*pi*turns)) + Nxy/2;
    z = round(t/c) + Nz_off;

    xyz = sub2ind(size(obj3d), x, y, z);

    obj3d(xyz) = 1;
end
