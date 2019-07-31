classdef LinOpWavePropKernel < LinOp
%{
    LinOpFreeSpaceProp: Free space propagation implemented by Rayleigh sommerfeld diffraction with 
    angular spectrum approach.
    
    :param Ny, Nx: lateral size of the 3D volume
    :param ppxy: pixel pitch along the lateral axises
    :param zd: distances between the 3D object to the detector
    :param lambda: wave length
    
    All attributes of parent class :class: `Map` are inherited. 
    
    Copyright (C) 2019
    Ni Chen, nichen@snu.ac.kr
    
    This program is free software: you can redistribute it and/or modify it under the terms of the 
    GNU General Public License as published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
    without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
    the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along with this program. If 
    not, see <http://www.gnu.org/licenses/>.
%}
    
    properties
        Ny;
        Nx;
        ppxy;
        zd;
        lambda = 532e-9;  
        mtf;  % Fourier transform of the psf
        type = 0; % 0:'RS', 1:'Fresnel', 2:'Franherfor'
    end
%% Constructor    
    methods
        function this = LinOpWavePropKernel(lambda, Ny, Nx, ppxy, zd)
            Nz = length(zd);
            this.lambda = lambda;
            this.sizein = [Ny Nx Nz];
            this.sizeout = [Ny Nx];
            this.ppxy = ppxy;
            this.zd = zd;            
                        
            % Sampling at the spectrum coordinates   
            KY = ((1:Ny)-round(Ny/2))*(1/(Ny*ppxy));
            KX = ((1:Nx)-round(Nx/2))*(1/(Nx*ppxy));      
            [kx, ky] = meshgrid(KX, KY);
            
            otf3d = zeros(Ny, Nx, Nz);
            switch this.type
                case 0  % "Rayleigh Sommerfeld" diffraction kernel
                    k = 1/lambda;
                    r = k.^2 - (kx.^2 + ky.^2);
                    r(r<0) = 0;                    
                    
                    for iz = 1:Nz  % 3D
                        otf_rs = exp(1i*2*pi*zd(iz)*sqrt(r)); 
                        otf_rs = ifftshift(otf_rs);
                        otf3d(:,:,iz) = otf_rs;
                    end
                    
                case 1  % Fresnel diffraction
                    k = 1/lambda;                    
                    for iz = 1:Nz  % 3D
                        otf_fr = exp(1i*2*pi*k*zd(iz)).*exp(1i*pi*lambda*zd(iz)*(kx.^2 + ky.^2)); 
                        otf_fr = ifftshift(otf_fr);
                        otf3d(:,:,iz) = otf_fr;
                    end

                case 2  % Fraunhofer
                    k = 2*pi/lambda;                    
                    for iz = 1:Nz  % 3D
                        otf_fra = exp(1i*k*zd(iz)).*exp(1i*k/(2*zd(iz)).*(u.^2+v.^2))/(1i*lambda*zd(iz)); 
                        otf_fra = ifftshift(otf_fra);
                        otf3d(:,:,iz) = otf_fra;
                    end
            end
            
            this.mtf = otf3d;
        end       
    
%         % Complex to real vector
%         function x_vec = C2V(x_complex)
%             x_vec = [real(x_complex); imag(x_complex)];
%         end
%         
%         % real vector to complex
%         function x_complex = V2C(x_vec)
%             [Nyy, Nxx] = size(x_vec);
%             x_complex = zeros(Nyy/2, Nxx);
%             x_complex = x_vec(1:Nyy/2, :) + 1i*x_vec(Nyy/2+1:Nyy, :);
%         end        
    end
    
    methods (Access = protected)
        function out = apply_(this,input)
            % Reimplemented from parent class :class:`LinOp`, object to wavefront
            out = ifft2(fft2(input).*this.mtf);
            out = sum(out, 3);
        end
        function out = applyAdjoint_(this,input)
            % Reimplemented from parent class :class:`LinOp`, wavefront to object
            out = ifft2(fft2(input).*conj(this.mtf));
        end
    end
    
end