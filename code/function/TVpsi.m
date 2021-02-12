%{
  Total variation for 3D image denoising
  y = TVpsi(x, lambda, tau, iter, Nvy, Nvx, Nvz)
 
  The output y approximately minimizes the Rudin-Osher-Fatemi(ROF) denoising model
     arg min_y = TV(y) + lambda/2 || x - y ||^2_2,
  where TV(y) is the total variation of y. tau specifies the stopping tolerance (default 1e-2).

  The minimization is solved using Chambolle's method,
  - A. Chambolle, "An Algorithm for Total Variation Minimization and Applications," J. Math. Imaging
    and Vision 20(1-2):89-97, 2004.
%}  

function y = TVpsi(x, lambda, tau, iter, Nvy, Nvx, Nvz)
	X = reshape(x, Nvy, Nvx, Nvz);
	Y = ProjTV(X, tau, lambda*0.5, iter);
	y = reshape(Y, Nvx*Nvy*Nvz, 1);
end

function p = ProjTV(u, tau, lambda, iter)
    [ny, nx, nz] = size(u);
    pn = zeros(ny, nx, nz, 3);
    div_pn = zeros(ny, nx, nz);
    tmp = pn;

%     k=0;
%     cont = 1;
%     while cont
%         k = k+1;
        
    for i = 1:iter
        a = Grad3D(div_pn - u./lambda);
        tmp(:, :, :, 1) = sqrt(a(:, :, :, 1).^2 + a(:, :, :, 2).^2 + a(:, :, :, 3).^2);
        tmp(:, :, :, 2) = tmp(:, :, :, 1);
        tmp(:, :, :, 3) = tmp(:, :, :, 1);
        
        denorm = 1.0 + tau.*tmp;
        pn = (pn + tau.*a)./denorm;
        
%         err = sqrt((-a(:, 1) + tmp(:,1).*pn(:,1)).^2 + (-a(:, 2) + tmp(:,2).*pn(:,2)).^2 + (-a(:, 3) + tmp(:,3).*pn(:,3)).^2);
%         cont = ((k<iter) & (err>tau));
        
        div_pn = Diver3D(pn);
    end

    p = u - lambda.*Diver3D(pn);
end

