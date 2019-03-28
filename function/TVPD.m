function rec = TVPD(k_hat, dat, nfx, f_true, nhz, max_pd_iter, max_cg_iter, num_iter)
%TVPD   Total variation regularization with primal-dual Newton's method.
%   rec = TVPD(K, d, nfx, nhz) returns a N-by-nfx matrix, which minimizes
%   the function: 
%       T(u) = ||K * u - d||^2/2 + alpha * J(u),
%   where K is a discretized integral operator, d is a discrete data, 
%   ||.|| denotes the l^2 norm, alpha is a positive regularization 
%   parameter, J is a smooth approximation to the Total Variation
%   function, and nhz is the number of sections to be reconstructed.
%
%   The primal-dual system is solved by Newton's method.
%
%   CG method is used to solve the set of linear equations existing in the 
%   Newton iteration. In each iteration, the primal-dual Newton's method 
%   gives a preliminary guess of the solution, and the GPLS function updates    
%   the solution to make it nonnegative.
%
%   Class support for inputs K, d, nfx, nhz:
%   	cell; float: double, single
%
%   This function is developed by Tristan X. Zhang.
%   See also CG.
  
  cg_steptol = 1; %1e-7;
  cg_residtol = 1; %1e-7
  cg_out_flag = 0;  %if flag = 0, output CG convergence info.
  reset_flag = 1; %input(' Enter 1 to reset; else enter 0: ');
  if exist('f_alpha','var')
    e_pd = [];
  end
  
  if reset_flag == 1
    alpha = 500; %regularization parameter alpha
    beta = 1e-50; %TV smoothing parameter beta
  
    % Discretize first derivative operators
    n = nfx;
    Delta_x = 1 / n; %dx
    Delta_y = Delta_x/2; %dy
    Delta_xy = Delta_x * Delta_y; %dxdy

    Dx1 = spdiags([-ones(n-1,1) ones(n-1,1)], [0 1], n-1,n) / Delta_x;
    I_truncx1 = spdiags(ones(nhz*n-1,1), 0, nhz*n-1,nhz*n);
    Dx = kron(Dx1, I_truncx1); %discretization of the gradient operator x

    Dy1 = spdiags([-ones(nhz*n-1,1) ones(nhz*n-1,1)], [0 1], nhz*n-1,nhz*n) / Delta_y;
    I_truncy1 = spdiags(ones(n-1,1), 0, n-1,n);
    Dy = kron(I_truncy1, Dy1); %discretization of the gradient operator y
  
    % Initialization
    alph_Dxy = alpha * Delta_xy;

    %% Set up parameters
    regvalue = 0;
    regmtx = regvalue*[0 -1 0;-1 4 -1;0 -1 0];
    k_hat_sq = genemtx(k_hat, [2;nhz], regmtx); %k_hat'*k_hat+alpha*L
    Kstar_d = mnttimesv(k_hat, dat, [2;nhz]); %compute k_hat'*dat;
    Kstar_d = reshape(Kstar_d,nfx,[])';
    Kstar_d = Kstar_d(:);
    % Solve nonnegatively constrained least squares minimization problem
    f_pd = gpls(k_hat, dat, nfx*nfx*nhz, nhz, num_iter);
    f_pd = reshape(f_pd, nfx, [])'; %primal variable
    u = zeros(n*nhz-1,n-1); %dual variables u and v
    v = zeros(n*nhz-1,n-1);

  end
 
  pd_gradnorm = []; %gradient norm
  snorm_vec = []; %step norm vector
  
%% Primal-dual Newton iteration
  for pd_iter = 1:max_pd_iter
    Dxf = Dx * f_pd(:); %df/dx
    Dyf = Dy * f_pd(:); %df/dy
    Df_squared = Dxf.^2 + Dyf.^2;
    psi_1 = 1 ./ sqrt(Df_squared + beta^2);
    psi_2 = -.5 * (Df_squared + beta^2).^(-1.5);

    Binv = spdiags(psi_1, 0, (n-1)*(nhz*n-1),(n-1)*(nhz*n-1));
    
    % Construct matrix Lbar 
    uu = u(:); vv = v(:);
    w = 2 * psi_2 ./ psi_1.^2;
    E11 = spdiags(1 + w.*uu.*Dxf, 0, (n-1)*(nhz*n-1),(n-1)*(nhz*n-1));
    E12 = spdiags(w.*uu.*Dyf, 0, (n-1)*(nhz*n-1),(n-1)*(nhz*n-1));
    E21 = spdiags(w.*vv.*Dxf, 0, (n-1)*(nhz*n-1),(n-1)*(nhz*n-1));
    E22 = spdiags(1+ w.*vv.*Dyf, 0, (n-1)*(nhz*n-1),(n-1)*(nhz*n-1));
    Lbar = Dx'*Binv*E11*Dx + Dy'*Binv*E22*Dy ...
 	+ Dx'*Binv*E12*Dy + Dy'*Binv*E21*Dx;
    Lbar = (Lbar + Lbar')/2;

    % Solve primal-dual system
    fvec = f_pd(:); %primal variable in vector form
    fvec1 = reshape(fvec, n*nhz,[])';
    fvec1 = fvec1(:);
    KstarKf = mntimesv(k_hat_sq, fvec1, [nhz;nhz]); %k_hat_sq'*fvec1;
    KstarKf = reshape(KstarKf,nfx,[])';
    KstarKf = KstarKf(:);
    r = Kstar_d - KstarKf - alph_Dxy*(Dx'*Binv*Dxf + Dy'*Binv*Dyf);
    
    gradnorm = norm(r);
    pd_gradnorm = [pd_gradnorm; gradnorm]; %norm of pd gradient 
    
%% Use PCG iteration to solve linear system (K'*K + alpha*Lbar)*Delta_f = r
    fprintf(' ... solving linear system using cg iteration ... \n');
    [Delf,residnormvec,stepnormvec,cgiter] = ...
      cg(k_hat_sq,Lbar,alph_Dxy,r,n,nhz,max_cg_iter,cg_steptol,cg_residtol);

    Delta_f = reshape(Delf,nfx*nhz,[]); %result of the above linear system

%% Compute dual Newton steps. 
    DxDf = Dx*Delf;
    DyDf = Dy*Delf;
    Delta_u = -u(:) + Binv*(Dx*f_pd(:) + E11*DxDf + E12*DyDf); %du
    Delta_u = reshape(Delta_u,n*nhz-1,n-1);
    Delta_v = -v(:) + Binv*(Dy*f_pd(:) + E21*DxDf + E22*DyDf); %dv
    Delta_v = reshape(Delta_v,n*nhz-1,n-1);
    f_pd = f_pd + Delta_f; %update primal variable.
    
    if exist('f_alpha','var')
      e_pd = [e_pd; norm(f_pd - f_alpha,'fro')/norm(f_alpha,'fro')];
    end
    
    %  Perform line search in the dual variables u,v.  Calculate minimum 
    %  rho in (0,1] for which |V(i,j) + rho*Delta_V(i,j)| = 1 for all (i,j). 
    %  Here |.| denotes Euclidean norm on R^2.      
    VdotdV = u.*Delta_u + v.*Delta_v;
    absVsq = u.^2 + v.^2;
    absdVsq = Delta_u.^2 + Delta_v.^2;
    rho = (sqrt(VdotdV.^2 + absdVsq.*(1-absVsq)) - VdotdV) ./ (absdVsq + eps);
    rho_min = min(min(rho(:)),1);
    if rho_min < 1
      rho_min = .9*rho_min; %ensure max_i sqrt(Vx(i)^2 + Vy(i)^2)<1 is maintained
    end
    u = u + rho_min*Delta_u;
    v = v + rho_min*Delta_v;
      
    % Display results 
    snorm = sqrt(norm(Delta_f(:))^2 + norm(Delta_u(:))^2 + norm(Delta_v(:))^2);
    snorm_vec = [snorm_vec; snorm]; %setp norm vector
    
    % Output primal-dual Newton convergence information. 
    fprintf(' PD iter=%3.0f, ||grad||=%6.4e, ||step||=%6.4e, nCG=%3.0f\n', ...
       pd_iter, gradnorm, snorm, cgiter);
       
     % Leave out corner entry f_pd(n,n), since it is coupled to other entries 
     % only through K'*K, and not through Lbar.  
     U_PD = f_pd(1:nhz*n-1,1:n-1); 
     figure(3)
     imagesc(U_PD), colorbar
     suptitle('Reconstruction')
     drawnow
     pause(3)
     
    f_pd = f_pd';f_pd = f_pd(:);
    f_pd = gpls1d(k_hat, dat, nfx*nfx*nhz, nhz, f_pd, num_iter);
    f_pd = reshape(f_pd, nfx, [])'; %update f_pd

  end %for pd iteration

  rel_soln_error = norm(f_pd(:)-f_true(:))/norm(f_true(:))
  rec = f_pd;
  clear max_pd_iter;

end





