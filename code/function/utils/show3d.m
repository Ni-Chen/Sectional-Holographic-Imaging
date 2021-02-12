function show3d(obj_vol, alpham)
%     obj_vol = permute(obj_vol, [2 3 1]);  % Complex images
    [Ny, Nx, Nz] = size(obj_vol);
    
    axes('fontsize', 14, 'xtick', 0:floor(Nx/4):Nx, 'ytick', 0:floor(Ny/4):Ny, 'ztick', 0:floor(Nz/4):Nz);
    vol3d('CData', obj_vol, 'texture', '3D');
%     view(-195, 18);
    view(-25, 20);
%     view(3);   % -37.5, 30
    
    colormap('hot');
    colormap(1 - colormap);
    colorbar;
    
    xlabel('x', 'fontsize', 14);
    ylabel('y', 'fontsize', 14);
    zlabel('z', 'fontsize', 14);
    
    axis([0 Nx 0 Ny 0 Nz]);
    axis tight;
    axis equal;
   
    box on,
    ax = gca;
    ax.BoxStyle = 'full';
    
    grid on; 
    grid minor;
    % zoom(0.7)
    alphamap('decrease', alpham); 
end