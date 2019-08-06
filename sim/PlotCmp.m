%% -------------------------------------------- Display --------------------------------------------
% Orthoviews(im,[],'Input Image (GT)');
% figure; show3d(gather(im), 0.001); axis normal;
% imdisp(abs(y),'Convolved mag', 1); imdisp(angle(y),'Convolved phase', 1);
% 
% % Back-propagation reconstruction
% im_bp = LinOpAdjoint(H)*y;
% Orthoviews(abs(im_bp),[],'BP Image');

%

% Deconvolution reconstruction comparison
% solve_lst = dir(['./output/', obj_name, '_*.mat']);
img_num = length(solve_lst);

if img_num > 0    
%     reset(gpuDevice(1));
    legend_name = {};
    method_name = {};
    
    for imidx = 1:img_num 
        solve_name = solve_lst(imidx).name;
        load(['./output/', solve_name]);
%         dlmread(['./output/', solve_name]);
        solve_result{imidx} = optSolve;
        
%         temp = strrep(solve_name,  '');
        method_name{imidx} = strrep(solve_name, '.mat', '');
%         method_name{imidx} = solve_name;
        
        legend_name = [legend_name, method_name{imidx}];
        
%         temp = abs(gather(optSolve.xopt));  %temp = abs(gather(optSolve.xopt));
%         Orthoviews(temp,[], method_name{imidx});
%         temp = (temp-min(temp(:)))/(max(temp(:))-min(temp(:))); 
%         figure('Name', method_name{imidx}); show3d(temp, 0.05); axis normal;
    end
          
    figure; 
    grid;title('Cost evolution'); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
    for imidx = 1:img_num 
        plot(solve_result{imidx}.OutOp.iternum,solve_result{imidx}.OutOp.evolcost,'LineWidth',1.5);     
        hold all;
    end
    legend(legend_name); 
    set(gcf,'paperpositionmode','auto');
    print('-dpng', ['./output/', 'cost.png']);
%     hold off;
    
    % Show SNR
    figure; 
    grid; hold all; title('Evolution SNR'); set(gca,'FontSize', 10);
    for imidx = 1:img_num
imidx
        lineStyle = randLineStyle();        
        semilogy(solve_result{imidx}.OutOp.iternum,solve_result{imidx}.OutOp.evolsnr, lineStyle{1},'LineWidth',1.5);
    end
    legend(legend_name,'Location','southeast'); xlabel('Iterations');ylabel('SNR (dB)');
    set(gcf,'paperpositionmode','auto');
    print('-dpng', ['./output/', 'snr.png']);


%     figure();   
%     title('Time cost')
%     hold on; grid; title('Runing Time');set(gca,'FontSize',12);
%     orderCol = get(gca,'ColorOrder');
%     for imidx = 1:img_num   
%         bar(imidx,[solve_result{imidx}.time],'FaceColor',orderCol(imidx,:),'EdgeColor','k');
%     end
%     set(gca,'xtick',1:img_num);ylabel('Time (s)'); set(gca,'xticklabels', legend_name);
%     set(gca,'XTickLabelRotation',50);
%     set(gcf,'paperpositionmode','auto');
%     print('-dpng', ['./output/', 'time.png']);
end
