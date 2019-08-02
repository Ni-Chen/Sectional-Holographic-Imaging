%% -------------------------------------------- Display --------------------------------------------
% Orthoviews(im,[],'Input Image (GT)');
% figure; show3d(gather(im), 0.001); axis normal;
% imdisp(abs(y),'Convolved mag', 1); imdisp(angle(y),'Convolved phase', 1);
% 
% % Back-propagation reconstruction
% im_bp = LinOpAdjoint(H)*y;
% Orthoviews(abs(im_bp),[],'BP Image');

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
        solve_result{imidx} = optSolve;
        
        temp = strrep(solve_name, [obj_name, '_'], '');
        method_name{imidx} = strrep(temp, '.mat', '');
        
        legend_name = [legend_name, method_name{imidx}];
        
%         temp = abs(gather(optSolve.xopt));  %temp = abs(gather(optSolve.xopt));
%         Orthoviews(temp,[], method_name{imidx});
%         temp = (temp-min(temp(:)))/(max(temp(:))-min(temp(:))); 
%         figure('Name', method_name{imidx}); show3d(temp, 0.05); axis normal;
    end
          
    figure('Name', 'Cost evolution'); 
    grid;title('Cost evolution'); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
    for imidx = 1:img_num 
        plot(solve_result{imidx}.OutOp.iternum,solve_result{imidx}.OutOp.evolcost,'LineWidth',1.5);     
        hold all;
    end
    legend(legend_name); 
    set(gcf,'paperpositionmode','auto');
    print('-dpng', ['./output/', obj_name, '_cost.png']);
    
    % Show SNR
%     ls = {'+','o','*','x','v','d','^','s','>','<'};
%     lc = linspecer(4, 'sequential'); 
    ls = {'o-','v-','x-'};
    lc = linspecer(9, 'qualitative');  % qualitative, sequential
    figure('Name', 'SNR'); 
%     set(gca(), 'LineStyleOrder', ls, 'ColorOrder', [0 0 1], 'NextPlot','add');  
    set(gca,'ColorOrder', lc, 'LineStyleOrder', ls, 'NextPlot', 'replacechildren');  
    grid; hold all; title('Evolution SNR'); set(gca,'FontSize', 10);
    for imidx = 1:img_num
        semilogy(solve_result{imidx}.OutOp.iternum,solve_result{imidx}.OutOp.evolsnr,'LineWidth',1.5);
    end
    legend(legend_name,'Location','southeast'); xlabel('Iterations');ylabel('SNR (dB)');
    set(gcf,'paperpositionmode','auto');
    print('-dpng', ['./output/', obj_name, '_snr.png']);
    
%     figure('Name', 'Time cost');    
%     hold on; grid; title('Runing Time');set(gca,'FontSize',12);
%     orderCol = get(gca,'ColorOrder');
%     for imidx = 1:img_num   
%         bar(imidx,[solve_result{imidx}.time],'FaceColor',orderCol(imidx,:),'EdgeColor','k');
%     end
%     set(gca,'xtick',1:img_num);ylabel('Time (s)'); set(gca,'xticklabels', legend_name);
%     set(gca,'XTickLabelRotation',50);
%     set(gcf,'paperpositionmode','auto');
%     print('-dpng', ['./output/', obj_name, '_time.png']);
end
