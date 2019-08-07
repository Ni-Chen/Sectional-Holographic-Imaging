%% -------------------------------------------- Display --------------------------------------------
% Orthoviews(im,[], 'Input Image (GT)');
% figure; show3d(gather(im),  0.001); axis normal;
% imdisp(abs(y), 'Convolved mag',  1); imdisp(angle(y), 'Convolved phase',  1);
% 
% % Back-propagation reconstruction
% im_bp = LinOpAdjoint(H)*y;
% Orthoviews(abs(im_bp), [], 'BP Image');

%

% Deconvolution reconstruction comparison
% solve_lst = dir(['./output/',  obj_name, '_*.mat']);
img_num = length(solve_lst);

if img_num > 0    
%     reset(gpuDevice(1));
    legend_name = {};
    method_name = {};
    
    for imidx = 1:img_num 
        solve_name = solve_lst(imidx).name;
        load([data_dir, solve_name]);
        solve_result{imidx} = optSolve;
        
        method_name{imidx} = strrep(solve_name, '.mat',  '');
        
        legend_name = [legend_name, method_name{imidx}];
        
%         temp = abs(gather(optSolve.xopt));  %temp = abs(gather(optSolve.xopt));
%         Orthoviews(temp,[],  method_name{imidx});
%         temp = (temp-min(temp(:)))/(max(temp(:))-min(temp(:))); 
%         figure('Name',  method_name{imidx}); show3d(temp, 0.05); axis normal;
    end
          
    figure; 
    grid;title('Cost evolution'); set(gca,'FontSize',12);xlabel('Iterations');ylabel('Cost');
    for imidx = 1:img_num 
        plot(solve_result{imidx}.OutOp.iternum,solve_result{imidx}.OutOp.evolcost,'LineWidth',1.5);     
        hold all;
    end
    legend(legend_name); 
    set(gcf,'paperpositionmode', 'auto');
    print('-dpng',  [data_dir, out_name,'_cost.png']);
%     hold off;
    
    % Show SNR
    figure; 
    grid; hold all; title('Evolution SNR'); set(gca,'FontSize',  10);
    for imidx = 1:img_num
        lineStyle = randLineStyle();        
        semilogy(solve_result{imidx}.OutOp.iternum,solve_result{imidx}.OutOp.evolsnr, lineStyle{1}, 'LineWidth', 1.5);
    end
    legend(legend_name,'Location', 'southeast'); xlabel('Iterations');ylabel('SNR (dB)');
    set(gcf, 'paperpositionmode', 'auto');
    print('-dpng', [data_dir, out_name, '_snr.png']);


    figure;   
    hold on; grid; title('Runing Time'); set(gca,'FontSize',10);
    orderCol = linspecer(img_num, 'qualitative');
    for imidx = 1:img_num   
        bar(imidx,[solve_result{imidx}.time], 'FaceColor',orderCol(imidx,:), 'EdgeColor', 'k');
    end
    set(gca,'xtick',1:img_num);ylabel('Time (s)'); set(gca,'xticklabels',  legend_name);
    set(gca,'XTickLabelRotation',50);
    set(gcf,'paperpositionmode', 'auto');
    print('-dpng', [data_dir, out_name, '_time.png']);
end
