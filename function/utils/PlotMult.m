%% -------------------------------------------- Display --------------------------------------------
% Orthoviews(im,[], 'Input Image (GT)');
% figure; show3d(gather(im),  0.001); axis normal;
% imdisp(abs(y), 'Convolved mag',  1); imdisp(angle(y), 'Convolved phase',  1);
% 
% % Back-propagation reconstruction
% im_bp = LinOpAdjoint(H)*y;
% Orthoviews(abs(im_bp), [], 'BP Image');

img_num = length(solve_lst);
if img_num > 0    
    legend_name = {};
    method_name = {};
    
    for imidx = 1:img_num 
        solve_name = solve_lst(imidx).name;
        load([data_dir, solve_name]);
        solve_result{imidx} = optSolve;
        
        method_name{imidx} = strrep(solve_name, '.mat',  '');
        
        legend_name = [legend_name, method_name{imidx}];        
    end
        
    %% Show objective cost value
    figure('Name', 'Cost');
%     figure('units','normalized','outerposition',[0 0 1 1]); set(gca,'FontSize', 12);
    grid;title('Cost evolution'); xlabel('Iterations');ylabel('Cost');
    for imidx = 1:img_num 
        lineStyle = randLineStyle();    
        plot(solve_result{imidx}.OutOp.iternum,solve_result{imidx}.OutOp.evolcost, lineStyle{1}, 'LineWidth',1.5);     
        hold all;
    end
    legend(legend_name); 
    set(gcf,'paperpositionmode', 'auto');
    print('-dpng',  [data_dir, out_name,'_cost.png']);
    
    %% Show SNR
    if isSim     
        figure('Name', 'SNR');
%         figure('units','normalized','outerposition',[0 0 1 1]); set(gca,'FontSize', 12);
        grid; hold all; title('Evolution SNR'); xlabel('Iterations');ylabel('SNR (dB)');
        for imidx = 1:img_num
            lineStyle = randLineStyle();
            semilogy(solve_result{imidx}.OutOp.iternum,solve_result{imidx}.OutOp.evolsnr, lineStyle{1}, 'LineWidth', 1.5);
        end
        legend(legend_name,'Location', 'southeast');
        set(gcf, 'paperpositionmode', 'auto');
        print('-dpng', [data_dir, out_name, '_snr.png']);
    end

    
    %% show time cost
    figure('Name', 'Time');
%     figure('units','normalized','outerposition',[0 0 1 1]); set(gca,'FontSize', 12);
    hold on; grid; title('Runing Time'); set(gca,'xtick',1:img_num); ylabel('Time (s)');
    orderCol = linspecer(img_num, 'qualitative');
    for imidx = 1:img_num   
        bar(imidx,[solve_result{imidx}.time], 'FaceColor',orderCol(imidx,:), 'EdgeColor', 'k');
    end
    set(gca,'xticklabels',  legend_name);
    set(gca,'XTickLabelRotation', 50);
    set(gcf,'paperpositionmode', 'auto');
    print('-dpng', [data_dir, out_name, '_time.png']);
end
