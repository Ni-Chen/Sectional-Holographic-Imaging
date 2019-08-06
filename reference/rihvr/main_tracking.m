% Track results of GPU IIPE Code
% Written by: Kevin Mallery, July 2015
% Modified by Mallery and You, April 2017

close all;
clearvars();

%% Set Parameters
% Import the centroid data
master_pathn = 'Outputs/';
file_base = 'centroids_%04d.csv';

% Hologram data paraemters
dz = 30;
startz = 0;
reso = 1.1;
true_fps = 100;
num_files = 2000;

% Tracking parameters
maxdisp = 15;
param.mem = 5;
param.good = 200;
param.dim = 3;
param.quiet = 0;
maxsteps = 2000;

% Post processing filter parameters
filtersize = 9;
BinWidth = 0.0001;

% data_folders used to enable processing of multiple sequences at once
data_folders = {''};
start_files = [1, 1, 1, 1, 1, 1];

%% Track the data and apply filter
for f = 1:length(data_folders)
    pathn = [master_pathn, data_folders{f}, '/'];
    
    start_file = start_files(f);
    
    xyzs = [];
    num_particles = 0;
    
    %% Read file data
    fprintf('Opening file %04d of %04d', 0, 0);
    for n = start_file:start_file+num_files-1
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b');
        fprintf('%04d of %04d', n-start_file, num_files);
        
        m = csvread([pathn, sprintf(file_base, n)]);
        file_particles = size(m, 1);
        xyzs(num_particles+1:num_particles+file_particles, 1) = m(:,1)*reso;
        xyzs(num_particles+1:num_particles+file_particles, 2) = m(:,2)*reso;
        xyzs(num_particles+1:num_particles+file_particles, 3) = m(:,3)*dz + startz;
        xyzs(num_particles+1:num_particles+file_particles, 4) = n;
        num_particles = num_particles + file_particles;
    end
    fprintf('\n');
    
    % Save results to allow this step to be skipped in the future
    save([master_pathn, data_folders{f}, '/xyzs_raw.mat'], 'xyzs');
    
    %% Tracking
    load([master_pathn, data_folders{f}, '/xyzs_raw.mat'], 'xyzs');
    
    num_del=find(xyzs(:,3)<1200);
    xyzs(num_del,:)=[];
    
    % Run the tracking code
    tracks = track(xyzs,maxdisp,param, maxsteps);
    
    % Save results to allow this step to be skipped in the future
    save([master_pathn, data_folders{f}, '/tracks_nofilt.mat'], 'tracks');
    
    %% Filter the results
    load([master_pathn, data_folders{f}, '/tracks_nofilt.mat'], 'tracks');
    
    lengths = zeros(1, max(tracks(:,5)));
    for i = 1:max(tracks(:,5))
        idx = find(tracks(:,5) == i);
        x = tracks(idx,1);
        y = tracks(idx,2);
        z = tracks(idx,3);
        
        [x,y,z] = avefilter(x,y,z, filtersize);
        
        tracks(idx,1) = x;
        tracks(idx,2) = y;
        tracks(idx,3) = z;
        lengths(i) = length(idx);
    end
    
    %% Plot the outputs
    num_tracks = max(tracks(:,5));
    fprintf('Found %d tracks\n', num_tracks);
    figure;
    title(data_folders{f}, 'Interpreter', 'none');
    for i = 1:num_tracks
        parts = tracks(tracks(:,5)==i, :);
        plot3(parts(:,1), parts(:,2), parts(:,3),'.-');
        grid on;
        hold on
    end
    axis equal;
    xlabel('x (um)');
    ylabel('y (um)');
    zlabel('z (um)');
    set(gca, 'ydir', 'reverse');
    
    multi_tracks{f} = tracks;
    
    fprintf('Saving to %s\n', [master_pathn, data_folders{f}, 'tracks.mat']);
    save([master_pathn, data_folders{f}, '/tracks.mat'], 'tracks');
end

%% Plots of Velocity Components

% Velocity component scatter plot
% Plot is done within KPM_analyze_tracks which also determines velocity
figure;
hold all;
for f = 1:length(data_folders)
    load([master_pathn, data_folders{f}, '/xyzs_raw.mat']);
    fps = true_fps*ones(num_files, 1);
    XYZUVW{f} = KPM_analyze_tracks(multi_tracks{f}, fps);
end
legend(data_folders);

% Another velocity component scatter plot
hold all;
figure;
for f = 1:length(data_folders)
    plot3(XYZUVW{f}(:,4), XYZUVW{f}(:,5), XYZUVW{f}(:,6), '.');
    hold on
end
xlabel('U (um/sec)');
ylabel('V (um/sec)');
zlabel('W (um/sec)');
title('Velocities', 'Interpreter', 'none');

% Determine velocity component statistics
for f = 1:length(data_folders)
    fprintf('for %s:\n', data_folders{f});
    aveU = mean(XYZUVW{f}(:, 4));
    stdU = std(XYZUVW{f}(:, 4));
    U_RMS=sqrt(sum(XYZUVW{f}(:,4).^2)/length(XYZUVW{f}));
    fprintf('  aveU = %f, stdU = %f, U_RMS = %f\n', aveU, stdU, U_RMS);
    aveV = mean(XYZUVW{f}(:, 5));
    stdV = std(XYZUVW{f}(:, 5));
    V_RMS=sqrt(sum(XYZUVW{f}(:,5).^2)/length(XYZUVW{f}));
    fprintf('  aveV = %f, stdV = %f, V_RMS = %f\n', aveV, stdV, V_RMS);
    aveW = mean(XYZUVW{f}(:, 6));
    stdW = std(XYZUVW{f}(:, 6));
    W_RMS=sqrt(sum(XYZUVW{f}(:,6).^2)/length(XYZUVW{f}));
    fprintf('  aveW = %f, stdW =%f, W_RMS = %f\n', aveW, stdW, W_RMS);
end

% Histogram of speeds
figure;
title('Speed Histograms');
hold on;
speeds = cell(size(XYZUVW));
xyspeeds = cell(size(XYZUVW));
for f = 1:length(data_folders)
    speeds{f} = zeros(size(XYZUVW{f}(:,1)));
    for i = 1:length(speeds{f})
        speeds{f}(i) = norm(XYZUVW{f}(i, 4:6));
        xyspeeds{f}(i) = norm(XYZUVW{f}(i,4:5));
    end
    
    h = histogram(speeds{f}, 100, 'Normalization', 'cdf');
    hist_y{f} = h.Values;
    hist_x{f} = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
end
legend(data_folders, 'Interpreter', 'none');
xlabel('Speed (um/sec)');
ylabel('CDF');

% Histogram of speeds as a line plot
figure;
title('Speed Histograms');
hold on;
for f = 1:length(data_folders)
    plot(hist_x{f}, hist_y{f});
end
legend(data_folders, 'Interpreter', 'none');
xlabel('Speed (um/sec)');
ylabel('CDF');

% Speed histogram with PDF normalization
figure;
title('Speed Histograms');
hold on;
speeds = cell(size(XYZUVW));
xyspeeds = cell(size(XYZUVW));
for f = 1:length(data_folders)
    speeds{f} = zeros(size(XYZUVW{f}(:,1)));
    for i = 1:length(speeds{f})
        speeds{f}(i) = norm(XYZUVW{f}(i, 4:6));
        xyspeeds{f}(i) = norm(XYZUVW{f}(i,4:5));
    end
    
    h = histogram(speeds{f}, 100, 'Normalization', 'pdf');
    hist_y{f} = h.Values;
    hist_x{f} = (h.BinEdges(1:end-1)+h.BinEdges(2:end))/2;
end
legend(data_folders, 'Interpreter', 'none');
xlabel('Speed (um/sec)');
ylabel('PDF');

% Speed histogram line plot with PDF normalization
figure;
title('Speed Histograms');
hold on;
for f = 1:length(data_folders)
    plot(hist_x{f}, hist_y{f});
end
legend(data_folders, 'Interpreter', 'none');
xlabel('Speed (um/sec)');
ylabel('PDF');

% Histogram of speeds in XY only (exclude Z component)
figure;
title('XY Speed Histograms');
hold on;
for f = 1:length(data_folders)
    histogram(xyspeeds{f}, 'Normalization', 'pdf');
end
legend(data_folders, 'Interpreter', 'none');
xlabel('Speed (um/sec)');
ylabel('Probability');

% Histograms of velocity in each component direction
figure;
title('X-Velocity Histogram');
hold on;
for f = 1:length(data_folders)
    histogram(XYZUVW{f}(:,4), 'Normalization', 'pdf');
end
xlabel('Velocity (um/sec)');
ylabel('Probability');
legend(data_folders, 'Interpreter', 'none');

figure;
title('Y-Velocity Histogram');
hold on;
for f = 1:length(data_folders)
    histogram(XYZUVW{f}(:,5), 'Normalization', 'pdf');
end
xlabel('Velocity (um/sec)');
ylabel('Probability');
legend(data_folders, 'Interpreter', 'none');

figure;
title('Z-Velocity Histogram');
hold on;
for f = 1:length(data_folders)
    histogram(XYZUVW{f}(:,6), 'Normalization', 'pdf', 'BinWidth', BinWidth);
end
xlabel('Velocity (um/sec)');
ylabel('Probability');
legend(data_folders, 'Interpreter', 'none');

%% Plots to search for specific features
displacements = cell(size(XYZUVW));
travel_dists = cell(size(XYZUVW));

% Compare distance travelled to net displacement 
figure;
hold on;
for f = 1:length(data_folders)
    for pid = 1:max(multi_tracks{f}(:,5))
        idx = multi_tracks{f}(:,5) == pid;
        part = multi_tracks{f}(idx, :);
        displacements{f}(pid) = norm(part(end,1:3) - part(1,1:3));
        
        travel_dists{f}(pid) = 0;
        for i = 2:length(part(:,1))
            travel_dists{f}(pid) = travel_dists{f}(pid) + ...
                norm(part(i,1:3)-part(i-1,1:3));
        end
    end
    
    dist_ratios = travel_dists{f} ./ displacements{f};
    
    histogram(dist_ratios, 100, 'Normalization', 'pdf', 'BinWidth', BinWidth);
    xlabel('Ratio of Distance Travelled to Displacement');
end
legend(data_folders, 'Interpreter', 'none');
ylabel('Probability');

% Look for spatial dependence on speed
figure;
hold on;
title('Speed vs. Wall');
for f = 1:length(data_folders)
    plot(XYZUVW{f}(:,3), speeds{f}, '.');
end
xlabel('Z (um)');
ylabel('Speed (um/sec)');
legend(data_folders, 'Interpreter', 'none');

% Look for preferential distribution of particles in space
figure;
title('X-Position Histogram');
hold on;
for f = 1:length(data_folders)
    histogram(XYZUVW{f}(:,1), 'Normalization', 'pdf', 'BinWidth', 10);
end
xlim([0, 2500]);
xlabel('X (um)');
ylabel('Probability');
legend(data_folders, 'Interpreter', 'none');

figure;
title('Y-Position Histogram');
hold on;
for f = 1:length(data_folders)
    histogram(XYZUVW{f}(:,2), 'Normalization', 'pdf', 'BinWidth', 10);
end
xlim([0, 2500]);
xlabel('Y (um)');
ylabel('Probability');
legend(data_folders, 'Interpreter', 'none');

figure;
title('Z-Position Histogram');
hold on;
for f = 1:length(data_folders)
    histogram(XYZUVW{f}(:,3), 'Normalization', 'pdf', 'BinWidth', 10);
end
xlabel('Z (um)');
ylabel('Probability');
legend(data_folders, 'Interpreter', 'none');