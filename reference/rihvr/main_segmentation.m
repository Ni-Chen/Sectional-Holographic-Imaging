close all;
clearvars();

%% Parameters

in_pathn = 'Outputs/';
in_fn = [in_pathn, 'inverse_voxels_%04d.csv'];
out_fn = [in_pathn, 'centroids_%04d.csv'];

frame_list = 1:1;

use_valid_area = false;
use_valid_shape = false;
use_valid_intensity = true;
plot_figures = false;
min_area = 100;
min_intensity = 1;

%% Do not modify below here

xyzs = [];
for fidx = 1:length(frame_list)
    fprintf('Frame %d of %d\n', fidx, length(frame_list));
    frame_id = frame_list(fidx);
    
    data = csvread(sprintf(in_fn, frame_id), 1, 0);
    data(:,1:3) = data(:,1:3)+1; % Switch from C to matlab indexing
    
    if use_valid_intensity
        idx = data(:,end) > min_intensity*max(data(:,end))/256;
        data = data(idx, :);
    end
    
    % Determine volume size from data rather than as input
    vs = max(data, [], 1);
    volume_size = [vs(2), vs(1), vs(3)];
    if (max(data(:,1) > volume_size(2))) || (max(data(:,2) > volume_size(1))) || (max(data(:,3) > volume_size(3)))
        error('Data is out of bounds, check volume_size parameter');
    end
    
    num_points = size(data,1);
    labels = (1:num_points)';
    
    % Giant volume
    volume = false(volume_size);
    
    fprintf('Point %08d of %08d',0,0);
    for i =  1:num_points
        if ~mod(i, 100)
            fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
            fprintf('%08d of %08d', i, num_points);
        end
        volume(data(i,2), data(i,1), data(i,3)) = true;
    end
    fprintf('\n');
    
    props = regionprops(volume, 'Centroid', 'BoundingBox', 'PixelIdxList', 'Area');
    num_objects = length(props);
    clear volume;
    
    %save(sprintf([in_pathn, 'regionprops_%04d.mat'], frame_id), 'props');
    
    cents = cat(1, props.Centroid);
    s = repmat(frame_id, size(cents, 1), 1);
    xyzs = [xyzs; [cents, s]];
    %save([in_pathn, 'xyzs.mat'], 'xyzs');
    
    if use_valid_area || use_valid_shape
        if use_valid_area
            valid_areas = find(cat(1, props.Area) > min_area);
        else
            valid_areas = 1:num_objects;
        end
    
        if use_valid_shape
            valid_shape = [];
            fprintf('prop %04d of %04d', 0, 0);
            for i = 1:length(props)
                fprintf('\b\b\b\b\b\b\b\b\b\b\b\b');
                fprintf('%04d of %04d', i, length(props));
                pix = props(i).PixelIdxList;
                [y,x,z] = ind2sub(volume_size, pix);
                
                bb = ceil(props(i).BoundingBox);
                img2d = false([bb(5), bb(4)]);
                xx = x - bb(1) + 1;
                yy = y - bb(2) + 1;
                idx = unique(sub2ind([bb(5), bb(4)], yy, xx));
                img2d(idx) = true;
                p2 = regionprops(img2d, 'Eccentricity');
                
                if p2(1).Eccentricity > 0.8
                    valid_shape = [valid_shape, i];
                end
            end
            fprintf('\n');
        else
            valid_shape = 1:num_objects;
        end
    
        valid_idx = intersect(valid_areas, valid_shape);
    else
        valid_idx = 1:num_objects;
    end
    
    if plot_figures
        particle_list = struct([]);
        figure;
        for i = 1:length(valid_idx)
            pix = props(valid_idx(i)).PixelIdxList;
            [y,x,z] =ind2sub(volume_size, pix);
            plot3(x, y, z, '.');
            hold all;

            [coef,~,latent,~,explained,mu] = pca([x,y,z]);
            dir = coef(:,1);
            center = mu;

            particle_list(i).centroid = props(valid_idx(i)).Centroid;
            particle_list(i).dir = dir;
        end
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
    end
    
    if plot_figures
        figure;
        for i = 1:length(props)
            pix = props(i).PixelIdxList;
            [y,x,z] = ind2sub(volume_size, pix);
            plot3(x, y, z, '.');
            hold all;
        end
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
    end
    
    props = props(valid_idx);
    cents = cat(1, props.Centroid);
    csvwrite(sprintf(out_fn, frame_id), cents);
end
