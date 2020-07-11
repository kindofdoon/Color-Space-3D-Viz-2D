% function show_photo_on_color_map

    % Provide a photo and color map. This program provides a "heatmap"
    % showing which colors the photo uses, and how often.

    clear
    clc
    
    %%
    
    addpath('C:\Users\Daniel\Desktop\PAINTERS_COPILOT\') % for access to color_difference.m
    
    %% Inputs
    
    Pho.res_max = 150;
    Map.res_max = Inf;
    
    opacity_if_unused = 0.35;
    
    Lab_precision = 10;
    
    upsample_target_width = 2000; % px, set to NaN to disable

    %% Load photos
    
    [Pho.filename, Pho.pathname, ~] = uigetfile([pwd '\*'],'Select Photo');
%     [Map.filename, Map.pathname, ~] = uigetfile([pwd '\*'],'Select Map');

%     Pho.pathname = 'C:\Users\Daniel\Desktop\PAINTERS_COPILOT\photos\';
%     Pho.filename = '_IMG_2649.jpg';
    Map.pathname = 'C:\Users\Daniel\Desktop\PAINTERS_COPILOT\tests\';
    Map.filename = 'color_gamut_Lab_8436_colors_4.6_dE_1_px_spacing_transparent_1x_upsample.png';
    
    Pho.I_RGB = imread([Pho.pathname, Pho.filename]);
    Map.I_RGB = imread([Map.pathname, Map.filename]);
    
    %% Process photos
    
    Map.I_RGB(Map.I_RGB==0) = 255; % set pure black to pure white for display purposes
    
    % Resize photos
    if max(size(Pho.I_RGB)) > Pho.res_max
        Pho.I_RGB = imresize(Pho.I_RGB, Pho.res_max / max(size(Pho.I_RGB)));
    end
    if max(size(Map.I_RGB)) > Map.res_max
        Map.I_RGB = imresize(Map.I_RGB, Map.res_max / max(size(Map.I_RGB)));
    end
    
    Pho.RGB = double(reshape(Pho.I_RGB, [numel(Pho.I_RGB)/3, 3])) / 255;
    Map.RGB = double(reshape(Map.I_RGB, [numel(Map.I_RGB)/3, 3])) / 255;
    
    Pho.I_Lab = rgb2lab(Pho.I_RGB, 'WhitePoint', 'D65');
    Map.I_Lab = rgb2lab(Map.I_RGB, 'WhitePoint', 'D65');
    
    Pho.L = Pho.I_Lab(:,:,1);
    Pho.a = Pho.I_Lab(:,:,2);
    Pho.b = Pho.I_Lab(:,:,3);
    Pho.Lab = [Pho.L(:), Pho.a(:), Pho.b(:)];
    
    Map.L = Map.I_Lab(:,:,1);
    Map.a = Map.I_Lab(:,:,2);
    Map.b = Map.I_Lab(:,:,3);
    Map.Lab = [Map.L(:), Map.a(:), Map.b(:)];
    
    Pho.hue = atan2d(Pho.Lab(:,3), Pho.Lab(:,2));
    Pho.sat = sqrt(Pho.Lab(:,2).^2 + Pho.Lab(:,3).^2);
    Pho.val = Pho.Lab(:,1);
    
    Map.hue = atan2d(Map.Lab(:,3), Map.Lab(:,2));
    Map.sat = sqrt(Map.Lab(:,2).^2 + Map.Lab(:,3).^2);
    Map.val = Map.Lab(:,1);
    
    %% Discretize data with finite precision for speed
    
    Map.Lab = round(Map.Lab/Lab_precision) * Lab_precision;
    Pho.Lab = round(Pho.Lab/Lab_precision) * Lab_precision;
    
    %%
    
    [Map.Lab_unique, Map.Lab_unique_ind, Map.Lab_unique_ID] = unique(Map.Lab, 'rows');
    
    Map.Lab_unique(1,:) = [Lab_precision, 0, 0];
    
    Map.qty_unique = size(Map.Lab_unique, 1);
    Map.RGB_unique = Map.RGB(Map.Lab_unique_ind,:);
    
    figure(3)
        clf
        hold on
        set(gcf,'color','white')
        scatter3(Map.Lab_unique(:,1), Map.Lab_unique(:,2), Map.Lab_unique(:,3), 10, Map.RGB_unique, 'filled')
        scatter3(Pho.Lab(:,1), Pho.Lab(:,2), Pho.Lab(:,3), 1, 'k')
        axis tight
        axis equal
        axis vis3d
        xlabel('L')
        ylabel('a')
        zlabel('b')
        grid on
        grid minor
    
    %% Match each color in the image to the map
    
    Usage = zeros(size(Map.I_RGB,1), size(Map.I_RGB,2)); % heatmap
    Pho.ID = zeros(size(Pho.Lab,1),1);
    
    for p = 1 : size(Pho.I_Lab,1)*size(Pho.I_Lab,2) % for each photo pixel
        
        % Match color against map
        [Lia, Pho.ID(p)] = ismember(Pho.Lab(p,:), Map.Lab_unique, 'rows');
        
        % Find where this color occurs in map
        ind_lin = find(Map.Lab_unique_ID == Pho.ID(p));
        
        % Increment the usage matrix
        Usage(ind_lin) = Usage(ind_lin) + 1;
        
    end
    
    match_actual = sum(Pho.ID~=0);
    match_goal   = size(Pho.Lab,1);
    
    disp([num2str(match_actual) ' of ' num2str(match_goal) ' (' num2str(round(match_actual/match_goal*1000)/10) '%) photo pixels were matched to the map'])
    
    %% Prepare result image(s)
    
    Usage(Usage>0) = 1;
    I_Usage = repmat(Usage, 1, 1, 3);
    I_Usage = I_Usage ./ max(I_Usage(:)); % normalize
    
    I_Composite = Map.I_RGB;
    is_unused = find(Usage==0);
    is_filled = find(sum(Map.I_RGB,3)~=255*3);
    is_inactive = intersect(is_unused, is_filled);
    
    for cc = 1 : 3
        
        CC = I_Composite(:,:,cc);
        CC(is_inactive) = 128;%CC(is_inactive) * opacity_if_unused;
        I_Composite(:,:,cc) = CC;
        
    end
    
    %% Show results
    
    figure(1)
    clf
    set(gcf,'color','white')
    image(Pho.I_RGB)
    axis tight
    axis equal
    
    figure(2)
    clf
    set(gcf,'color','white')
    subplot(3,1,1)
        image(Map.I_RGB)
        axis tight
        axis equal
    subplot(3,1,2)
        image(I_Usage)
        axis tight
        axis equal
    subplot(3,1,3)
        image(I_Composite)
        axis tight
        axis equal
        
    %% Export result(s)
    
    if size(I_Composite,2) < upsample_target_width
        upsample_scale = floor(upsample_target_width / size(I_Composite,2));
        I_Composite = imresize(I_Composite, upsample_scale, 'nearest');
    else
        upsample_scale = 1;
    end
    
    imwrite(I_Composite, [ regexprep(Pho.filename,'\..+','') '_usage_' num2str(Lab_precision) '_dE_' num2str(upsample_scale) 'x_scale' '.png'])
    
    
    
% end



















































