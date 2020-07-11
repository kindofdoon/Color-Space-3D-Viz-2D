% function color_space_3D_viz_2D

    %%
    
    clear
    clc
    
    %% Constants
    
    JND = 2.3; % Delta_E
    
    %% Inputs
    
    dE = JND * 1; % delta_E color resolution in Lab space; 2.3 = 1.0 JND
    
    spacing = 1; % px, between sections
    border  = 1; % px, around sections
    background = 'transparent'; % 'white', 'black', or 'transparent'
    upsample_target_width = 750; % px, set to NaN to disable
    
    %% Create color gamut
    
    disp('Generating gamut...')
    
%     [RGB, Lab, Gamut] = generate_Lab_gamut(25, 0.9, dE, 0); % medium res
%     [RGB, Lab, Gamut] = generate_Lab_gamut(30, 0.9, dE, 0); % higher res
    [RGB, Lab, Gamut] = generate_Lab_gamut(32, 0.95, dE, 0); % high res
    
    RGB(RGB<0) = 0;
    RGB(RGB>1) = 1;
    
    L = Lab(:,1);
    a = Lab(:,2);
    b = Lab(:,3);
    
    Hue = atan2d(b, a);
    Sat = sqrt(a.^2 + b.^2);
    Val = L;
    
    % Modify gamut as needed
    Hue = deg2rad(Hue + 180);
    LCh = [Val, Sat, Hue];
    
    qty_colors = size(Lab,1);
    qty_values = length(unique(Val));
    
    %% Prepare main image
    
    switch background
        case 'black'
            bac_col = 0;
        case 'white'
            bac_col = 1;
        case 'transparent'
            bac_col = 0;
        otherwise
            error('Unrecognized background')
    end
    
    vals = unique(Val);
    
    % Estimate the size of the image
    sub_image_qty    = length(vals);
    sub_image_width  = length(unique(a));
    sub_image_height = length(unique(b));
    
    margin = 100; % px, to be cropped later
    I = nan([sub_image_height+2*margin, sub_image_width*sub_image_qty+2*margin+spacing*qty_values, 3]);
    
    %% Prepare sub-images
    
    disp('Generating and arranging section sub-images...')
    
    cursor = [1, 1+border]; % row, col
    
    for s = 1 : length(vals)
        
        ind_match = find(Val == vals(s));
        
        L_ = Lab(ind_match, 1);
        a_ = Lab(ind_match, 2);
        b_ = Lab(ind_match, 3);
        
        % Convert color value to pixel location using histogram binning,
        % leveraging uniform spacing
        a_range = [min(a_), max(a_)];
        b_range = [min(b_), max(b_)];
        
        a_bins = linspace(min(a_range)-dE/2, max(a_range)+dE/2, round(abs(diff(a_range))/dE)+2);
        b_bins = linspace(min(b_range)-dE/2, max(b_range)+dE/2, round(abs(diff(b_range))/dE)+2);
        
        dims = [length(b_bins)-1, length(a_bins)-1];
        
        % Initialize sub-image
        I_sub = nan(dims(1), dims(2), 3);
        
        [~, ~, y_ind] = histcounts(b_, b_bins);
        [~, ~, x_ind] = histcounts(a_, a_bins);
        
        y_ind = -y_ind + max(y_ind) + 1; % invert vertical indices
        
        b_mat = nan(size(I_sub,1), size(I_sub,2));
        
        for p = 1 : length(ind_match) % for each pixel
            for cc = 1 : 3
                I_sub(y_ind(p), x_ind(p), cc) = RGB(ind_match(p), cc);
                b_mat(y_ind(p), x_ind(p))     = b_(p);
            end
        end
        
        % Generate vertical offset to align neutral axes of sub-images
        [~, ind] = min(abs(b_mat(:)));
        [offset, ~] = ind2sub(size(b_mat), ind);
        
        px_occupied = sum(sum(~isnan(I_sub(:,:,1))));
        if px_occupied ~= length(ind_match)
            warning(['Mismatch on section #' num2str(s) ': ' num2str(length(ind_match)) ' colors, but ' num2str(px_occupied) ' pixels occupied'])
        end
        
        % Determine horizontal position of sub-image; push as little to the
        % right as possible without overlapping any previous sub-image
        ind_horiz = cursor(2) + 2;
        is_overlapping = 1;
        
        ind_y = round(size(I,1)/2 - offset + [0:dims(1)-1]);
        
        while is_overlapping
        
            Mask_Main = ~isnan(I(:,:,1));
            
            ind_x = ind_horiz+[0:dims(2)-1];
            
            I_copy = nan(size(I));
            I_copy(ind_y, ind_x, :) = I_sub; % patch in the sub-image
            Mask_Sub = ~isnan(I_copy(:,:,1));
            
            Mask_Sum = Mask_Main + Mask_Sub;
            if max(Mask_Sum(:)) == 1
                is_overlapping = 0;
            end
            
            ind_horiz = ind_horiz + 1;
            
%             figure(2)
%                 clf
%                 set(gcf,'color','white')
%                 pcolor(flipud(Mask_Sum))
%                 shading flat
%                 axis tight
%                 axis equal
%                 grid on
%                 colorbar
%                 drawnow
            
        end
        
        ind_horiz = ind_horiz - 1 + spacing;
        cursor(2) = ind_horiz;
        
        % Write sub-image to main image
        ind_x = cursor(2)+[0:dims(2)-1];

        [Y, X] = meshgrid(ind_y, ind_x);
        X = X';
        Y = Y';
        R = I_sub(:,:,1);
        G = I_sub(:,:,2);
        B = I_sub(:,:,3);
        R = R(:);
        G = G(:);
        B = B(:);
        
        for p = 1 : length(ind_match)
            if ~isnan(R(p))
                I(Y(p), X(p), 1) = R(p);
                I(Y(p), X(p), 2) = G(p);
                I(Y(p), X(p), 3) = B(p);
            end
        end
        
%         figure(3)
%             clf
%             set(gcf,'color','white')
%             image(I)
%             axis equal
%             axis tight
%             drawnow
%             pause
        
        % Simple patching - does not support overlapping layers
%         I(ind_y, ind_x, :) = I_sub;
        
    end
    
    %% Crop image
    
    % Adapted from snap.m
    
    CC = I(:,:,1); % pick one color channel arbitrarily
    
    % Scan cols
    scan_col = zeros(1,size(CC,1));
    for a = 1:size(CC,1)
        extract = CC( max([1,a-border]) : min([size(CC,1),a+border]),:);
        scan_col(a) = min(min(isnan(extract)));
    end

    % Scan rows
    scan_row = zeros(1,size(CC,2));
    for a = 1:size(CC,2)
        extract = CC(:, max([1,a-border]) : min([size(CC,2),a+border]));
        scan_row(a) = min(min(isnan(extract)));
    end
    
    I = I(scan_col==0,scan_row==0,:); % crop image
    
    %% Upsample image
    
    if size(I,2) < upsample_target_width
        upsample_scale = floor(upsample_target_width / size(I,2));
        I = imresize(I, upsample_scale, 'nearest');
    else
        upsample_scale = 1;
    end
    
    %% Generate transparency layer
    
    if strcmp(background, 'transparent')
        Alpha = double(~isnan(I(:,:,1)));
    else
        Alpha = ones(size(I(:,:,1)));
    end
    
    %% Fill in background color
    
    for y = 1 : size(I,1)
        for x = 1 : size(I,2)
            if isnan(I(y,x,1))
                I(y,x,:) = bac_col;
            end
        end
    end
    
    %% Show and save final image
    
    figure(1)
        clf
        set(gcf,'color','white')
        set(gcf,'GraphicsSmoothing','off') % no anti-aliasing, preserve hard edges
        image(I)
        axis tight
        axis equal
    
    imwrite(I, ['color_gamut_Lab_' num2str(qty_colors) '_colors_' num2str(dE/JND) '_JND_' num2str(dE) '_dE_' num2str(upsample_scale) 'x_scale_' background '.png'], 'Alpha', Alpha)
    
    
    
    
%     end



















































