function [mrk] = f_EOS_mrkID_T(IMG,fig,r_cut)
% Feb. 2018 - X.Gasparutto - HUG
% Identify mrk with radio-opaque balls
% apply threshold in a moving windows

% 1 - Filter
    F.img = wiener2(IMG,[10 10]); 
    
    % 2 - Threshold with moving windows
    % everything below .8*max is put to zeros
    % 0.8 is arbitrary, works fine at the moment (16.02.18)
        % 2.1. Initialisation
        Fmax =max(max(F.img));
%         T.img = F.img;

        % 2.2. Max in moving windows
        % 2.2.1 Initialisation of the window
            % Image Dimension
            ydim = size(F.img,1); 
            xdim = size(F.img,2);
            % windows size in function of image dimension
            dy = floor(ydim/5);
            dx = floor(xdim/2); % xdim/4 does not work better
        
        % 2.2.2. Apply Threshold
        FT.img = F.img;
        for i = 1:2
            for j = 1:5
                yj = (j-1)*dy+1:j*dy;
                xi = (i-1)*dx+1:i*dx;
                % a. Select subpart of image
                tmp_imgF = F.img(yj, xi);
                % b. Local maxima
                tmp_maxF = max(max(tmp_imgF));
                % c. Threshold
                tmp_imgF(tmp_imgF < tmp_maxF * r_cut) = uint16(0);
                % d. Fill T again
                FT.img(yj,xi) = tmp_imgF;
                clear tmp_*
            end
        end
  
    % 3 - Clean Image After Threshold
        fFT  = imfill(double(FT.img),'hole');
        afFT = bwareaopen(fFT,10);
        stats2 = regionprops(afFT,'Area','Centroid','Perimeter','MajorAxisLength','MinorAxisLength','BoundingBox');

        % TMP: loop to get param from stat in tables
        nstat = size(stats2,1);
        for i =1:nstat % find clean way to do that 
            Area(i) = stats2(i).Area;
            local_centre(i,:) = stats2(i).Centroid;
            MajorAxisLength(i) = stats2(i).MajorAxisLength;
            MinorAxisLength(i) = stats2(i).MinorAxisLength;
            Perimeter(i) = stats2(i).Perimeter;
            BoundingBox(i,:) = stats2(i).BoundingBox;
        end
    
    % 4 - Identify fishing weights
        % 4.1. constraint on minor and major axis length
        %should be equal if circle, 20% difference is acceptable
        tmp1 = MinorAxisLength./MajorAxisLength;
        tmp2 = find(tmp1>0.8);
        tmp3 = find(MinorAxisLength(tmp2)>8); % 8px arbitrary
        mrk = round(local_centre(tmp2(tmp3),:));
        
        % RMK: moving windows is ok but at the moment does not get all fishing
        % weights, windows might be too large so some weights are in the same
        % zone as brigther weights (with thigh thickness + lead mrk)
        % reducing window could work but will add much noise 

    switch fig
        case 'on'
        figure;
%2016         montage(IMG,'size',[1 1]);hold on;plot(mrk(:,1),mrk(:,2),'g*')
        imshow(IMG);hold on;plot(mrk(:,1),mrk(:,2),'g*')
        plot([dx dx],ylim,'y')
        plot([dx*2 dx*2],ylim,'y')
        for i =1:4
            plot(xlim, [dy*i dy*i],'r')
        end
    end