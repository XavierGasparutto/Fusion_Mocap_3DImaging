function [mrk_out, l] = f_EOS_mrkID_H(IMG,mrk_in,test_fig)
% Feb. 2018 - X.Gasparutto - HUG
% Identify markers position on one view of EOS based on the height of 
% the markers identified on the other view

mrk_out = [];

n_mrk = size(mrk_in,1);
cpt = 0;  % n mrk
dy2 = 15; % band size

[hIMG, wIMG] = size(IMG); % Height & Width IMG

% Filter
F_IMG = wiener2(IMG,[10 10]);

for i = 1:n_mrk
    % b. Band Selection 
    % take +/- dy2 around height of mrk in frontal view
    y_mrk = mrk_in(i,2)-dy2 : mrk_in(i,2)+dy2;
    % is y_mrk out of bound?
    if y_mrk(end) > hIMG % Too High
        y_mrk = y_mrk(1):hIMG;
    end
    if y_mrk(1) < 1 % Too Low
       y_mrk = 1:y_mrk(end);
    end
    l     = F_IMG(y_mrk,:); 
    % c. Image Processing
    el   = edge(l);             
    fel  = imfill(el,'holes');  
    ofel = bwareaopen(fel,180); 
    
    % d. test fig
    switch test_fig
        case 'on'
        figure;
        subplot(4,1,1);imshow(l);    title('l') %montage(l);
        subplot(4,1,2);imshow(el);   title('el') %montage(el);
        subplot(4,1,3);imshow(fel);  title('fel') %montage(fel);
        subplot(4,1,4);imshow(ofel); title('ofel') %montage(ofel)
    end
        
    % e. Mrk Identification
        % e.1 - Region properties
        rp =  regionprops(ofel);
        % e.2 - Is there Marker?
        n_mrk = size(rp,1);
        switch n_mrk 
            case 0    % No Marker
                mrk_out(i,:) = [0 0];
                % e.3 - Get local marker
            case 1    % One Marker
                mrk(i,:) = rp.Centroid;                                    % Local
                mrk_out(i,:) = [mrk(i,1) mrk(i,2)+y_mrk(1)];               % EOS Coordinate System
                switch test_fig
                    case 'on'
                     subplot(4,1,1);hold on; plot(mrk(i,1),mrk(i,2),'+r');
                     subplot(4,1,4);hold on; plot(mrk(i,1),mrk(i,2),'+r');
                end
                
            otherwise % More than one marker
                % Get marker closer to the center
                % Put y position in a vector

                % >> TESTS    
                for i = 1:n_mrk
                    tmp(i) = rp(i).Centroid(2) - dy2;
%                     bboxX(i) = rp(i).BoundingBox(1);
                    bboxY(i) = rp(i).BoundingBox(2);
%                     bboxW(i) = rp(i).BoundingBox(3);
                    bboxH(i) = rp(i).BoundingBox(4);
                end
                % Centre of the box
                tmp2 = [bboxY; bboxY + bboxH];
                tmp3 = mean(tmp2,1) - 15;
                % Distance between top image and top box should be close to
                % bottom image and bottom box, that means the shape of mrk
                % is roughly centered (NOT SURE ABOUT HYPOTHESIS)
                    tmp_a = bboxY;
                    tmp_b = 2*dy2 - (bboxY + bboxH);
                    tmp_c = abs(tmp_a - tmp_b);
                % TESTS <<
                [~,id] = min(tmp);
                [~,id] = min(tmp3);
%                 [~,id] = min(tmp_c);
                mrk(i,:) = rp(id).Centroid;
                
                mrk_out(i,:) = [mrk(i,1) mrk(i,2)+y_mrk(1)];
                switch test_fig
                    case 'on'
                        subplot(4,1,1);hold on; plot(mrk(i,1),mrk(i,2),'+r');
                        subplot(4,1,4);hold on; plot(mrk(i,1),mrk(i,2),'+r');
                end
        end
        
        % e.4 - Get mrk in EOS reference frame
        
end
