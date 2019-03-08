function [mrk, mrk2] = f_EOS_mrkID(front,side,r_cut,test_fig)
% X.Gasparutto - Feb. 2018 - HUG

% Algorithm:
% 1 - Get radio opaque mrk position in frontal view:
    % based on a threshold in a moving window that removes everything below 
    % r_cut * Max of window
% 2 - Based on the height of the markers identified:
    % take a band of image of +/- dy around the height, process image and 
    % get mrk
% 3 - 3D position: Merge the positions of side and front view 
% 4 - Check for missing markers:
    % do 1,2,3 by starting with the side view to get markers not identified
    % on frontal view
    
% Input: front: front DICOM image from EOS imaging
%        side:  side DICOM image from EOS imaging
%        r_cut: Threshold between 0 and 1 - default 0.8
%        test_fig: option to plot or not additional figures

% Output : mrk:  all markers identified
%          mrk2: Only markers identified in part 4(Optional)           
%% 

h = waitbar(0,'Processing');

% 1 - Get Radio-Opaque mrk in frontal view
mrk_f = f_EOS_mrkID_T(front,test_fig,r_cut);

waitbar(0.25,h,'Processing');

% 2 - Marker identification based on the mrk heights from frontal view 
mrk_l = f_EOS_mrkID_H(side,mrk_f,test_fig);

waitbar(0.5,h,'Processing');

% 3 - 3D Position - Merge front and side view
cpt = 0;
for i = 1:size(mrk_f)
   % III.1 - False positive? 
   % if a mrk was not identified on lat views then false positive
   %    if mrk was not identified on lat view then lat view coordinates are 0
   if sum(mrk_l(i,:)) ~= 0 % marker
       cpt = cpt + 1;
       % III.2 - Merge 
       % Coordinate system convention:
       %  - X is lateral axis pointing left
       %  - Y is vertical axis pointing upward
       %  - Z is postero anterior axis
       % => Lateral view is (Z,Y)
       % => Frontal view is (X,Y)
       x = mrk_f(i,1);
       y = (mrk_f(i,2) + mrk_l(i,2)) / 2;
       z = mrk_l(i,1);
       %
       mrk(cpt,:) = [x y z];
   end
end

% IV - Check for missing Markers on Lateral view 
    % IV.1 - Get markers on lateral view with Threshold
    [mrk_l2] = f_EOS_mrkID_T(side,test_fig,r_cut);
    
    waitbar(0.75,h,'Processing');
    
    % IV.2 - Are markers new markers? was it already identified?
    % compare Z and Y value with other markers
    % a. Initialisation
    cpt = 0;
    d_px = 10;
    new_mrk_l =[];

    for i =1:size(mrk_l2)
        % b. Coordinate of mrk
        z = mrk_l2(i,1);
        y = mrk_l2(i,2);
        % c. Compare with identified mrk
        cmp_z = mrk_l(:,1) - z;
        cmp_y = mrk_l(:,2) - y;
        % d. Find difference < d_px on Z
        [a,~] = find( cmp_z < d_px);
        [b,~] = find( cmp_z(a) > -d_px);
        % e. Find difference < d_px on Y
        [c,~] = find( cmp_y < d_px);
        [d,~] = find( cmp_y(c) > - d_px);
        % f. If no differences < d_px on Z and Y then New Marker
        if isempty(d) && isempty(b)
            cpt = cpt + 1;
            new_mrk_l(cpt,:) = mrk_l2(i,:);
        end
        % e. Test plot
        switch test_fig
            case 'on'
            figure; plot(cmp_z,'*r');hold on;plot(cmp_y,'+b');grid minor;hold off
        end
    end
    
    % IV.3 - Identify New Markers on Frontal view
    if isempty(new_mrk_l)== 0 % If additional markers have been detected
       cpt = 0;
       [new_mrk_f] = f_EOS_mrkID_H(front,new_mrk_l,test_fig);
        mrk2 = [];
       % IV.4 - Merge
        for i = 1:size(new_mrk_f)
           % IV.4.1 - False positive? 
           % if a mrk was not identified on lat views then false positive
           %    if mrk was not identified on lat view then lat view coordinates are 0
           if sum(new_mrk_f(i,:)) ~= 0 % marker
               cpt = cpt + 1;
               % IV.4.2 - Merge 
               % Coordinate system convention:
               %  - X is lateral axis pointing left
               %  - Y is vertical axis pointing upward
               %  - Z is postero anterior axis
               % => Lateral view is (Z,Y)
               % => Frontal view is (X,Y)
               x = new_mrk_f(i,1);
               y = (new_mrk_f(i,2) + new_mrk_l(i,2)) / 2;
               z = new_mrk_l(i,1);
               %
               mrk2(cpt,:) = round([x y z]);
           end
        end
        % IV.5 - Add to the list of markers
        if isempty(mrk2) == 0
            tmp = mrk;
        mrk = round([tmp;mrk2]);
        end
    end
waitbar(1,h,'Processing');    
close (h)