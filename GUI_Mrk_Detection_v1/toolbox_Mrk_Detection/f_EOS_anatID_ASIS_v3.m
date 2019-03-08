function [ASIS]= f_EOS_anatID_ASIS_v3(Front,Sagit,side,HJC)
% Author: v1 - X Gasparutto - 03.18 - HUG

% - Need to add a contrast tool

% Find Right or Left Atero Superior Iliac Spines (Radiographic definition)
% 1 - Frontal view - radiographic definition
% 2 - Sagittal view - (bit more difficult)
%
% v3 - question : do you see the medial and lateral edge of the superior
%                 part of the iliac bone?
% -> list dialog to select the method used to identify?

%% Step 1: Frontal View
% image size
[H,W]=size(Front);
% Plot 
Yimg = 1:round(H/2);
fig1 = figure; hold on
img = imshow(Front(Yimg,:));
set(fig1,'units','normalized','outerposition',[0 0 1 1])

text(100, 100,'RIGHT','color',[1 1 1],'HorizontalAlignment','left')
text(W - 100, 100,'LEFT','color',[1 1 1],'HorizontalAlignment','right')

switch side
    case 'right'
    title({'Select the zone around the Right Iliac Bone'},'FontSize',12,'FontWeight','demi');
    case 'left'
    title({'Select the zone around the Left Iliac Bone'}, 'FontSize',12,'FontWeight','demi');
end

% Selection of Sub Image - Zoom
[J,H,W] = f_subIMG(Front);
% Inverse Grayscale
% K=imcomplement(J); % inversion des nuances de gris 
% Sub Image full screen
% Update Image
delete(fig1.Children.Children(1:end-1)); 
img.CData = J; hold on
xx = xlim;
yy = ylim;
axis([xx yy])
%img.CData = K; 

% Choose the algo (2 choices so far)
q1 = questdlg('Is the medial edge of the superior iliac bone visible?');
cond1 = 'No';
while strcmp(cond1,'No')% Multiple tries if needed
    switch q1
        case 'Yes'
            inst1 = 'CLICK points on the lateral edge of the superior ilium';
            inst2 = 'CLICK points on the medial edge of the superior ilium';
        case 'No'
            inst1 = 'CLICK points on the zone between superior and inferior iliac spines';
            inst2 = 'CLICK points on the anterior part of the iliac crest';
    end
        
    % Instruction 
    title({inst1}, 'FontSize',16, 'FontWeight','demi');
    % Select zone between inferior anterior iliac spine and ASIS
    [x,y]=getpts;
    [ zM, yM, rM]=approximation_cercle(x,y);
    trace_cercle(zM,yM,rM,'g',':'); %plot(x,y,'+b')
    axis([xx yy])
    % Select anterior part of the iliac crest
    title({inst2}, 'FontSize',16, 'FontWeight','demi');
    [x2,y2]=getpts;
    [z2M,y2M,r2M]=approximation_cercle(x2,y2);
    trace_cercle(z2M,y2M,r2M,'c',':');%plot(x2,y2,'*g')
    axis([xx yy])
    % Circle intersection - work with triangle O1,O2,I w/ I = intersection point
    oM  = [ zM;   yM];
    o2M = [z2M;  y2M];
    [I] = f_intersect_circle(oM,rM,o2M,r2M);

    if isreal(I) % Circle intersection
        switch q1 
            case 'Yes'
            % Take lowest point 
                [~,idx] = max(I(:,2)); % y = 0 is top of image
            case 'No'
            % Take most lateral point
            switch side
                case 'right'
                    [~,idx] = min(I(:,1));
                case 'left'
                    [~,idx] = max(I(:,1));
            end
        end
        asi = I(idx,:);
        p_asi = plot(asi(1),asi(2),'+r');
        
        % Is it ok?
        cond1 = questdlg('Is the fit right?');
        delete(fig1.Children.Children(1:end-1))
    else
        waitfor(msgbox('Circles do not overlap'))
        % Remove lines from image
        delete(fig1.Children.Children(1:end-1))
%         for  k =size(fig1.Children.Children):-1:1
%             switch fig1.Children.Children(k).Type
%                 case 'line'
%                     delete(fig1.Children.Children(k).Type)
%             end
%         end
    end
end
% Position in EOS coordinate system
x_eos = asi(1)+W; 
y_eos = asi(2)+H; 
% y0_eos = y0 + H;
% x0_eos = x0 + W;


%% Step 2: Sagittal view
% Update Image
delete(fig1.Children.Children(1:end-1)); 
img.CData = Sagit; 
[yy,xx] = size(Sagit);
axis([1 xx 1 yy])
% Plot Hips As a guide
c_h = 'w';
trace_cercle(HJC.R(3),HJC.R(2),HJC.R(4),c_h,'-');
trace_cercle(HJC.L(3),HJC.L(2),HJC.L(4),c_h,'--');

title('Select The Pelvis and Hips Joint Centre','FontSize',16, 'FontWeight','demi')
% Select sub image
[J2,H2,W2] = f_subIMG(Sagit);

switch side
    case 'right'
        Title1 = 'click point on line that crosses the RIGHT asis';
        line = '-';
    case 'left'
        Title1 = 'click point on line that crosses the LEFT asis';
        line = '--';
end

cond2 = 'No';
c_h = 'b';
% Contrast
edgeThreshold = 0.5;
amount = 0.51;
K2 = localcontrast(J2, edgeThreshold, amount);
% Sub Image Full Screen
delete(fig1.Children.Children(1:end-1)); 
img.CData = K2; hold on
[yy,xx] = size(K2);
axis([1 xx 1 yy])

% Plot R & L HJC As a guide
trace_cercle(HJC.R(3)-W2,HJC.R(2)-H2,HJC.R(4),c_h,'-');
trace_cercle(HJC.L(3)-W2,HJC.L(2)-H2,HJC.L(4),c_h,'--');
% Notify which is which
text(HJC.R(3)-W2,HJC.R(2)-H2,{'RIGHT';'HJC'},'color',c_h,'HorizontalAlignment','Center')
text(HJC.L(3)-W2,HJC.L(2)-H2,{'LEFT';'HJC'},'color',c_h,'HorizontalAlignment','Center')
% Plot height of Superior Anterio Iliac Spine from frontal view
plot(xlim,[y_eos-H2 y_eos-H2],'-','color',c_h)
title(Title1,'FontSize',16, 'FontWeight','demi')
%
while strcmp(cond2,'No')
    [z,~] = getpts;
    z = mean(z,1);
    plot(z,y_eos-H2,'+r')
    cond2 = questdlg('Is the fit right?');
    switch cond2
        case 'Yes'
            z_eos  = z + W2;
    end
    delete(fig1.Children.Children(1)); % delete last line plotted on image
end
close
%%
ASIS = [x_eos; y_eos; z_eos];

end
