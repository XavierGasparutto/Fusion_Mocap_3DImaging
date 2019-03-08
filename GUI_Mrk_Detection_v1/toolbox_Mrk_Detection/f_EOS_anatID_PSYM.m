function [psym]=f_EOS_anatID_PSYM(Front,Sagit)
% Author: v1 - X.Gasparutto (Updated) 

%% 1 - get position of symphysis on frontal view
[I_HB,I_WB] = size(Front);

% Pre-Zoom
h_IMG = 1:round(I_HB/2); 
I=Front(h_IMG,:);
fig1 = figure;hold on;
% montage(I,[1 1]);
img = imshow(I);set(fig1,'units','normalized','outerposition',[0 0 1 1])
text(100, 100,'RIGHT','color',[1 1 1],'HorizontalAlignment','left')
text(I_WB - 100, 100,'LEFT','color',[1 1 1],'HorizontalAlignment','right')
title('Select the zone around the Pubic Symphysis','FontSize',16,'FontWeight','demi');
% Select zone of the Pelvis (zoom)
[J,h_f,w_f] = f_subIMG(I); 
[J_HB,J_WB] = size(J);

% Update Image
delete(fig1.Children.Children(1:end-1)); 
img.CData = J;
x_txt = round(0.05 * J_WB);
y_txt = round(0.05 * J_HB);
text(J_WB/10  , J_HB/10,'RIGHT','color',[1 1 1],'HorizontalAlignment','left')
text(J_WB/10*9, J_HB/10,'LEFT','color',[1 1 1],'HorizontalAlignment','right')
% Get symphysis position in lateral medial direction
% get center of Symphysis
title('Click points along the pubic SYMPHYSIS','FontSize',16)
cond1 = 'No';
while strcmp(cond1,'No')
    [x,y] = getpts; symF = mean([x,y],1); 
    Ymax_f = min(y); % 0 is max
    Ymin_f = max(y); 
    % Plot 
    local_f.PSYM = symF;
    hold on; plot(local_f.PSYM(1),local_f.PSYM(2),'*r')
    plot(x,y,'+b'); clear x y 
    plot(xlim,[Ymax_f Ymax_f],'r')
    plot(xlim,[Ymin_f Ymin_f],'r')
% Was it good?
    cond1 = questdlg('Is the Lateral-Medial position right?');
    delete(fig1.Children.Children(1:4)) % Clear Lines
end


% Points in global
f_local = fieldnames(local_f);
for j =1:size(f_local,1)
    h_tmp = local_f.(f_local{j})(2) + h_f + h_IMG(1);
    w_tmp = local_f.(f_local{j})(1) + w_f;
    tmp_eos_f.(f_local{j}) = round([w_tmp h_tmp]); clear *_tmp;
end    
    
Ymin = Ymin_f + h_f + h_IMG(1);
Ymax = Ymax_f + h_f + h_IMG(1);

%% 2 - get position of symphysis on Sagittal view
clear I*
[I_HB,I_WB] = size(Sagit);

% Pre-Zoom
h_IMG_s = 1:round(I_HB/2); 
I=Sagit(h_IMG_s,:);

% Update Figure
delete(fig1.Children.Children(1:end-1)) % clear all figure elements apart image
img.CData = I; % update image
text(100, 100,'FRONT','color',[1 1 1],'HorizontalAlignment','left')
text(I_WB - 100, 100,'BACK','color',[1 1 1],'HorizontalAlignment','right')
title('Select the zone around the Pubic Symphysis','FontSize',16,'FontWeight','demi');

% plot lines from frontal image
hold on
plot([0 I_WB],[Ymin-h_IMG_s(1) Ymin-h_IMG_s(1)],'r')
plot([0 I_WB],[Ymax-h_IMG_s(1) Ymax-h_IMG_s(1)],'r')

% Select sub part of image 
[J,h_s,w_s] = f_subIMG(I);
[J_HB,J_WB] = size(J);

% Update Image
delete(fig1.Children.Children(1:end-1))
img.CData = J; 
x_txt = round(0.05 * J_WB);
y_txt = round(0.05 * J_HB);
text(x_txt, y_txt,'ANTERIOR','color',[1 1 1],'HorizontalAlignment','left')
text(J_WB - x_txt, y_txt,'POSTERIOR','color',[1 1 1],'HorizontalAlignment','right')

% plot lines from frontal image
hold on
plot([0 J_WB],[Ymin-h_IMG_s(1)-h_s Ymin-h_IMG_s(1)-h_s],'r')
plot([0 J_WB],[Ymax-h_IMG_s(1)-h_s Ymax-h_IMG_s(1)-h_s],'r')

cond2 = 'No';
while strcmp(cond2,'No')
    % Select Symphysis
    title('Click Points on the ANTERIOR part of the SYMPHYSIS','FontSize',16)
    [z,y]=getpts;
    % Fit with polynomial
    zm = mean(z); ym = mean(y);
    z_tmp = z - zm; y_tmp = y - ym;
    p3 = polyfit(y_tmp,z_tmp,3);
    yy = linspace(min(y_tmp),max(y_tmp),100);
    zz = polyval(p3,yy);
    % Most anterior point <=> dz/dt = 0
    % dz/dt = 3.p1.y^2 + 2.p2.y + p3
    % 2nd order polynomial = 0; 2 solutions
    y1 = 1/(6*p3(1)) * (-2*p3(2) + sqrt(4*p3(2)*p3(2)-12*p3(1)*p3(3)));
    y2 = 1/(6*p3(1)) * (-2*p3(2) - sqrt(4*p3(2)*p3(2)-12*p3(1)*p3(3)));
    % Select most anterior point
    if polyval(p3,y1) < polyval(p3,y2); y0 = y1; else y0 = y2; end
    z0 = polyval(p3,y0); 
    %
    PS_z = z0 + zm;
    PS_y = y0 + ym;
    % Check
    plot(z,y,'+g')
    plot(zz+zm,yy+ym,'y')
    plot(PS_z,PS_y,'ro')
    % Was the fit rigth
    cond2 = questdlg('Is the fit right?');
    delete(fig1.Children.Children(1:3))
end
close
%
local_s.PSYM = [PS_z PS_y];

% Points in global
f_local_s = fieldnames(local_s);
for j =1:size(f_local_s,1)
    h_tmp = local_s.(f_local_s{j})(2) + h_s + h_IMG_s(1);
    w_tmp = local_s.(f_local_s{j})(1) + w_s;
    tmp_eos_s.(f_local_s{j}) = round([w_tmp h_tmp]); clear *_tmp;
end    

%% Condyle center position in EOS
% x (medio-lat) is from frontal view
% y (height) is from sagittal view
% z (forward backward) is from sagittal view
psym = [tmp_eos_f.PSYM(1) tmp_eos_s.PSYM(2) tmp_eos_s.PSYM(1)];


