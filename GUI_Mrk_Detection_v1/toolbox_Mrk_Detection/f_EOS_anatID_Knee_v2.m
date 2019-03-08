function [cM,rM,cL,rL]=f_EOS_anatID_Knee_v2(Front,Sagit,side)
% Author: X.Gasparutto - HUG
% v2: Everything in same figure

switch side
    case 'right'
        i=1;
        Title = 'Select the zone around the RIGHT Knee';
        subLeft = 'LATERAL';
        subRight= 'MEDIAL';
    case 'left'
        i=2;
        Title = 'Select the zone around the LEFT Knee';
        subLeft = 'MEDIAL';
        subRight= 'LATERAL';
end


%% 1 - get dimensions of condyles on frontal view
[I_HB,I_WB] = size(Front);

% Pre-Zoom
h_IMG = round(I_HB/2) : round(3/4*I_HB); 
I=Front(h_IMG,:);
fig = figure(21);hold on;
img = imshow(I);
set(fig,'units','normalized','outerposition',[0 0 1 1])

text(100, 100,'RIGHT','color',[1 1 1],'HorizontalAlignment','left')
text(I_WB - 100, 100,'LEFT','color',[1 1 1],'HorizontalAlignment','right')
title(Title,'FontSize',16,'FontWeight','demi');

% Zoom on the Knee
[J,h_f,w_f] = f_subIMG(I); 
[J_HB,J_WB] = size(J);

% Update image
delete(fig.Children.Children(1:end-1)); 
img.CData = J; hold on

x_txt = round(0.05 * J_WB);
y_txt = round(0.05 * J_HB);
text(x_txt, y_txt,subLeft,'color',[1 1 1],'HorizontalAlignment','left')
text(J_WB - x_txt, y_txt,subRight,'color',[1 1 1],'HorizontalAlignment','right')

cond1 = 'No';
while strcmp(cond1,'No')
% Get condyle position in lateral medial direction
    % Medial Condyle
    % get Medial part of Medial Condyle
    title('Click points on MEDIAL part of MEDIAL Condyle','FontSize',16)
    [x,y] = getpts; medM = mean([x,y],1); clear x y 
    % get Lat    part of Medial Condyle
    title('Click points on LATERAL part of MEDIAL Condyle','FontSize',16)
    [x,y] = getpts; medL = mean([x,y],1); clear x y 

        % Center of medial condyle
        % -> The height computed here is irrelevant, use sphere fitting in
        % sagittal view to get it
        local_f.medC = (medM + medL) / 2;
%         hold on; 
        plot(local_f.medC(1),local_f.medC(2),'*r')

    % get Medial part of Medial Condyle
    title('Click points on MEDIAL part of LATERAL Condyle','FontSize',16)
    [x,y] = getpts; latM = mean([x,y],1); clear x y 
    % get Lat    part of Medial Condyle
    title('Select points on LATERAL part of LATERAL Condyle','FontSize',16)
    [x,y] = getpts; latL = mean([x,y],1); clear x y 

        % Center of medial condyle
        % -> The height computed here is irrelevant, use sphere fitting in
        % sagittal view to get it
        local_f.latC = (latM + latL) / 2;
%         hold on; 
        plot(local_f.latC(1),local_f.latC(2),'*c')
    
% Get Band to help identification of Right left, medial lateral on sagittal image
% NOTE: min is for top Part as point (0,0) is top left of image
    
    % get Top part of Medial Condyle
    title('Click points on TOP part of MEDIAL Condyle','FontSize',16)
    [x,y] = getpts; 
    if size(x,1) > 1 ; local_f.medT = min([x,y]); else; local_f.medT = [x,y]; end
    clear x y; plot(xlim,[local_f.medT(2) local_f.medT(2)],'r') 
    % get Bottom Part of Medial Condyle
    title('Click points on BOTTOM part of MEDIAL Condyle','FontSize',16)
    [x,y] = getpts; 
    if size(x,1) > 1 ; local_f.medB = max([x,y]); else; local_f.medB = [x,y]; end
    clear x y; plot(xlim,[local_f.medB(2) local_f.medB(2)],'--r') 

    % get Top part of Lateral Condyle
    title('Click points on the TOP part of MEDIAL Condyle','FontSize',16)
    [x,y] = getpts; 
    if size(x,1) > 1 ; local_f.latT = min([x,y]); else; local_f.latT = [x,y]; end
    clear x y; plot(xlim,[local_f.latT(2) local_f.latT(2)],'c') 
    % get Bottom Part of Lateral Condyle
    title('Click points on the BOTTOM part of MEDIAL Condyle','FontSize',16)
    [x,y] = getpts; 
    if size(x,1) > 1 ; local_f.latB = max([x,y]); else; local_f.latB = [x,y]; end
    clear x y; plot(xlim,[local_f.latB(2) local_f.latB(2)],'--c') 
    
% Was it good?
    cond1 = questdlg('Is the identification right?');
    delete(fig.Children.Children(1:6))
end

% Points in global
f_local = fieldnames(local_f);
for j =1:size(f_local,1)
    h_tmp = local_f.(f_local{j})(2) + h_f + h_IMG(1);
    w_tmp = local_f.(f_local{j})(1) + w_f;
    tmp_eos_f.(f_local{j}) = round([w_tmp h_tmp]); clear *_tmp;
end    
    
%% 2 - get dimensions of condyles on Sagittal view
clear I*
[I_HB,I_WB] = size(Sagit);

% Pre-Zoom
h_IMG_s = round(I_HB/2) : round(3/4*I_HB); 
I=Sagit(h_IMG_s,:);

% Update image
delete(fig.Children.Children(1:end-1)); 
img.CData = I; hold on
text(100, 100,'FRONT','color',[1 1 1],'HorizontalAlignment','left')
text(I_WB - 100, 100,'BACK','color',[1 1 1],'HorizontalAlignment','right')
title(Title,'FontSize',16,'FontWeight','demi');

% plot lines from frontal image
plot(xlim,[tmp_eos_f.medT(2)-h_IMG_s(1) tmp_eos_f.medT(2)-h_IMG_s(1)],'r')
plot(xlim,[tmp_eos_f.medB(2)-h_IMG_s(1) tmp_eos_f.medB(2)-h_IMG_s(1)],'--r')
plot(xlim,[tmp_eos_f.latT(2)-h_IMG_s(1) tmp_eos_f.latT(2)-h_IMG_s(1)],'c')
plot(xlim,[tmp_eos_f.latB(2)-h_IMG_s(1) tmp_eos_f.latB(2)-h_IMG_s(1)],'--c')

% Zoom on Condyles 
[J,h_s,w_s] = f_subIMG(I);
[J_HB,J_WB] = size(J);

% Update image
delete(fig.Children.Children(1:end-1)); 
img.CData = J; hold on
x_txt = round(0.05 * J_WB);
y_txt = round(0.05 * J_HB);
text(x_txt, y_txt,'ANTERIOR','color',[1 1 1],'HorizontalAlignment','left')
text(J_WB - x_txt, y_txt,'POSTERIOR','color',[1 1 1],'HorizontalAlignment','right')

% Plot lines from frontal image / Height of condyles
hold on
plot(xlim,[tmp_eos_f.medT(2)-h_IMG_s(1)-h_s tmp_eos_f.medT(2)-h_IMG_s(1)-h_s],'r')
plot(xlim,[tmp_eos_f.latT(2)-h_IMG_s(1)-h_s tmp_eos_f.latT(2)-h_IMG_s(1)-h_s],'c');
plot(xlim,[tmp_eos_f.medB(2)-h_IMG_s(1)-h_s tmp_eos_f.medB(2)-h_IMG_s(1)-h_s],'--r')
plot(xlim,[tmp_eos_f.latB(2)-h_IMG_s(1)-h_s tmp_eos_f.latB(2)-h_IMG_s(1)-h_s],'--c')
legend('MEDIAL','LATERAL','Location','southeast');

cond2 = 'No';
while strcmp(cond2,'No')
    % Select Condyles
    % Medial Condyle
    title('Click Points on the POSTERIOR part of the MEDIAL condyle','FontSize',16)
    [z,y]=getpts;
    % Détermination des cercles des moindres carrés passant par ces points 
    [zM,yM,rM]=approximation_cercle(z,y);
    % Check
    plot(z,y,'+g')
    trace_cercle(zM,yM,rM,'r','-');
    % Was the fit rigth
    cond2 = questdlg('Is the fit right?');
    switch cond2
        case 'No'
            delete(fig.Children(2).Children(1:3))
    end
end

% Lateral Condyle
title('Click points on the POSTERIOR part of the LATERAL condyle','FontSize',16)

cond3 = 'No';
while strcmp(cond3,'No')
    [z,y]=getpts;
    % Détermination des cercles des moindres carrés passant par ces points 
    [zL,yL,rL]=approximation_cercle(z,y);
    % Check
    plot(z,y,'+y')
    trace_cercle(zL,yL,rL,'c','-');
    % Was the fit rigth
    cond3 = questdlg('Is the fit right?');
    switch cond3
        case 'No' % Plot figure rm circle + points
             delete(fig.Children(2).Children(1:3))
    end
end
close

% Save Points of Condyles
local_s.medC = [zM yM];
local_s.latC = [zL yL];

% Condyles in global
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
cM = [tmp_eos_f.medC(1) tmp_eos_s.medC(2) tmp_eos_s.medC(1)];
rM = round(rM);
cL = [tmp_eos_f.latC(1) tmp_eos_s.latC(2) tmp_eos_s.latC(1)];
rL = round(rL);
