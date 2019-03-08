% Author: Xavier Gasparutto 
% Institution: Willy Taillard Kinesiology Laboratory, HUG, UNIGE, Geneva
% Release: March 2019
% This GUI needs the image processing toolbox of matlab
% it was developped on Matlab R2016b and tested on R2015b and R2018a

% -------------------------------------------------------------------------
% This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 
% International License. To view a copy of this license, visit 
% http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to 
% Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
% -------------------------------------------------------------------------

function varargout = GUI_Mrk_Detection(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_Mrk_Detection_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_Mrk_Detection_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

end

% --- Executes just before GUI_Mrk_Detection is made visible.
function GUI_Mrk_Detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_Mrk_Detection (see VARARGIN)

addpath('toolbox_Mrk_Detection')

% Choose default command line output for GUI_Mrk_Detection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
if nargin == 3
    initial_dir = pwd;
elseif nargin > 4
    if strcmpi(varargin{1},'dir')
        if exist(varargin{2},'dir')
            initial_dir = varargin{2};
        else
            errordlg('Input argument must be a valid directory','Input Argument Error!')
            return
        end
    else
        errordlg('Unrecognized input argument','Input Argument Error!');
        return;
    end
end
end

% --- Outputs from this function are returned to the command line.
function varargout = GUI_Mrk_Detection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

%% Load DICOMs
% --- Executes on button press in loadfile.
function loadfile_Callback(hObject, eventdata, handles)
% hObject    handle to loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dir_path files Front Sagit AnatPt Info
% Clear figure tables and previous variables
AnatPt = [];
% Initialise tables
set(handles.Table_SkinMrk,'Data',[])
set(handles.Table_AnatPt,'Data',[])
for i =1:24; R{i} = ['Mrk',num2str(i)];end
set(handles.Table_SkinMrk,'RowName',R)

cla(handles.axes3)

[files,dir_path] = uigetfile('data\*.*','Select Frontal and Sagittal view','MultiSelect','on');

% Check if both files have been selected
if size(files,2) == 2
    % Edit Text
    % Directory
    set(handles.Directory,'string',dir_path)
    % Files
    set(handles.FrontFile,'string',files{1})
    set(handles.SagittalFile,'string',files{2})

    % Load Image
    h = waitbar(0,'Processing');
    % set(handles.axes3,'Selected','on')
    Front = dicomread([dir_path,'\',files{1}]);
    Info.Front = dicominfo([dir_path,'\',files{1}]);
    waitbar(0.5,h,'Processing');
    Sagit  = dicomread([dir_path,'\',files{2}]);
    Info.Sagit = dicominfo([dir_path,'\',files{1}]);
    waitbar(1,h,'Processing');
    close(h)
    imshow([Front Sagit])%montage([Front Sagit],'size',[1 1]) % Ok in 2016a
    
    % Plot X Y Z
    W = size(Front,2);
    r = size(Sagit,2)/ 5;
    hold on;
    font = 16;
    quiver(1,1,r,0,'linewidth',2,'color','r','MaxHeadSize',2); text(r/3,r/8,'x','color','r','FontSize',font) % X
    quiver(1,1,0,r,'linewidth',2,'color','g','MaxHeadSize',2); text(r/8,3*r/5,'y','color','g','FontSize',font) % Y
    quiver(W + 1,1,r,0,'linewidth',2,'color','b','MaxHeadSize',2); text(W + r/3,r/8,'z','color','b','FontSize',font) % Z
    quiver(W + 1,1,0,r,'linewidth',2,'color','g','MaxHeadSize',2); text(W + r/8,3*r/5,'y','color','g','FontSize',font) % Y
    % set(handles.axes3,'Selected','off')
else 
    msgbox('Select the FRONTAL and SAGITTAL files')
end
end

%% Skin Markers Detection
% --- Executes on button press in SkinMrk.
function SkinMrk_Callback(hObject, eventdata, handles)
% hObject    handle to SkinMrk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Front Sagit Info

r_cut = str2double(handles.r_cut.String); % add selection in GUI

test_fig = 'off'; % on / off
 
% Clear Markers in Figure 
axes(handles.axes3)
h_tmp = findobj('type','line'); % Find Markers on plot
delete(h_tmp(:));               % Delete Markers on plot

% Initialise Row Name
for i = 1:size(handles.Table_SkinMrk.RowName,1)
    Row{i} = ['Mrk',num2str(i)];
end
set(handles.Table_SkinMrk,'RowName',Row)
set(handles.Table_SkinMrk,'Data',[])

[mrk_eos,mrk_eos2] = f_EOS_mrkID(Front,Sagit,r_cut,test_fig);

% Plot MRK -> In Fig GUI
xdim = size(Front,2);
hold on
% Add Markers
    x = mrk_eos(:,1);
    y = mrk_eos(:,2);
    z = mrk_eos(:,3);
% Front view (X,Y)
    plot(x,y,'*g','MarkerSize',5)
% Lateral view (Z,Y)
    plot(z + xdim,y,'*g','MarkerSize',5)
    
% Markers identified on side view
if exist('mrk_eos2')
    if isempty(mrk_eos2)== 0 % If additional markers have been detected on side view
    clear x y z
        x = mrk_eos2(:,1);
        y = mrk_eos2(:,2);
        z = mrk_eos2(:,3);
    % Front view (X,Y)
        plot(x,y,'*r','MarkerSize',5)
    % Lateral view (Z,Y)
        plot(z + xdim,y,'*r','MarkerSize',5)
    end
end

% Fill Marker Table with values
% Mrk Name

% Table_SkinMrk
px2mm = Info.Sagit.PixelSpacing (1);
mrk_eos_mm = mrk_eos * px2mm; % EOS to mm

% mrk_eos_mm(2) = ; 
set(handles.Table_SkinMrk,'Data',mrk_eos_mm);

% Missing Data?
% Are there missing markers? no -> do nothing, yes -> next step
cond1 = questdlg('Are there unidentified markers?');
if strcmp(cond1,'Yes')
    % 1 - Plot Frontal view with small black band on sides to see markers at edge
        % 1.1 - Add black band to see markers at the edge
        dimY = size(Front,1); 
        dx = 50;
        tmpIMG = [zeros(dimY,dx),Front,zeros(dimY,dx)];
        % 1.2 - Height and Width boundaries of image - for rectangle
        [HB,WB] = size(tmpIMG); 
        % 1.3 - Plot frontal view
        select_mrk = figure; imshow(tmpIMG); hold on;
        set(select_mrk,'units','normalized','outerposition',[0 0 1 1]);
        % 1.4 - Plot markers identified
        plot(mrk_eos(:,1)+dx,mrk_eos(:,2),'+g','MarkerSize',5)
        title('Select Zone Around the Mrk')
        %
    n_mrk = size(mrk_eos,1);
    cpt = 0;
    while strcmp(cond1,'Yes') % ask question after each mrk
    
    % 2 - Select unidentified Marker on Frontal View
        % 2.1 - Manual ID of Rectangle around marker 
        clear mrk_tmp
        figure(select_mrk.Number)
        cpt = cpt +1;
        [subIMG,h,w] = f_subIMG(tmpIMG);
        
        % 2.2 - ID mrk with threshold on sub image
        r_cut2 = 0.9;
        local_mrk(cpt,:) = f_EOS_mrkID_subT(subIMG,r_cut2);
        
            % >> TO DO -- Need to add mrk that were previously identified
            % on subIMG, happens quite often <<
            
            % 2.3.1 - Check visually if it worked
            figure; imshow(subIMG); hold on; 
            plot(local_mrk(cpt,1),local_mrk(cpt,2),'*r')
            % Successful? Y/N
            cond2 = questdlg('Was the marker identified properly?');
            close
            
        % 2.4 - Get Mrk X and Y Position in Global
            mrk_tmp(cpt,1:2) = round([local_mrk(cpt,1) + round(w) - dx,...
                               local_mrk(cpt,2) + round(h)]); 
            
        % 2.5 - Is this New Marker?
        % marker is not new if it is within a box of 10px in X and Y around other mrk
        x = [mrk_tmp(cpt,1) - 10, mrk_tmp(cpt,1) + 10];
        y = [mrk_tmp(cpt,2) - 10, mrk_tmp(cpt,2) + 10];
        
        newmrk = find(mrk_eos(:,1)> x(1) & mrk_eos(:,1)< x(2)  ...
                  & mrk_eos(:,2)> y(1) & mrk_eos(:,2)< y(2));
        
        if isempty(newmrk) % New Marker
            
    % 3 - Get Marker Position on Lateral View (If Work)
            if strcmp(cond2,'Yes')
                clear mrk_band
                % 3.1 - Try Band - Add Height constraint on Lat view
                [mrk_eos3,band] = f_EOS_mrkID_H(Sagit,mrk_tmp(cpt,:),'off');
                % Band dimension
                [hb,~] = size(band);

                % 3.1' - If Band did not work - Select Marker Manually
                if sum(mrk_eos3) == 0 % band did not work
                    figure; 
                    imshow(band); title('Select Zone Around the Mrk');
                    % 3.1'.1 - Mrk Selection
                    [subBand,~,w2] = f_subIMG(band);
                    % 3.1'.2 -ID marker with Threshold
                    mrk_subBand = f_EOS_mrkID_subT(subBand,r_cut);
                    mrk_band = [round(mrk_subBand(1) + w2) , round(mrk_subBand(2))];
                    % 3.1'.4 -Mrk in EOS
                    mrk_eos3 = [mrk_band(1) ,round((mrk_band(2) + mrk_tmp(cpt,2) - hb/2))];

                else % Band did work - mrk in local for visualisation
                    mrk_band = round( [mrk_eos3(1), mrk_eos3(2) - mrk_tmp(cpt,2) + hb /2 ]);

                end

                % 3.2 - Check If marker ID is correct
                figure; 
                imshow(band); hold on;plot(mrk_band(:,1),mrk_band(:,2),'+r')
                cond3 = questdlg('Was the marker identified properly?');
                close
                %
                if strcmp(cond3,'Yes') % marker identified on both views
                    % Get Position in EOS
                    x = mrk_tmp(cpt,1);
                    y = round((mrk_tmp(cpt,2) + mrk_eos3(2))/2);
                    z = round(mrk_eos3(1));
                    %
                    mrk_eos(n_mrk+ cpt,:) = [x y z];
                else % Band worked but did not get the right mrk 
                    figure; 
                    imshow(band); hold on;
                     % Select zone around marker
                     [subBand,~,w2] = f_subIMG(band);
                     % 3.1'.2 -ID marker with Threshold
                     mrk_subBand = f_EOS_mrkID_subT(subBand,r_cut);
                     mrk_band = [round(mrk_subBand(1) + w2) , round(mrk_subBand(2))];
                     % 3.1'.4 -Mrk in EOS
                     mrk_eos3 = [mrk_band(1) ,round((mrk_band(2) + mrk_tmp(cpt,2) - hb/2))];
                     plot(mrk_band(:,1),mrk_band(:,2),'+r')
                     %
                     cond4 = questdlg('Was the marker identified properly?');
                     if strcmp(cond4,'Yes')
                         close
                         % Get Position in EOS
                         x = mrk_tmp(cpt,1);
                         y = round((mrk_tmp(cpt,2) + mrk_eos3(2))/2);
                         z = round(mrk_eos3(1));
                         mrk_eos(n_mrk+ cpt,:) = [x y z];
                     end
                end 
                % 3.3 - Update Figures & Table
                % Update Figure: mrk ID
                    figure(select_mrk.Number); hold on
                    plot(mrk_eos(end,1) + dx,mrk_eos(end,2),'*r')
                % Update Figure GUI
                    axes(handles.axes3)      
                    plot(mrk_eos(end,1),mrk_eos(end,2),'*c')
                    plot(mrk_eos(end,3) + xdim,mrk_eos(end,2),'*c')
                % Update Table
                    mrk_eos_mm = mrk_eos * px2mm;
                    set(handles.Table_SkinMrk,'Data',mrk_eos_mm);
            end                 
            cond1 = questdlg('Are there unidentified markers?');
        
        else % Not A new Marker
            cond1 = questdlg({'Marker Already Identified';'';'Are there unidentified markers?'});
            switch cond1
                case 'No'
                    close
            end
        end
    end
else % Do Nothing

end
end

%% Anatomical Points
% --- Executes on button press in RightHip.
function RightHip_Callback(hObject, eventdata, handles)
% hObject    handle to RightHip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Front Sagit AnatPt Info
[r_HJC, r_R] = f_EOS_anatID_Hip_v3(Front,Sagit,'right');

% Plot On fig
xdim = size(Front,2);
axes(handles.axes3); hold on
trace_cercle(r_HJC(1),r_HJC(2),r_R,'c','-')
trace_cercle(r_HJC(3)+xdim,r_HJC(2),r_R,'c','-')

% Put in Table
px2mm = Info.Sagit.PixelSpacing (1);
% AnatPt(1,1:4) = round([r_HJC r_R] * px2mm);
AnatPt(1,1:4) = [r_HJC r_R] * px2mm;
set(handles.Table_AnatPt,'Data',AnatPt);
end

% --- Executes on button press in LeftHip.
function LeftHip_Callback(hObject, eventdata, handles)
% hObject    handle to LeftHip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Front Sagit AnatPt Info

[l_HJC, l_R] = f_EOS_anatID_Hip_v3(Front,Sagit,'left');
% Plot On fig
xdim = size(Front,2);
axes(handles.axes3); hold on
trace_cercle(l_HJC(1),l_HJC(2),l_R,'c','--')
trace_cercle(l_HJC(3)+xdim,l_HJC(2),l_R,'c','--')

% Put in Table
px2mm = Info.Sagit.PixelSpacing (1);
% AnatPt(2,1:4) = round([l_HJC l_R] * px2mm);
AnatPt(2,1:4) = [l_HJC l_R] * px2mm;
set(handles.Table_AnatPt,'Data',AnatPt);
end

% --- Executes on button press in RightKnee.
% Knee algo was updated X.Gasparutto 2018
function RightKnee_Callback(hObject, eventdata, handles)
% hObject    handle to RightKnee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Front Sagit AnatPt Info

[r_cM, r_rcM, r_cL, r_rcL] = f_EOS_anatID_Knee_v2(Front,Sagit,'right');

% Plot on main figure
% xdim = handles.axes3.XLim(2) / 2;
axes(handles.axes3);hold on
xdim = size(Front,2);
% Frontal Plane
trace_cercle(r_cM(1),r_cM(2),r_rcM,'y','-');
trace_cercle(r_cL(1),r_cL(2),r_rcL,'g','-');
% Sagittal Plane
trace_cercle(r_cM(3)+xdim, r_cM(2),r_rcM,'y','-');
trace_cercle(r_cL(3)+xdim, r_cL(2),r_rcL,'g','-');

% Add to table
px2mm = Info.Sagit.PixelSpacing (1);
% AnatPt(3,1:4) = round([r_cM r_rcM] * px2mm);
% AnatPt(4,1:4) = round([r_cL r_rcL] * px2mm);
AnatPt(3,1:4) = [r_cM r_rcM] * px2mm;
AnatPt(4,1:4) = [r_cL r_rcL] * px2mm;
set(handles.Table_AnatPt,'Data',AnatPt);
end

% --- Executes on button press in LeftKnee.
function LeftKnee_Callback(hObject, eventdata, handles)
% hObject    handle to LeftKnee (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Front Sagit AnatPt Info

[l_cM, l_rcM, l_cL, l_rcL] = f_EOS_anatID_Knee_v2(Front,Sagit,'left');

% Plot on main figure
xdim = size(Front,2);
axes(handles.axes3);hold on
% Frontal Plane
trace_cercle(l_cM(1),l_cM(2),l_rcM,'y','--');
trace_cercle(l_cL(1),l_cL(2),l_rcL,'g','--');
% Sagittal Plane
trace_cercle(l_cM(3)+xdim,l_cM(2),l_rcM,'y','--');
trace_cercle(l_cL(3)+xdim,l_cL(2),l_rcL,'g','--');

% Add to table
px2mm = Info.Sagit.PixelSpacing (1);
% AnatPt(5,1:4) = round([l_cM l_rcM] * px2mm);
% AnatPt(6,1:4) = round([l_cL l_rcL] * px2mm);
AnatPt(5,1:4) = [l_cM l_rcM] * px2mm;
AnatPt(6,1:4) = [l_cL l_rcL] * px2mm;
set(handles.Table_AnatPt,'Data',AnatPt);
end

% --- Executes on button press in RASIS.
function RASIS_Callback(hObject, eventdata, handles)
% hObject    handle to RASIS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Front Sagit Info AnatPt

% Where Anat Point identified?
if isempty(handles.Table_AnatPt.Data) == 0
    % Load Hips Position & put in pixels
    mm2px = 1 / Info.Sagit.PixelSpacing (1);
%     HJC.R = round(handles.Table_AnatPt.Data(1,:) * mm2px);
%     HJC.L = round(handles.Table_AnatPt.Data(2,:) * mm2px);
    HJC.R = handles.Table_AnatPt.Data(1,:) * mm2px;
    HJC.L = handles.Table_AnatPt.Data(2,:) * mm2px;
    
    % Where Hip Identified? Need it for the lateral view
    if (sum(HJC.R) == 0 || sum(HJC.L) == 0) == 0

        [RASIS]= f_EOS_anatID_ASIS_v3(Front,Sagit,'right',HJC);
        
        % Plot on main figure
        xdim = size(Front,2);
        axes(handles.axes3);hold on
        % Frontal Plane
        plot(RASIS(1),RASIS(2),'vk','MarkerFaceColor',[255 255 0]/ 255)
        % Sagittal Plane
        plot(RASIS(3)+ xdim,RASIS(2),'vk','MarkerFaceColor',[255 255 0]/ 255)
        % Add to table
        px2mm = Info.Sagit.PixelSpacing (1);
%         AnatPt(8,1:3) = round(RASIS * px2mm);
        AnatPt(8,1:3) = RASIS * px2mm;
        set(handles.Table_AnatPt,'Data',AnatPt);
    else
        msgbox('Identify Right and Left Hip Joint Centres First')

    end
else
    msgbox('Identify Right and Left Hip Joint Centres First')
end
%
% Add To Table - TO DO
end

% --- Executes on button press in LASIS.
function LASIS_Callback(hObject, eventdata, handles)
% hObject    handle to LASIS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Front Sagit Info AnatPt

% Where Anat Points identified? 
if isempty(handles.Table_AnatPt.Data) == 0
    % Load Hips Position & put in pixels
    mm2px = 1 / Info.Sagit.PixelSpacing (1);
%     HJC.R = round(handles.Table_AnatPt.Data(1,:) * mm2px);
%     HJC.L = round(handles.Table_AnatPt.Data(2,:) * mm2px);
    HJC.R = handles.Table_AnatPt.Data(1,:) * mm2px;
    HJC.L = handles.Table_AnatPt.Data(2,:) * mm2px;
    % Where Hip Identified? Need it for the lateral view
    if (sum(HJC.R) == 0 || sum(HJC.L) == 0) == 0

        [LASIS]= f_EOS_anatID_ASIS_v3(Front,Sagit,'left',HJC);
        
    % Plot on main figure
        xdim = size(Front,2);
        axes(handles.axes3);hold on
        % Frontal Plane
        plot(LASIS(1),LASIS(2),'^k','MarkerFaceColor',[255 255 0]/ 255)
        % Sagittal Plane
        plot(LASIS(3)+ xdim,LASIS(2),'^k','MarkerFaceColor',[255 255 0]/ 255)
    
    % Add to table
        px2mm = Info.Sagit.PixelSpacing (1);
%         AnatPt(9,1:3) = round(LASIS * px2mm);
        AnatPt(9,1:3) = LASIS * px2mm;
        set(handles.Table_AnatPt,'Data',AnatPt);    
    else
        msgbox('Identify Right and Left Hips First')
    end
else
    msgbox('Identify Right and Left Hip Joint Centres First')
end

end

% --- Executes on button press in PSymphysis.
function PSymphysis_Callback(hObject, eventdata, handles)
% hObject    handle to PSymphysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Front Sagit Info AnatPt

% 1 - Identify PSYM
[PSYM]=f_EOS_anatID_PSYM(Front,Sagit);

% 2 - Plot PSYM on main figure
    xdim = size(Front,2);
    axes(handles.axes3);hold on
    % Frontal Plane
    plot(PSYM(1),PSYM(2),'dk','MarkerFaceColor',[255 255 0]/ 255)
    % Sagittal Plane
    plot(PSYM(3)+ xdim,PSYM(2),'dk','MarkerFaceColor',[255 255 0]/ 255)
% 3 - Write in Table 
px2mm = Info.Sagit.PixelSpacing (1);
% AnatPt(7,1:3) = round(PSYM * px2mm);
AnatPt(7,1:3) = PSYM * px2mm;
set(handles.Table_AnatPt,'Data',AnatPt);
end

%% Labelling Skin Markers
% --- Executes on button press in Manual_Label.
function Manual_Label_Callback(hObject, eventdata, handles)
% hObject    handle to Manual_Label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Front Info
px2mm = Info.Sagit.PixelSpacing (1);
% Get All points
xyz = handles.Table_SkinMrk.Data;
% Sort by height
[~,id_x] = sort(xyz(:,2));
xyz_s = xyz(id_x,:);
% Identify from top to bottom
h = figure;
scr_size = get(groot,'ScreenSize');
h.Position = [scr_size(3)/6 1 scr_size(3)/3 scr_size(4)/10*9];
imshow(Front);hold on
% add left & right
text(10,100,'Right','Color','w','FontSize',14)
text(size(Front,2)-10,100,'Left','Color','w','FontSize',14,'HorizontalAlignment','Right')
% Put image on right part of screen
plot(xyz_s(:,1)/px2mm,xyz_s(:,2)/px2mm,'r+')
b = plot(xyz_s(1,1)/px2mm,xyz_s(1,2)/px2mm,'co');
tmp = inputdlg('Name of marker?  ');
C{1} = tmp{1};
for i =2:size(xyz_s,1)
    b.XData = xyz_s(i,1)/px2mm;
    b.YData = xyz_s(i,2)/px2mm;
    tmp = inputdlg('Name of marker?  ');
    C{i} = tmp{1};
end
close

set(handles.Table_SkinMrk,'RowName',C)
set(handles.Table_SkinMrk,'Data',xyz_s)

end

%% EXPORT
% --- Executes on button press in ExportMAT.
function ExportMAT_Callback(hObject, eventdata, handles)
% hObject    handle to ExportMAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dir_path files

% 1 - Get Data
    % 1.1 - Skin Markers
    % Load from table
    data_tmp = handles.Table_SkinMrk.Data;
    row_tmp = handles.Table_SkinMrk.RowName;

    if ~isempty(data_tmp)
        % Structure form
        for i = 1:size(data_tmp,1)
            SkinMrk.([row_tmp{i}]) = data_tmp(i,:);
        end
        clear *_tmp
    end
    
    % 1.2 - Anatomical Points 
    % Load from table
    data_tmp = handles.Table_AnatPt.Data;
    row_tmp  = handles.Table_AnatPt.RowName;
    
    % Structure form
    if ~isempty(data_tmp)
        for i = 1:size(data_tmp,1)
            AnatPt.([row_tmp{i}]) = data_tmp(i,1:4);
        end
        clear *_tmp
    end
    
% 2 - Save File
    % 2.1 - File Name - Work if EOS files are saved as 'project_number_visit_*' 
        tmp1 = files{1};
        tmp2 = strfind(tmp1,'_');
    op = upper(handles.Operator.String); % Operator 
    file_name = tmp1(1:tmp2(3)-1); 
    % 2.2 - Is there an export file already?
    if isempty(op)
        tmp3 = dir([dir_path,files{1}(1:8),'*.mat']);
    else
        tmp3 = dir([dir_path,files{1}(1:8),'*',op,'*.mat']);
    end
    %
    if isempty(tmp3)
            idx = 1; 
    else
        for k=1:size(tmp3,1)
            tmp5 = strfind(tmp3(k).name,'_');
            tmp6 = strfind(tmp3(k).name,'.');
            tmp7(k) = str2num(tmp3(k).name(tmp5(end)+1:tmp6-1));
        end
        idx = max(tmp7) + 1; 
    end
    % save data to directory of DICOM files
    if isempty(op)
        out_name = [file_name,'_mrkEOS_',num2str(idx),'.csv'];
    else
        out_name = [file_name,'_mrkEOS_',op,'_',num2str(idx),'.csv'];       
    end
    
    if exist('SkinMrk','var')
        if exist('AnatPt','var')
            save([dir_path,out_name],'SkinMrk','AnatPt')
            msgbox('Data was succesfully exported')
        else
            save([dir_path,out_name],'SkinMrk')
            msgbox('Data was succesfully exported')
        end
    else
        if exist('AnatPt')
            save([dir_path,out_name],'AnatPt')
            msgbox('Data was succesfully exported')
        else
            msgbox('No data to export')
        end
    end
end

% --- Executes on button press in ExportCSV.
function ExportCSV_Callback(hObject, eventdata, handles)
% hObject    handle to ExportCSV (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dir_path files    
    
% 1 - Get Data    
    % 1.1 - Skin Markers
    % Load from table
    if ~isempty(handles.Table_SkinMrk.Data)
        data_tmp = handles.Table_SkinMrk.Data;
        row_tmp = handles.Table_SkinMrk.RowName;
        % Prepare Table form
        for i = 1:size(data_tmp,1)
            X(i,1) = data_tmp(i,1);
            Y(i,1) = data_tmp(i,2);
            Z(i,1) = data_tmp(i,3);
            R(i,1) = 0;
            row{i} = [row_tmp{i}];
        end
        clear *_tmp
    else 
        i = 0;
    end
    
    % 1.2 - Anatomical Points
    % Load from table
    if ~isempty(handles.Table_AnatPt.Data)
        data_tmp = handles.Table_AnatPt.Data;
        row_tmp = handles.Table_AnatPt.RowName;
        % Prepare Table form
        for j = 1:size(data_tmp,1)
            X(j+i,1) = data_tmp(j,1);
            Y(j+i,1) = data_tmp(j,2);
            Z(j+i,1) = data_tmp(j,3);
            R(j+i,1) = data_tmp(j,4);
            row{j+i} = row_tmp{j};
        end
        clear *_tmp
    end
    
    % Were there false positive in mrk detection?
    tmp = strcmp(row,'');
    if isempty(tmp) == 1  % Fake Mrk Identified
        row{strcmp(row,'')} = 'notMrk';
    end

    % Export
    if exist('X','var')
        T = table(X,Y,Z,R,'RowNames',row);
        % File Name - Work if EOS files are saved as 'project_number_visit_*' 
            tmp1 = files{1};
            tmp2 = strfind(tmp1,'_');
            op = upper(handles.Operator.String); % Operator
        file_name = tmp1(1:tmp2(3)-1); 
        % Is there an export file already?
        if isempty(op)
            tmp3 = dir([dir_path,files{1}(1:8),'*.mat']);
        else
            tmp3 = dir([dir_path,files{1}(1:8),'*',op,'*.mat']);
        end
        %
        if isempty(tmp3)
            idx = 1; 
        else
            for k=1:size(tmp3,1)
                tmp5 = strfind(tmp3(k).name,'_');
                tmp6 = strfind(tmp3(k).name,'.');
                tmp7(k) = str2num(tmp3(k).name(tmp5(end)+1:tmp6-1));
            end
            idx = max(tmp7) + 1; 
        end
        % Save Data
        if isempty(op)
            out_name = [file_name,'_mrkEOS_',num2str(idx),'.csv'];
        else
            out_name = [file_name,'_mrkEOS_',op,'_',num2str(idx),'.csv'];       
        end
        writetable(T,[dir_path,out_name],'WriteRowNames',true)
        msgbox('Data was succesfully exported')
    else
        msgbox('No data to export')
    end
end


% --- Executes on button press in ExportXLS.
function ExportXLS_Callback(hObject, eventdata, handles)
% hObject    handle to ExportXLS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dir_path files    
    
% 1 - Get Data    
    % 1.1 - Skin Markers
    % Load from table
    data_tmp = handles.Table_SkinMrk.Data;
    row_tmp = handles.Table_SkinMrk.RowName;
    if ~isempty(data_tmp)
        % Prepare Table form
        for i = 1:size(data_tmp,1)
            X(i,1) = data_tmp(i,1);
            Y(i,1) = data_tmp(i,2);
            Z(i,1) = data_tmp(i,3);
            R(i,1) = 0;
            row{i} = [row_tmp{i}];
        end
        clear *_tmp
    else 
        i = 0;
    end 
    
    % 1.2 - Anatomical Points 
    % Load from table
    data_tmp = handles.Table_AnatPt.Data;
    row_tmp = handles.Table_AnatPt.RowName;
    if ~isempty(data_tmp)
        % Prepare Table form
        for j = 1:size(data_tmp,1)
            X(j+i,1) = data_tmp(j,1);
            Y(j+i,1) = data_tmp(j,2);
            Z(j+i,1) = data_tmp(j,3);
            R(j+i,1) = data_tmp(j,4);
            row{j+i} = row_tmp{j};
        end
        clear *_tmp
    end
    
    % Export
    if exist('X','var')
        T = table(X,Y,Z,R,'RowNames',row);
        % File Name - Work if EOS files are saved as 'project_number_visit_*' 
            tmp1 = files{1};
            tmp2 = strfind(tmp1,'_');
        op = upper(handles.Operator.String); % Operator
        file_name = [tmp1(1:tmp2(3)-1)]; 
        % Is there an export file already?
        if isempty(op)
            tmp3 = dir([dir_path,files{1}(1:8),'*.mat']);
        else
            tmp3 = dir([dir_path,files{1}(1:8),'*',op,'*.mat']);
        end
        %
        if isempty(tmp3)
            idx = 1; 
        else
            for k=1:size(tmp3,1)
                tmp5 = strfind(tmp3(k).name,'_');
                tmp6 = strfind(tmp3(k).name,'.');
                tmp7(k) = str2num(tmp3(k).name(tmp5(end)+1:tmp6-1));
            end
            idx = max(tmp7) + 1; 
        end
        
        % Export
        if isempty(op)
            out_name = [file_name,'_mrkEOS_',num2str(idx),'.xls'];
        else
            out_name = [file_name,'_mrkEOS_',op,'_',num2str(idx),'.xls'];       
        end
        writetable(T,[dir_path,out_name],'WriteRowNames',true)
        msgbox('Data was succesfully exported')
    else
        msgbox('No data to export')
    end
end

%% LOAD (.mat, .csv, .xls)
% --- Executes on button press in LoadProcessedData.
function LoadProcessedData_Callback(hObject, eventdata, handles)
% hObject    handle to LoadProcessedData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Info Front

% Load Data
[file,dir_path] = uigetfile('*.*','Select Processed Data');
ext = file(end-2:end);

switch ext
    case 'csv'
        csv = importdata([dir_path,file],',');
        % a - Anat Point
        % Index Anatomical Points 
        idx_anat = find(strcmp(csv.textdata,'RHJC')) - 1;
        % Fill table
        if isempty(idx_anat) ~= 1
            set(handles.Table_AnatPt,'data',csv.data(idx_anat:end,:))
        end
        % b - Skin Markers
        if isempty(idx_anat) ~=1 % are there anatomical points?
            set(handles.Table_SkinMrk,'data',csv.data(1:idx_anat-1,1:3))
            set(handles.Table_SkinMrk,'RowName',csv.textdata(2:idx_anat,1))
        else
            set(handles.Table_SkinMrk,'data',cav.data(1:end,1:3))
            set(handles.Table_SkinMrk,'RowName',csv.textdata(2:end,1))
        end
        
    case 'mat'
        load([dir_path,file])
        % Fill Tables
        % a - Skin Markers
        f_sm = fieldnames(SkinMrk);
        for i = 1:size(f_sm,1)
            data_tmp(i,:) = SkinMrk.(f_sm{i});
            row_tmp{i,:} = f_sm{i};
        end
        set(handles.Table_SkinMrk,'data',data_tmp)
        set(handles.Table_SkinMrk,'RowName',row_tmp); clear *_tmp

        % b - Anat Point
        f_ap = fieldnames(AnatPt);
        for i = 1:size(f_ap,1)
            data_tmp(i,:) = AnatPt.(f_ap{i});
            row_tmp{i,:} = f_ap{i};
        end
        set(handles.Table_AnatPt,'data',data_tmp)
        set(handles.Table_AnatPt,'RowName',row_tmp)
        
    case 'xls'
        [NUM,TXT] = xlsread([dir_path,file]);
        % Fill Tables
        % a - Anat Point
        % Index Anatomical Points 
        idx_anat = find(strcmp(TXT,'R_HJC')) - 1;
        % Fill table
        if isempty(idx_anat) ~= 1 % are there anatomical points?
            set(handles.Table_AnatPt,'data',NUM(idx_anat:end,:))
        end
        % b - Skin Markers
        if isempty(idx_anat) ~=1 % are there anatomical points?
            set(handles.Table_SkinMrk,'data',NUM(1:idx_anat-1,1:3))
            set(handles.Table_SkinMrk,'RowName',TXT(2:idx_anat,1))
        else % no anatomical points
            set(handles.Table_SkinMrk,'data',NUM(1:end,1:3))
            set(handles.Table_SkinMrk,'RowName',TXT(2:end,1))
        end
end

% Plot Data
px2mm = Info.Sagit.PixelSpacing(1); % from EOS - check PixelSpacing
xdim = size(Front,2);    % px size of frontal view
axes(handles.axes3);hold on
% - Anatomical points
am = handles.Table_AnatPt.Data / px2mm;
ar = handles.Table_AnatPt.RowName;

% /!\ There is a +1 in the plots, not sure why but it works
if isempty(am) == 0 
    for i=1:size(am,1)
        side = ar{i,1}(1);
        med  = ar{i,1}(2);
        switch side
            case 'R'
                line = '-';
            case 'L'
                line = '--';
        end
        switch med
            case 'M'
                c = 'g';
            otherwise
                c = 'w';
        end
        trace_cercle(am(i,1)+1, am(i,2)+1,am(i,4)+1,c,line)
        trace_cercle(am(i,3)+1+xdim,am(i,2)+1,am(i,4)+1,c,line)
    end
end
% - Skin Markers
sm = handles.Table_SkinMrk.Data / px2mm;
if isempty(sm) == 0 
plot(sm(:,1)+1,sm(:,2)+1,'*r')
plot(sm(:,3)+1+xdim,sm(:,2)+1,'*r')
end

end

%% other nothing interesting here
% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
end

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
end

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
end

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Directory_Callback(hObject, eventdata, handles)
% hObject    handle to Directory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Directory as text
%        str2double(get(hObject,'String')) returns contents of Directory as a double
end
% --- Executes during object creation, after setting all properties.
function Directory_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Directory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function FrontFile_Callback(hObject, eventdata, handles)
% hObject    handle to FrontFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrontFile as text
%        str2double(get(hObject,'String')) returns contents of FrontFile as a double
end

% --- Executes during object creation, after setting all properties.
function FrontFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrontFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function SagittalFile_Callback(hObject, eventdata, handles)
% hObject    handle to SagittalFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SagittalFile as text
%        str2double(get(hObject,'String')) returns contents of SagittalFile as a double
end

% --- Executes during object creation, after setting all properties.
function SagittalFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SagittalFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function r_cut_Callback(hObject, eventdata, handles)
% hObject    handle to r_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r_cut as text
%        str2double(get(hObject,'String')) returns contents of r_cut as a double
end

% --- Executes during object creation, after setting all properties.
function r_cut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function Operator_Callback(hObject, eventdata, handles)
% hObject    handle to Operator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Operator as text
%        str2double(get(hObject,'String')) returns contents of Operator as a double
end

% --- Executes during object creation, after setting all properties.
function Operator_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Operator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
