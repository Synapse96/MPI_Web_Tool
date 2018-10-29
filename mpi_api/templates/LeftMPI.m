function varargout = LeftMPI(varargin)
% LEFTMPI MATLAB code for LeftMPI.fig
%      LEFTMPI, by itself, creates a new LEFTMPI or raises the existing
%      singleton*.
%
%      H = LEFTMPI returns the handle to a new LEFTMPI or the handle to
%      the existing singleton*.
%
%      LEFTMPI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEFTMPI.M with the given input arguments.
%
%      LEFTMPI('Property','Value',...) creates a new LEFTMPI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LeftMPI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LeftMPI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LeftMPI

% Last Modified by GUIDE v2.5 26-Jul-2016 15:27:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LeftMPI_OpeningFcn, ...
    'gui_OutputFcn',  @LeftMPI_OutputFcn, ...
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


% --- Executes just before LeftMPI is made visible.
function LeftMPI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LeftMPI (see VARARGIN)

% Choose default command line output for LeftMPI
handles.output = hObject;
setappdata(0,'h',gcf);


handles.DataLoaded=0;

handles.SignalFilePath=NaN;
currentfolder = pwd;
handles.DatabasePath=currentfolder(1:(end-4));

handles.halfID = nan;
handles.TimeStamp = nan;


handles.AO = nan;
handles.AC = nan;
handles.MO = nan;
handles.MC = nan;

axes(handles.axes1);

axes(handles.axes2);

% axes(handles.axes3);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LeftMPI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Executes when GUI is resized.
function GUI_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Outputs from this function are returned to the command line.
function varargout = LeftMPI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Load_image.
function Load_image_Callback(hObject, eventdata, handles)
% hObject    handle to Load_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Please not in this function there is a middle line set to 64, which can
% be changed depending on the image

handles=guidata(gca);
% hMainGui = getappdata(0, 'hMainGui');

SignalFilePath=handles.SignalFilePath;
DatabasePath=handles.DatabasePath;


if isnan(SignalFilePath)
    
    [FileName,PathName, FilterIndex]=uigetfile('*','Choose the signal file', DatabasePath);
    if isequal(FileName,0)
        msgbox('No signal has been selected');
        return;
    end
    SignalFilePath=[PathName FileName];
    set(handles.SignalFileNameText,'String',FileName);
    handles.SignalFileName=FileName;
    
    
else
    
    choice = questdlg('Want to choose a different signal file?','Loading Video','Yes','No, use current one','No, use current one');
    switch choice
        case 'Yes'
            axes(handles.axes1);cla;
            axes(handles.axes2);cla;
            
            handles.AO = nan;
            handles.AC = nan;
            handles.MO = nan;
            handles.MC = nan;
            
            
            [FileName,PathName, FilterIndex]=uigetfile('*','Choose the signal file', DatabasePath);
            if isequal(FileName,0)
                msgbox('No signal has been selected');
                return;
            end
            SignalFilePath=[PathName FileName];
            set(handles.SignalFileNameText,'String',FileName);
            handles.SignalFileName=FileName;
        case 'No, use current one'
            SignalFilePath=handles.SignalFilePath;
    end
end

handles.SignalFilePath=SignalFilePath;



fileformat = SignalFilePath((length(SignalFilePath)-3):end);

if (strcmp(fileformat,'.mha'));
    info = getInformation(SignalFilePath);
    imageSize = info.DimSize;
    str = sprintf('[%d %d]',imageSize(1),imageSize(2));
%     set(handles.ImageSize,'String',str);
    LoadedData = ReadMetaFile(info);
    
else
    
    dataTemp = h5read(SignalFilePath, '/MovieGroup2/AcqPWCW/RawData/RawDataUnit');
    time = h5read(SignalFilePath, '/MovieGroup2/AcqPWCW/RawData/TimeStamp');
    
    dataNew = zeros(size(dataTemp));
    dataNew(1:64,:)=dataTemp(65:128,:);
    dataNew(65:128,:)=dataTemp(1:64,:);
    imageSize = size(dataTemp);
    str = sprintf('[%d %d]',imageSize(1),imageSize(2));
    %     set(handles.ImageSize,'String',str);
    
    LoadedData = dataNew;
    handles.TimeStamp = time;
    
end

if isstruct(LoadedData)
    data=LoadedData.data;
else
    data=LoadedData;
end

handles.originalData = data;
axes(handles.axes1)
imshow(data,'Border','loose');

disp(FileName);

clc
guidata(hObject,handles);



% --- Executes on button press in Re_arrange.
function Re_arrange_Callback(hObject, eventdata, handles)
% hObject    handle to Re_arrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%
% Please note in this function there is a baseline (middle line) for
% image rearrange set to 64, which can be changed depending on the image;
% Better implementation can be allowing user to select baseline
%%

handles = guidata(gca);
data = handles.originalData;
raw = data;

choice = questdlg('Rearragne the image?','Rearrange image','Flip','Shift','Shift then flip','Shift then flip');
switch choice
    case 'Flip'
        data1 = flipdim(raw(1:64,:),1);
        data2 = flipdim(raw(65:128,:),1);
        data(1:64,:)=data2;
        data(65:128,:)=data1;
        
    case 'Shift'
        data = zeros(size(raw));
        data(1:64,:)=raw(65:128,:);
        data(65:128,:)=raw(1:64,:);
    case 'Shift then flip'
        
        data = zeros(size(raw));
        data1 = flipdim(raw(1:64,:),1);
        data2 = flipdim(raw(65:128,:),1);
        data(1:64,:)=data1;
        data(65:128,:)=data2;
        
end
axes(handles.axes1);cla;
axes(handles.axes1)
imshow(data,'Border','loose');
handles.originalData = data;

guidata(hObject,handles);


% --- Executes on button press in Valve_detection.
function Valve_detection_Callback(hObject, eventdata, handles)
% hObject    handle to Valve_detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gca);
data= handles.originalData;
time = handles.TimeStamp;
[result,resultNumR,MC,valMC,AO,valAO,AC,valAC,MO,valMO]=mainfile1028F(data,time);

ultra = replicateUltrasound(data);
ultra = imadjust(ultra,[0.05 0.7],[0 0.69],1.2);

plotLine = 1:2:floor(size(ultra,1)/2)*2;
if ~isempty([AO MO])
    mask = ultra<0;
    mask(plotLine,[AO MO]) = 1;
    addOP = imoverlay(ultra,mask,[0 1 0]);
else addOP = ultra;
    
end
if ~isempty([MC AC])
    mask = ultra<0;
    mask(plotLine,[MC AC]) = 1;
    addED = imoverlay(addOP,mask,[1 1 0]);
else   addED = ultra ;
end
axes(handles.axes2)
imshow(addED, 'Border','loose');



set(handles.table,'data',result); 
set(handles.table,'Visible','on'); 
set(handles.copy,'Visible','on'); 

assignin('base', 'result', result)

% axes(handles.axes3)
% imshow(addED,'Border','loose');
guidata(hObject,handles);


% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gca);

close(handles.output)


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gca);
axes(handles.axes1);
cla

axes(handles.axes2);
cla


guidata(hObject,handles);


% --- Executes on button press in copy.
function copy_Callback(hObject, eventdata, handles)
% hObject    handle to copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
resultTable = get(handles.table,'data');
d = cell2mat(resultTable(2:end,:));
% d = rand(5,2);
size_d = size(d);
str = '';
for i=1:size_d(1)
    for j=1:size_d(2)
        str = sprintf ( '%s%f\t', str, d(i,j) );
    end
    str = sprintf ( '%s\n', str );
end
clipboard ('copy', str );
guidata(hObject,handles);
