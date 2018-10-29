function varargout = RightMPI_b(varargin)
% RIGHTMPI_B MATLAB code for RightMPI_b.fig
%      RIGHTMPI_B, by itself, creates a new RIGHTMPI_B or raises the existing
%      singleton*.
%
%      H = RIGHTMPI_B returns the handle to a new RIGHTMPI_B or the handle to
%      the existing singleton*.
%
%      RIGHTMPI_B('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RIGHTMPI_B.M with the given input arguments.
%
%      RIGHTMPI_B('Property','Value',...) creates a new RIGHTMPI_B or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RightMPI_b_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RightMPI_b_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RightMPI_b

% Last Modified by GUIDE v2.5 26-Jul-2016 11:07:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RightMPI_b_OpeningFcn, ...
                   'gui_OutputFcn',  @RightMPI_b_OutputFcn, ...
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


% --- Executes just before RightMPI_b is made visible.
function RightMPI_b_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RightMPI_b (see VARARGIN)

% Choose default command line output for RightMPI_b
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

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes RightMPI_b wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RightMPI_b_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadImage.
function LoadImage_Callback(hObject, eventdata, handles)
% hObject    handle to LoadImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles=guidata(gca)
% hMainGui = getappdata(0, 'hMainGui');
clc
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
    handles.SignalFilePath= SignalFilePath;
    
else
    
    choice = questdlg('Want to choose a different signal file?','Loading Video','Yes','No, use current one','No, use current one');
    switch choice
        case 'Yes'
            axes(handles.axes1);cla;
            axes(handles.axes2);cla;
%             axes(handles.axes3);cla;
            %             set(handles.axes3, 'Visible', 'off');
            %             set(handles.table,'Data',[]);
            
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
%             set(handles.SignalFileNameText,'String',FileName);
            handles.SignalFileName=FileName;
        case 'No, use current one'
            SignalFilePath=handles.SignalFileNameText;
    end
end

handles.SignalFileNameText=SignalFilePath;



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
    dataNew = dataTemp;
     dataNew = zeros(size(dataTemp));
     dataNew(1:64,:)=dataTemp(65:128,:);
     dataNew(65:128,:)=dataTemp(1:64,:);
    imageSize = size(dataTemp);
    str = sprintf('[%d %d]',imageSize(1),imageSize(2));
    %     set(handles.ImageSize,'String',str);
    
    LoadedData = dataNew;
    %     handles.TimeStamp = time';
    
end

if isstruct(LoadedData)
    data=LoadedData.data;
else
    data=LoadedData;
end

handles.originalData = data;
handles.timeStamp = time;

axes(handles.axes1)
imshow(data,'Border','loose');

disp(FileName);

clc
guidata(hObject,handles);

% --- Executes on button press in rearrange.
function rearrange_Callback(hObject, eventdata, handles)
% hObject    handle to rearrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gca);
data = handles.originalData;
raw = data;

choice = questdlg('Rearragne the image?','Rearrange image','Flip','Shift','Shift then flip','Shift then flip');
switch choice
%     case 'Flip'
%         data1 = flipdim(raw(1:64,:),1);
%         data2 = flipdim(raw(65:128,:),1);
%         data(1:64,:)=data2;
%         data(65:128,:)=data1;
%         
%     case 'Shift'
%         data = zeros(size(raw));
%         data(1:64,:)=raw(65:128,:);
%         data(65:128,:)=raw(1:64,:);
%     case 'Shift then flip'
%         
%         data = zeros(size(raw));
%         data1 = flipdim(raw(1:64,:),1);
%         data2 = flipdim(raw(65:128,:),1);
%         data(1:64,:)=data1;
%         data(65:128,:)=data2;
% 
% 
    case 'Flip'
        data1 = flipdim(raw(1:64,:),1);
        data2 = flipdim(raw(65:128,:),1);
        data(1:64,:)=data2;
        data(65:128,:)=data1;
        
    case 'Shift'
        data = zeros(size(raw));
        data(1:80,:)=raw(49:128,:);
        data(81:128,:)=raw(1:48,:);
    case 'Shift then flip'
        
        data = zeros(size(raw));
        data1 = flipdim(raw(1:80,:),1);
        data2 = flipdim(raw(81:128,:),1);
        data(1:80,:)=data1;
        data(81:128,:)=data2;
end

axes(handles.axes1);cla;
axes(handles.axes1)
imshow(data,'Border','loose');
handles.originalData = data;

guidata(hObject,handles);


% --- Executes on button press in autoDetection.
function autoDetection_Callback(hObject, eventdata, handles)
% hObject    handle to autoDetection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gca);
data= handles.originalData;
time = handles.timeStamp;

[addevent,bInterval]=RightMPI_b_wave(data,time);
axes(handles.axes2)
imshow(addevent, 'Border','loose');

resultNum = num2cell(bInterval');
labelTit = {'Wave begin' 'Wave end' 'Interval' 'HR'};
result = [labelTit;resultNum];
assignin('base', 'result', result);

set(handles.table,'data',result); 
set(handles.table,'Visible','on'); 
set(handles.copy,'Visible','on'); 

guidata(hObject,handles);

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
% set(handles.copy,'Visible','off'); 

guidata(hObject,handles);

% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
% hObject    handle to exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(gca);

close(handles.output)


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
