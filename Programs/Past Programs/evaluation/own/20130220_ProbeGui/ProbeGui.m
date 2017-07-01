function varargout = ProbeGui(varargin)
% PROBEGUI M-file for ProbeGui.fig
%      PROBEGUI, by itself, creates a new PROBEGUI or raises the existing
%      singleton*.
%
%      H = PROBEGUI returns the handle to a new PROBEGUI or the handle to
%      the existing singleton*.
%
%      PROBEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROBEGUI.M with the given input arguments.
%
%      PROBEGUI('Property','Value',...) creates a new PROBEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProbeGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProbeGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProbeGui

% Last Modified by GUIDE v2.5 21-Feb-2013 14:54:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProbeGui_OpeningFcn, ...
                   'gui_OutputFcn',  @ProbeGui_OutputFcn, ...
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


% --- Executes just before ProbeGui is made visible.
function ProbeGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProbeGui (see VARARGIN)

% Choose default command line output for ProbeGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ProbeGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ProbeGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in PushbuttonOpen.
function PushbuttonOpen_Callback(hObject, eventdata, handles)
% hObject    handle to PushbuttonOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- File Dialog Box to open characteristic
FilterSpec = {'*.dat;*.raw;*.txt'};
[fn,path,FilterIndex] = uigetfile(FilterSpec);
handles.data.FileName = [path fn];
disp(['load ' handles.data.FileName])

A = load(handles.data.FileName);
plot(A(:,1), A(:,2))
set(gca,'xlim',[min(A(:,1)) max(A(:,1))])

guidata(hObject, handles)


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.data.ProbeTheory = get(get(handles.uipanel1, ...
  'SelectedObject'), 'String');

switch handles.data.ProbeTheory
 case 'Langmuir'
  set(handles.EditB,'Enable','off')
 case 'Demidov'
  set(handles.EditB,'Enable','on')
end




% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditB_Callback(hObject, eventdata, handles)
% hObject    handle to EditB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditB as text
%        str2double(get(hObject,'String')) returns contents of EditB as a double


% --- Executes during object creation, after setting all properties.
function EditB_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditProbeL_Callback(hObject, eventdata, handles)
% hObject    handle to EditProbeL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditProbeL as text
%        str2double(get(hObject,'String')) returns contents of EditProbeL as a double


% --- Executes during object creation, after setting all properties.
function EditProbeL_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditProbeL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EditProbeR_Callback(hObject, eventdata, handles)
% hObject    handle to EditProbeR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditProbeR as text
%        str2double(get(hObject,'String')) returns contents of EditProbeR as a double


% --- Executes during object creation, after setting all properties.
function EditProbeR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditProbeR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over PushbuttonOpen.
function PushbuttonOpen_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to PushbuttonOpen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in PushbuttonEval.
function PushbuttonEval_Callback(hObject, eventdata, handles)
% hObject    handle to PushbuttonEval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get string of selected radio button
handles.data.ProbeTheory = get(get(handles.uipanel1, ...
  'SelectedObject'), 'String');

ProbeGuiEvaluationA



function EditScanStepsize_Callback(hObject, eventdata, handles)
% hObject    handle to EditScanStepsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EditScanStepsize as text
%        str2double(get(hObject,'String')) returns contents of EditScanStepsize as a double


% --- Executes during object creation, after setting all properties.
function EditScanStepsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EditScanStepsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editVolAvg_Callback(hObject, eventdata, handles)
% hObject    handle to editVolAvg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editVolAvg as text
%        str2double(get(hObject,'String')) returns contents of editVolAvg as a double


% --- Executes during object creation, after setting all properties.
function editVolAvg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editVolAvg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenuFitGoal.
function popupmenuFitGoal_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuFitGoal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuFitGoal contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuFitGoal


% --- Executes during object creation, after setting all properties.
function popupmenuFitGoal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuFitGoal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
