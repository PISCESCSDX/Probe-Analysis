function varargout = Etcher_L_Probe_Gui_v2(varargin)
% ETCHER_L_PROBE_GUI_V2 MATLAB code for Etcher_L_Probe_Gui_v2.fig
%      ETCHER_L_PROBE_GUI_V2, by itself, creates a new ETCHER_L_PROBE_GUI_V2 or raises the existing
%      singleton*.
%
%      H = ETCHER_L_PROBE_GUI_V2 returns the handle to a new ETCHER_L_PROBE_GUI_V2 or the handle to
%      the existing singleton*.
%
%      ETCHER_L_PROBE_GUI_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ETCHER_L_PROBE_GUI_V2.M with the given input arguments.
%
%      ETCHER_L_PROBE_GUI_V2('Property','Value',...) creates a new ETCHER_L_PROBE_GUI_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Etcher_L_Probe_Gui_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Etcher_L_Probe_Gui_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Etcher_L_Probe_Gui_v2

% Last Modified by GUIDE v2.5 15-Jun-2017 16:33:43

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Etcher_L_Probe_Gui_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @Etcher_L_Probe_Gui_v2_OutputFcn, ...
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


% --- Executes just before Etcher_L_Probe_Gui_v2 is made visible.
function Etcher_L_Probe_Gui_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Etcher_L_Probe_Gui_v2 (see VARARGIN)

%Create tab group
handles.tgroup = uitabgroup('Parent', handles.figure1,'TabLocation', 'top');
handles.tab1 = uitab('Parent', handles.tgroup, 'Title', 'Main');
handles.tab2 = uitab('Parent', handles.tgroup, 'Title', 'Final Plot');

%Place panels into each tab
set(handles.P1,'Parent',handles.tab1)
set(handles.P2,'Parent',handles.tab2)


%Reposition each panel to same location as panel 1
set(handles.P2,'position',get(handles.P1,'position'));



%Set default Etcher Properties
set(handles.species_txt,'String','1');
set(handles.mass_txt,'String','2.014');
set(handles.rf_power_txt,'String','1300');
set(handles.magnet_txt,'String','50');
set(handles.flow_rate_txt,'String','12');
set(handles.pressure_txt,'String','7.0');
set(handles.bias_volt_txt,'String','-100.0');

%Set default probe properties
set(handles.res_txt,'String','50');
set(handles.tip_length_txt,'String','3.2e-3');
set(handles.tip_diameter_txt,'String','4.6e-4');

% Choose default command line output for Etcher_L_Probe_Gui_v2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Etcher_L_Probe_Gui_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Etcher_L_Probe_Gui_v2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function vin_path_txt_Callback(hObject, eventdata, handles)
if ~isempty(get(handles.vout_path_txt,'String'))
    set(handles.fit_data_button,'Enable','on');
end

% --- Executes during object creation, after setting all properties.
function vin_path_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vin_path_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function vout_path_txt_Callback(hObject, eventdata, handles)
if ~isempty(get(handles.vin_path_txt,'String'))
    set(handles.fit_data_button,'Enable','on');
end


% --- Executes during object creation, after setting all properties.
function vout_path_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vout_path_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in vin_button.
function vin_button_Callback(hObject, eventdata, handles)
% hObject    handle to vin_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename Pathname] = uigetfile('*.txt','Select Input File',...
    'C:\Users\Admin\Documents\MATLAB\probe');
handles.filename{1} = Filename;
handles.pathname{1} = Pathname;

if ischar(handles.filename{1}) && ischar(handles.pathname{1})
    set(handles.vin_path_txt,'String',[Pathname Filename]);
end

if ~isempty(get(handles.vout_path_txt,'String'))
    set(handles.fit_data_button,'Enable','on');
end
guidata(hObject,handles);


% --- Executes on button press in vout_button.
function vout_button_Callback(hObject, eventdata, handles)
% hObject    handle to vout_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename Pathname] = uigetfile('*.txt','Select Input File',...
    'C:\Users\Admin\Documents\MATLAB\probe');
handles.filename{2} = Filename;
handles.pathname{2} = Pathname;

if ischar(handles.filename{2}) && ischar(handles.pathname{2})
    set(handles.vout_path_txt,'String',[Pathname Filename]);
end

if ~isempty(get(handles.vin_path_txt,'String'))
    set(handles.fit_data_button,'Enable','on');
end
guidata(hObject,handles);

function save_path_txt_Callback(hObject, eventdata, handles)
% hObject    handle to save_path_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of save_path_txt as text
%        str2double(get(hObject,'String')) returns contents of save_path_txt as a double


% --- Executes during object creation, after setting all properties.
function save_path_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to save_path_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in save_dir_button.
function save_dir_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_dir_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Filename Pathname] = uiputfile('*.png','Save As',...
    'C:\Users\Admin\Documents\MATLAB\probe');

if ischar(Filename) && ischar(Pathname)
    set(handles.save_path_txt,'String',[Pathname Filename]);
end

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.tgroup.SelectedTab = handles.tab2;

fig = gcf;
fig.PaperPositionMode = 'auto';
print(get(handles.save_path_txt,'String'),'-dpng','-r0')

handles.tgroup.SelectedTab = handles.tab1;



function species_txt_Callback(hObject, eventdata, handles)
% hObject    handle to species_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of species_txt as text
%        str2double(get(hObject,'String')) returns contents of species_txt as a double


% --- Executes during object creation, after setting all properties.
function species_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to species_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mass_txt_Callback(hObject, eventdata, handles)
% hObject    handle to mass_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mass_txt as text
%        str2double(get(hObject,'String')) returns contents of mass_txt as a double


% --- Executes during object creation, after setting all properties.
function mass_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mass_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in H_button.
function H_button_Callback(hObject, eventdata, handles)
% hObject    handle to H_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of H_button
if get(hObject,'Value')== get(hObject,'Max')
    set(handles.species_txt,'String','1');
    set(handles.mass_txt,'String','1.008');
end

% --- Executes on button press in D_button.
function D_button_Callback(hObject, eventdata, handles)
% hObject    handle to D_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of D_button
if get(hObject,'Value')== get(hObject,'Max')
    set(handles.species_txt,'String','1');
    set(handles.mass_txt,'String','2.014');
end

% --- Executes on button press in He_button.
function He_button_Callback(hObject, eventdata, handles)
% hObject    handle to He_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of He_button
if get(hObject,'Value')== get(hObject,'Max')
    set(handles.species_txt,'String','2');
    set(handles.mass_txt,'String','4.003');
end

% --- Executes on button press in Ar_button.
function Ar_button_Callback(hObject, eventdata, handles)
% hObject    handle to Ar_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Ar_button
if get(hObject,'Value')== get(hObject,'Max')
    set(handles.species_txt,'String','18');
    set(handles.mass_txt,'String','39.948');
end

function uibuttongroup1_SelectionChangedFcn(hObject,eventdata,handles)
% hObject    handle to the selected object in uibuttongroup1 
% eventdata  structure with the following fields
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag')
    case 'H_button'
        
    case 'D_button'
        
    case 'He_button'
        
    case 'Ar_button'
end



function flow_rate_txt_Callback(hObject, eventdata, handles)
% hObject    handle to flow_rate_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of flow_rate_txt as text
%        str2double(get(hObject,'String')) returns contents of flow_rate_txt as a double


% --- Executes during object creation, after setting all properties.
function flow_rate_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to flow_rate_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function magnet_txt_Callback(hObject, eventdata, handles)
% hObject    handle to magnet_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of magnet_txt as text
%        str2double(get(hObject,'String')) returns contents of magnet_txt as a double


% --- Executes during object creation, after setting all properties.
function magnet_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to magnet_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



function rf_power_txt_Callback(hObject, eventdata, handles)
% hObject    handle to rf_power_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rf_power_txt as text
%        str2double(get(hObject,'String')) returns contents of rf_power_txt as a double


% --- Executes during object creation, after setting all properties.
function rf_power_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rf_power_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function pressure_txt_Callback(hObject, eventdata, handles)
% hObject    handle to pressure_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pressure_txt as text
%        str2double(get(hObject,'String')) returns contents of pressure_txt as a double


% --- Executes during object creation, after setting all properties.
function pressure_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pressure_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bias_volt_txt_Callback(hObject, eventdata, handles)
% hObject    handle to bias_volt_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bias_volt_txt as text
%        str2double(get(hObject,'String')) returns contents of bias_volt_txt as a double


% --- Executes during object creation, after setting all properties.
function bias_volt_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bias_volt_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function res_txt_Callback(hObject, eventdata, handles)
% hObject    handle to res_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of res_txt as text
%        str2double(get(hObject,'String')) returns contents of res_txt as a double


% --- Executes during object creation, after setting all properties.
function res_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to res_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tip_length_txt_Callback(hObject, eventdata, handles)
% hObject    handle to tip_length_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tip_length_txt as text
%        str2double(get(hObject,'String')) returns contents of tip_length_txt as a double


% --- Executes during object creation, after setting all properties.
function tip_length_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tip_length_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tip_diameter_txt_Callback(hObject, eventdata, handles)
% hObject    handle to tip_diameter_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tip_diameter_txt as text
%        str2double(get(hObject,'String')) returns contents of tip_diameter_txt as a double


% --- Executes during object creation, after setting all properties.
function tip_diameter_txt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tip_diameter_txt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit27_Callback(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit27 as text
%        str2double(get(hObject,'String')) returns contents of edit27 as a double


% --- Executes during object creation, after setting all properties.
function edit27_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fit_data_button.
function fit_data_button_Callback(hObject, eventdata, handles)
% hObject    handle to fit_data_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% ------------------------------------------------------
% define fitting defaults as percentages
nCut   = 3;             % Always cut off first and last n points
  
minI_s = 0.05;          % Start 05% between min and floating pot
maxI_s = 0.25;          % End   25%   " "
minN_i = 0.05;          % Start 05%
maxN_i = 0.25;          % End   25%
minT_e = 0.05;          % Start 05%
maxT_e = 0.17;          % End   17%
minE_s = 0.65;          % Start 65% 
maxE_s = 0.80;          % End   80%

%% ------------------------------------------------------
% define physical constants needed
q_e = 1.60217657e-19;                 % electron charge [C]
m_e = 9.10938291e-31;                 % electron mass   [kg]
m_n = 1.660538921e-27;                % nucleon mass    [kg]

%% ------------------------------------------------------
% plotting variables
promptFit.mark1 =  5.0;
promptFit.mark2 =  8.5;
promptFit.line1 =  1.5;

backColor = [0.95 0.95 0.95];
legColor  = [0.85 0.85 0.85];

%% ------------------------------------------------------
% get the experimental parameters from the user
Z        = str2double(get(handles.species_txt,'String'));
amu      = str2double(get(handles.mass_txt,'String'));
rfpower  = str2double(get(handles.rf_power_txt,'String'));
magnet   = str2double(get(handles.magnet_txt,'String'));
flow     = str2double(get(handles.flow_rate_txt,'String'));
pressure = str2double(get(handles.pressure_txt,'String'));
bias     = str2double(get(handles.bias_volt_txt,'String'));
auto     = get(handles.checkbox2,'Value');
Res      = str2double(get(handles.res_txt,'String'));
l_t      = str2double(get(handles.tip_length_txt,'String'));
d_t      = str2double(get(handles.tip_diameter_txt,'String'));

%fprintf('%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n%d\n',Z,amu,rfpower,magnet,...
%    flow,pressure,auto,Res,l_t,d_t)

%% ------------------------------------------------------
% get the files for analysis and put them in arrays
%load data from textboxes
file_vin  = get(handles.vin_path_txt,'String');
file_vout = get(handles.vout_path_txt,'String');

%load VP
fid = fopen(file_vin);
if fid == -1
    error(['Error occured when opening ' handles.filename{1}]);
end
V_p = fscanf(fid,'%g');
fclose(fid);

%load IP
fid = fopen(file_vout);
if fid == -1
    error(['Error occured when opening ' handles.filename{2}]);
end
I_p = fscanf(fid,'%g');         % Still in Volts, must convert to current
fclose(fid);

%placeholder for if I want to include this
%parsed_vin  = parseName(file_vin ); %naming purposes
%parsed_vout = parseName(file_vout);

%% ------------------------------------------------------
% Clean up data for plotting 
% cut off ends of the data
V_p =  V_p(nCut:end-nCut);
I_p =  I_p(nCut:end-nCut);

[~,iSort] = sort(V_p);
V_p = V_p(iSort);              %can't I just do this on the line before?
I_p = I_p(iSort);
I_p = -I_p/Res;                     % Current [A]

winSize  = round(length(I_p)*0.01);	% determine size of window for running average
winSize2 = round(winSize/2);

I_a = tsmovavg(I_p,'t',winSize,1);
I_a = circshift(I_a,-winSize2);

V_p = V_p(winSize2+1:end-winSize2-1);
I_p = I_p(winSize2+1:end-winSize2-1);
I_a = I_a(winSize2+1:end-winSize2-1);

% I_a = smooth(V_p,I_a,0.01,'rlowess'); % slows too much, skews 3rd fit

%% ---------------------------------------------------------
% plot the full I-V curve
% get min and max Vp to fit to

%{
if exist('savedFit','var')
    promptFit = savedFit;       
else
    promptFit.v0a = min(V_p); %sets the endpoints of the plot
    promptFit.v0b = max(V_p);
end
%}

promptFit.v0a = min(V_p); %sets the endpoints of the plot
promptFit.v0b = max(V_p);

promptFit.i0a = iFind(V_p,promptFit.v0a);
promptFit.i0b = iFind(V_p,promptFit.v0b);

axes(handles.axes1);
cla
set(gca,'XMinorTick','on')
grid on
plotIV_p(V_p,I_p,promptFit);

choice = auto;
while choice == 0
 vals = getRange(V_p,promptFit.i0a,promptFit.i0b,...
                 {'Enter the min V_p :',...
                  'Enter the max V_p :'});
 promptFit.v0a = vals(1);
 promptFit.v0b = vals(2);
 promptFit.i0a = iFind(V_p,promptFit.v0a);
 promptFit.i0b = iFind(V_p,promptFit.v0b);
 
 % replot the adjusted I-V curve
 
 cla
 set(gca,'XMinorTick','on')
 grid on
 plotIV_p(V_p,I_p,promptFit);
 
 choice = menu('Continue if plot is reasonable or change',...
               'change',...          % results in choice = 0
               'continue') - 1;      % results in choice = 1
end

V_p = V_p(promptFit.i0a:promptFit.i0b);        % save the newly cut vectors
I_p = I_p(promptFit.i0a:promptFit.i0b);
I_a = I_a(promptFit.i0a:promptFit.i0b);

%% ------------------------------------------------------
% limits for plotting
Vmin = min(V_p);
Vmax = max(V_p);
Imin = min(I_p);
Imax = max(I_p);
Vmm  = [Vmin,Vmax];
Imm  = [Imin,Imax];
promptFit.Vmin = Vmin;
promptFit.Vmax = Vmax;
promptFit.Imin = Imin;
promptFit.Imax = Imax;

%% ------------------------------------------------------
% fit points utilizing floating potential and above percentages
% find V_float and i_float
x0  = intersections(Vmm,[0,0],V_p,I_a);
V_f = mean(x0);
i_f = iFind(V_p,V_f);

promptFit.V_f = V_f;
promptFit.i_f = i_f;

% Initial guesses for fitting of data are either saved or default
%if exist('savedFit','var')
%else
    promptFit.v1a = Vmin + (V_f-Vmin)*minI_s;
    promptFit.v1b = Vmin + (V_f-Vmin)*maxI_s;
    promptFit.v2a = Vmin + (V_f-Vmin)*minN_i;
    promptFit.v2b = Vmin + (V_f-Vmin)*maxN_i;
    promptFit.v3a = V_f  + (Vmax-V_f)*minT_e;
    promptFit.v3b = V_f  + (Vmax-V_f)*maxT_e;
    promptFit.v4a = V_f  + (Vmax-V_f)*minE_s;
    promptFit.v4b = V_f  + (Vmax-V_f)*maxE_s;
%end

%% ------------------------------------------------------
% Initial guess for n_ion fitting
% promptFit.i2a = iFind(V_p,-90.0);
% promptFit.i2b = iFind(V_p,-75.0);

promptFit.i2a = iFind(V_p,promptFit.v2a);
promptFit.i2b = iFind(V_p,promptFit.v2b);

cla
set(gca,'XMinorTick','on')
grid on
dI2 = plotFitN_i(V_p,I_p,I_a,promptFit);

choice = auto;
while choice == 0 
 % get n_ion fitting parameters from user       
 vals = getRange(V_p,promptFit.i2a,promptFit.i2b,...
                 {'Enter the min V_p for n_ion fit:',...
                  'Enter the max V_p for n_ion fit:'});
 promptFit.v2a = vals(1);
 promptFit.v2b = vals(2);
 promptFit.i2a = iFind(V_p,promptFit.v2a);
 promptFit.i2b = iFind(V_p,promptFit.v2b);
 
 cla
 set(gca,'XMinorTick','on')
 grid on
 dI2 = plotFitN_i(V_p,I_p,I_a,promptFit);
 
 choice = menu('Continue if plot is reasonable or change',...
               'change',...          % results in choice = 0
               'continue') - 1;      % results in choice = 1
end


%% ------------------------------------------------------
% Initial guess for I_sat fitting
% promptFit.i1a = iFind(V_p,-90.0);
% promptFit.i1b = iFind(V_p,-75.0);

promptFit.i1a = iFind(V_p,promptFit.v1a);
promptFit.i1b = iFind(V_p,promptFit.v1b);

cla
set(gca,'XMinorTick','on')
grid on
I_sat = plotFitI_s(V_p,I_p,I_a,promptFit);

choice = auto;
while choice == 0 
 % get I_sat fitting parameters from user       
 vals = getRange(V_p,promptFit.i1a,promptFit.i1b,...
                 {'Enter the min V_p for I_sat fit:',...
                  'Enter the max V_p for I_sat fit:'});
 promptFit.v1a = vals(1);
 promptFit.v1b = vals(2);
 promptFit.i1a = iFind(V_p,promptFit.v1a);
 promptFit.i1b = iFind(V_p,promptFit.v1b);
 
 
 cla
 set(gca,'XMinorTick','on')
 grid on
 I_sat = plotFitI_s(V_p,I_p,I_a,promptFit);
 
 choice = menu('Continue if plot is reasonable or change',...
               'change',...          % results in choice = 0
               'continue') - 1;      % results in choice = 1
end


%% ------------------------------------------------------
% begin fit for T_e and E_sat
LI_p = -abs(log(I_p-I_sat));
LI_a = -abs(log(I_a-I_sat));

% initial guesses for Te and Esat fit
% promptFit.i3a = iFind(V_p,11.0);
% promptFit.i3b = iFind(V_p,14.0);
% promptFit.i4a = iFind(V_p,28.0);
% promptFit.i4b = iFind(V_p,35.0);

promptFit.i3a = iFind(V_p,promptFit.v3a);
promptFit.i3b = iFind(V_p,promptFit.v3b);
promptFit.i4a = iFind(V_p,promptFit.v4a);
promptFit.i4b = iFind(V_p,promptFit.v4b);

subplot(1,1,1);
handles.axes1 = findobj(gca,'type','axes');
cla
set(gca,'XMinorTick','on')
grid on
[V_s,T_e] = plotFitTE(V_p,LI_p,LI_a,promptFit);

choice = auto;
while choice == 0 
 % get T_e and E_sat fitting parameters from user       
 vals = getRange2(V_p,promptFit.i3a,promptFit.i3b,...
                      promptFit.i4a,promptFit.i4b,...
                  {'Enter the min V_p for T_e fit:',...
                   'Enter the max V_p for T_e fit:'...
                   'Enter the min V_p for E_sat fit',...
                   'Enter the max V_p for E_sat fit'});
 promptFit.v3a = vals(1);
 promptFit.v3b = vals(2);
 promptFit.i3a = iFind(V_p,promptFit.v3a);
 promptFit.i3b = iFind(V_p,promptFit.v3b);
 promptFit.v4a = vals(3);
 promptFit.v4b = vals(4);
 promptFit.i4a = iFind(V_p,promptFit.v4a);
 promptFit.i4b = iFind(V_p,promptFit.v4b);
 
 cla
 set(gca,'XMinorTick','on')
 grid on
 [V_s,T_e] = plotFitTE(V_p,LI_p,LI_a,promptFit);
 
 choice = menu('Continue if plot is reasonable or change',...
               'change',...          % results in choice = 0
               'continue') - 1;      % results in choice = 1
end

%% ------------------------------------------------------
% calculations and final plotting

% probe area
A_p = pi*(d_t/2)^2 + pi*d_t*l_t;  % [m^2]

% Ion flux and n_e
phi = I_sat/-q_e/A_p;             % [ions/m^2/s]
m_i = amu*m_n;                    % [kg] (Atomic mass)
n_e = phi/sqrt(q_e*T_e/m_i);      % [m^-3] (plasma density n_0) 
% exp(-Z/2)/sqrt(e*Z*T_e/m_i);    % [m^-3]

% Ion energy
E_ion = V_s-bias;                 % [eV]

% Ion density
n_ion = sqrt(pi^2*m_i/2/A_p^2/q_e^3*abs(dI2));  % from OML theory, m-3

% defining info for legends
IsatInfo = ['I_{sat} fit:'...
            '   V_{min} = ' num2str(promptFit.v1a,3)...
            '   V_{max} = ' num2str(promptFit.v1b,3)];
I2pInfo  = ['n_{ion} fit:'...
            '   V_{min} = ' num2str(promptFit.v2a,3)...
            '   V_{max} = ' num2str(promptFit.v2b,3)];                                          
TeInfo   = ['T_{e  } fit:'...
            '   V_{min} = ' num2str(promptFit.v3a,3)...
            '   V_{max} = ' num2str(promptFit.v3b,3)];                 
EsatInfo = ['E_{sat} fit:'...
            '   V_{min} = ' num2str(promptFit.v4a,3)...
            '   V_{max} = ' num2str(promptFit.v4b,3)];
TeEsatInfo = {TeInfo,EsatInfo};

%repeatFit = repeatPath;

%% ------------------------------------------------------
%------------------- Plotting -------------------------
%---------------- Final I-V curve ---------------------
%------------------------------------------------------

axes(handles.axes2);
title(IsatInfo)
plotI_s(V_p,I_p,I_a,promptFit);
h=legend('Raw Data',...
         'Run Avg',...
         'I_{sat} fit');
set(h,'Location','northwest','color',legColor)

axes(handles.axes3);
title(I2pInfo)
plotFitN_i(V_p,I_p,I_a,promptFit);
h=legend('Raw Data',...
         'Run Avg',...
         'n_{ion} fit');
set(h,'Location','southwest','color',legColor)

axes(handles.axes4);
title(TeEsatInfo)
plotFitTE(V_p,LI_p,LI_a,promptFit);
h = legend('Raw Data',...
           'Run Avg',...
           'T_e fit',...
           'E_{sat} fit',...
           'V_{plasma}');
set(h,'Location','southeast','color',legColor)

set(handles.axes5,'visible','off');
axes(handles.axes5);
% summary of measurements

phiInfo = '\phi_{ion} [m^{-2} s^{-1}] = ';

hrz1 = 0.2;
hrz2 = 0.6;
hght = 0.98;
delh = 0.1; 

text(hrz1,hght-delh*0,     'species [Z] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*1, 'atomic mass [u] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*2,    'rf power [W] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*3, 'main magnet [A] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*4, 'gas flow [sccm] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*5,   'pressure [mT] = ','HorizontalAlignment','right');
text(hrz1,hght-delh*6,        'bias [V] = ','HorizontalAlignment','right');

text(hrz2,hght-delh*0,   'V_{float} [V] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*1,  'V_{plasma} [V] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*2,        'T_e [eV] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*3,    'n_0 [m^{-3}] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*4,'n_{ion} [m^{-3}] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*5,    'E_{ion} [eV] = ','HorizontalAlignment','right');
text(hrz2,hght-delh*6,             phiInfo ,'HorizontalAlignment','right');

text(hrz1,hght-delh*0, num2str(       Z));
text(hrz1,hght-delh*1, num2str(     amu));
text(hrz1,hght-delh*2, num2str( rfpower));
text(hrz1,hght-delh*3, num2str(  magnet));
text(hrz1,hght-delh*4, num2str(    flow));
text(hrz1,hght-delh*5, num2str(pressure));
text(hrz1,hght-delh*6, num2str(    bias));

text(hrz2,hght-delh*0, num2str(   V_f,4));
text(hrz2,hght-delh*1, num2str(   V_s,4));
text(hrz2,hght-delh*2, num2str(   T_e,4));
text(hrz2,hght-delh*3, num2str(   n_e,4));
text(hrz2,hght-delh*4, num2str( n_ion,4));
text(hrz2,hght-delh*5, num2str( E_ion,4));
text(hrz2,hght-delh*6, num2str(   phi,4));

parsed_vin  = parseName(strcat(handles.pathname{1},handles.filename{1}));
parsed_vout  = parseName(strcat(handles.pathname{2},handles.filename{2}));

text(.10,.120,  'V_p path: ','HorizontalAlignment','right');
text(.10,.050, ' I_p path: ','HorizontalAlignment','right');
text(.10,.125,  parsed_vin ,'interpreter','none');
text(.10,.055,  parsed_vout,'interpreter','none');

handles.tgroup.SelectedTab = handles.tab2; %sets the active tab to the final plot tab
%disp('button press done')

function [Xvalue] = getRange( VP,iMin,iMax,pStr)
 prompt = pStr;
 name = 'SELECT DATA RANGE TO FIT';
 numlines = 1;
 minVi = num2str(VP(iMin),3);
 maxVi = num2str(VP(iMax),3);
 defAnswer = {minVi,maxVi};
 Xvalue = inputdlg(prompt,name,numlines,defAnswer);
 Xvalue = str2double(Xvalue);
 
function [Xvalue] = getRange2(VP,iMin,iMax,jMin,jMax,pStr)
 prompt = pStr;
 name = 'SELECT DATA RANGE TO FIT';
 numlines = 1;
 minVi = num2str(VP(iMin),3);
 maxVi = num2str(VP(iMax),3);
 minVj = num2str(VP(jMin),3);
 maxVj = num2str(VP(jMax),3);
 defAnswer = {minVi,maxVi,minVj,maxVj};
 Xvalue = inputdlg(prompt,name,numlines,defAnswer);
 Xvalue = str2double(Xvalue);
     
function [iVal] = iFind(Vi,xVal)
 iVal = 1;
 for i=1:size(Vi,1)
     if Vi(i)>xVal
        iVal = i;
        break
     else 
	iVal = length(Vi)-1;
     end
 end
 
 %% shorten the file name for display
function [fileName] = parseName(fileName)
 for i=1:20
    if i==1
        [t,rem] = strtok(fileName,'\');
        fileName = strcat(t,'\...');
    elseif i>2
        [t,rem] = strtok(rem,'\');
        if isempty(t)==1
            break
        else
            fileName = strcat(fileName,'\',t);
        end
    else
        [t,rem] = strtok(rem,'\');
    end
 end
 
function [] = plotIV_p(VP,IP,fitV)
  iMin = fitV.i0a;
  iMax = fitV.i0b;
  %figure(1)
  %clf
  title('I-V Curve')
  xlabel('V_p [V]')
  ylabel('I_p [mA]')
  axis([min(VP(iMin:iMax))     max(VP(iMin:iMax))...
        min(IP(iMin:iMax))*1e3 max(IP(iMin:iMax))*1e3])
  grid on
  hold on
  plot( VP(iMin:iMax),IP(iMin:iMax)*1e3,'db','MarkerSize', fitV.mark1)
  
% plot and fit n_ion curve
function [derI2,n_std] = plotFitN_i(VP,IP,IA,fitV)
 iMin = fitV.i2a;
 iMax = fitV.i2b;
 
 xMin = fitV.Vmin;
 xMax = fitV.V_f;
 yMin = 0;
 yMax = mean(IP(1:200))^2*1e6;
 
 p   = polyfit(VP(iMin:iMax),IA(iMin:iMax).^2,1);
 nF  = VP*p(1)+p(2);
 derI2 =   p(1);
 diffV = IP(iMin:iMax).^2 - nF(iMin:iMax).^2;
 n_std = std(diffV);
 
 V1 = [VP(iMin),VP(iMax)];
 I1 = [nF(iMin),nF(iMax)];
 
 xlabel('V_p [V]')
 ylabel('I^2_p [mA]^2')
 axis([xMin xMax yMin yMax])
 hold on
 plot( VP, (IP.^2)*1e6,  'xb',...
       VP, (IA.^2)*1e6,  '.c',...
       VP,      nF*1e6,  '-r',...
          'MarkerSize', fitV.mark1,...
           'LineWidth', fitV.line1);
 plot( V1,      I1*1e6,  'or',...
     'MarkerFaceColor',   'r',...
          'MarkerSize', fitV.mark2);

% plot and fit I_sat curve
function [I_s,I_std] = plotFitI_s(VP,IP,IA,fitV)
iMin = fitV.i1a;
iMax = fitV.i1b;

xMin = fitV.Vmin;
xMax = fitV.V_f;
yMin = fitV.Imin*1e3;
yMax = 0;

p   = polyfit(VP(iMin:iMax),IA(iMin:iMax),1);
IF  =  VP*p(1)+p(2);
IPS = -VP*p(1)+IP;
IAS = -VP*p(1)+IA;
I_s =  p(2);
Iss = [I_s,I_s];
VV  = [fitV.Vmin,fitV.V_f];
V1 = [VP(iMin),VP(iMax)];
I1 = [IF(iMin),IF(iMax)];

diffV = IP(iMin:iMax) - IF(iMin:iMax);
I_std = std(diffV);


ax1=subplot(2,1,1);
cla(ax1,'reset');
set(gca,'XMinorTick','on')
xlabel('V_p [V]')
ylabel('I_p [mA]')
axis([xMin xMax yMin yMax])
grid on
hold on
plot( VP, IP*1e3,  'xb',...
      VP, IA*1e3,  '.c',...
      VP, IF*1e3,  '-r',...
    'MarkerSize', fitV.mark1,...
     'LineWidth', fitV.line1)     
plot( V1, I1*1e3,  'or',...         % Show fitting limits
     'MarkerFaceColor', 'r',...
     'MarkerSize',fitV.mark2)
 
yMin = (I_s - (I_s-min(IA))*0.20)*1e3;
yMax = (I_s + (I_s-min(IA))*0.20)*1e3;

ax2=subplot(2,1,2);
cla(ax2,'reset');
set(gca,'XMinorTick','on')
xlabel('V_p [V]')
ylabel('I_p - I_{slope of sat} [mA]')
axis([xMin xMax yMin yMax])
grid on
hold on
plot( VP, IPS*1e3, 'xb',...
      VP, IAS*1e3, '.c',...
      VV, Iss*1e3,'--r',...
     'MarkerSize', fitV.mark1,...
      'LineWidth', fitV.line1)
  
 %% plot only the I_sat curve
function [I_s,I_std] = plotI_s(VP,IP,IA,fitV)
iMin = fitV.i1a;
iMax = fitV.i1b;

xMin = fitV.Vmin;
xMax = fitV.V_f;
yMin = fitV.Imin*1e3;
yMax = 0;

p   = polyfit(VP(iMin:iMax),IA(iMin:iMax),1);
IF  =  VP*p(1)+p(2);
I_s =  p(2);
V1 = [VP(iMin),VP(iMax)];
I1 = [IF(iMin),IF(iMax)];

diffV = IP(iMin:iMax) - IF(iMin:iMax);
I_std = std(diffV);

xlabel('V_p [V]')
ylabel('I_p [mA]')
axis([xMin xMax yMin yMax])
hold on
plot( VP, IP*1e3,  'xb',...
      VP, IA*1e3,  '.c',...
      VP, IF*1e3,  '-r',...
    'MarkerSize', fitV.mark1,...
     'LineWidth', fitV.line1)     
plot( V1, I1*1e3,  'or',...         % Show fitting limits
     'MarkerFaceColor', 'r',...
     'MarkerSize',    fitV.mark2)


function [V_s,T_e,T_std] = plotFitTE(VP,IP,IA,fitV)
iMin = fitV.i3a;
iMax = fitV.i3b;
jMin = fitV.i4a;
jMax = fitV.i4b;

xMin = fitV.V_f;
xMax = fitV.Vmax;
yMin = IA(fitV.i_f);
yMax = IA(end);

p1 = polyfit(VP(iMin:iMax),IA(iMin:iMax),1);
TF = VP*p1(1)+ p1(2);
T_e   = 1/p1(1);
diffV = IP(iMin:iMax) - TF(iMin:iMax);
T_std = std(diffV);

Vi = [VP(iMin),VP(iMax)];
Ii = [TF(iMin),TF(iMax)];

p2 = polyfit(VP(jMin:jMax),IA(jMin:jMax),1);
EF = VP*p2(1) + p2(2);
diffV = IP(jMin:jMax) - EF(jMin:jMax);
E_std = std(diffV);

Vj = [VP(jMin),VP(jMax)];
Ij = [EF(jMin),EF(jMax)];

% find V_plasma at intersection of T_e and n_ion
V_s = -(p1(2)-p2(2))/(p1(1)-p2(1));
Vss = [V_s,V_s];
Imm = [min(IA),max(IA)];

xlabel('V_p [V]')
ylabel('ln(I_p-I_{sat}) (A)')
axis([xMin xMax yMin yMax])
hold on
plot(     VP,  IP,  'xb',...
          VP,  IA,  '.c',...
          VP,  TF,  '-r',...               % fit with the initial guess
          VP,  EF,  '-g',...
         Vss, Imm, '--k',...
     'MarkerSize', fitV.mark1,...
      'LineWidth', fitV.line1)
plot(     Vi,  Ii,       'or',...
     'MarkerFaceColor',   'r',...
          'MarkerSize', fitV.mark2)
plot(     Vj,  Ij,       'og',...
     'MarkerFaceColor',   'g',...
          'MarkerSize', fitV.mark2)      
