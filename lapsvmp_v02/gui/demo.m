function varargout = demo(varargin)
% {demo} runs the manifold regularization GUI.
%
% Author: Stefano Melacci (2012)
%         mela@dii.unisi.it
%         * based on the code of Vikas Sindhwani, vikas.sindhwani@gmail.com         

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @demo_OpeningFcn, ...
                   'gui_OutputFcn',  @demo_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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

% --- Executes just before demo is made visible.
function demo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for demo
handles.output = hObject;

% UIWAIT makes demo wait for user response (see UIRESUME)
% uiwait(handles.figure1);
handles.datafile='';
handles.kernel='rbf';
handles.kernelparam=0.35;
handles.algorithm='LapSVM (Newton)';
handles.labeled=2;
handles.gamma_A=0.014362;
handles.gamma_I=0.7852;
handles.nn=6;
handles.graphdistance='euclidean';
handles.graphweights='heat';
handles.graphweightparam=0;
handles.laplaciandeg=1;
handles.laplaciannormalized=0;
handles.data.x=[];
handles.data.y=[];
handles.classfier=[];
handles.svplot=[];
handles.alphathres=0;
handles.discardbias=0;

datasets = dir('*.mat');
list = cell(1,1+length(datasets));
list{1} = 'SELECT ONE';
for i = 1:length(datasets)
    list{i+1} = datasets(i).name;
end

set(handles.dataset_popup,'String',list);
set(handles.distance_popup,'String',{'euclidean','cosine'});
set(handles.weights_popup,'String',{'heat','binary','distance'});
set(handles.kernel_popup,'String',{'rbf','linear','poly'});
set(handles.algorithm_popup,'String',{'LapSVM (Newton)','LapSVM (PCG)', ...
    'LapSVM (Dual)', 'LapRLSC','OCLapSVM (Newton)', 'OCLapSVM (PCG)', ...
    'SVM (Newton)','SVM (PCG)','OCSVM (Newton)', 'OCSVM (PCG)', ...
    'SVM (Dual)','RLSC'});
set(handles.gamma_A_slider,'Value',(log10(handles.gamma_A+10^(-5))+5)/7); 
set(handles.gamma_A_value,'String',num2str(handles.gamma_A));
set(handles.gamma_I_slider,'Value',(log10(handles.gamma_I+10^(-5))+5)/7);
set(handles.gamma_I_value,'String',num2str(handles.gamma_I));
set(handles.nn_value,'String',num2str(handles.nn));
set(handles.weights_value,'String',num2str(handles.graphweightparam));
set(handles.laplacian_degree_value,'String',num2str(handles.laplaciandeg));
set(handles.laplaciannorm_checkbox,'Value',handles.laplaciannormalized);
set(handles.labeled_value,'String',num2str(handles.labeled));
set(handles.kernelparam_value,'String',num2str(handles.kernelparam));
set(handles.alpha_discard_value,'String',num2str(handles.alphathres));

set(handles.generate_button,'Enable','off');
set(handles.run_button,'Enable','off');
set(handles.hack_button,'Enable','off');
set(handles.sv_checkbox,'Enable','off');

set(handles.bias_text,'String','');
set(handles.sv_text,'String','');
set(handles.minalpha_text,'String','');
set(handles.meanalpha_text,'String','');
set(handles.maxalpha_text,'String','');

guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = demo_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function gamma_A_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_A_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 0;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on slider movement.
function gamma_A_slider_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_A_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p=get(hObject,'Value');
handles.gamma_A=10^(7*p-5)-10^(-5);
set(handles.gamma_A_value,'String',num2str(handles.gamma_A));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function gamma_I_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_I_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on slider movement.
function gamma_I_slider_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_I_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
p=get(hObject,'Value');
handles.gamma_I=10^(7*p-5)-10^(-5);
set(handles.gamma_I_value,'String',num2str(handles.gamma_I));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function algorithm_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to algorithm_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in algorithm_popup.
function algorithm_popup_Callback(hObject, eventdata, handles)
% hObject    handle to algorithm_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns algorithm_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from algorithm_popup

contents = get(hObject,'String');
handles.algorithm=contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles
options=make_options();
options.NN=handles.nn;
options.Kernel=handles.kernel;
options.KernelParam=handles.kernelparam;
options.gamma_A=handles.gamma_A;
options.gamma_I=handles.gamma_I;
options.GraphDistanceFunction=handles.graphdistance;
options.GraphWeights=handles.graphweights;
options.LaplacianNormalize=handles.laplaciannormalized;
options.NewtonLineSearch=0;
options.GraphWeightParam=handles.graphweightparam;
options.LaplacianDegree=handles.laplaciandeg;
options.Verbose=0;
options.UseBias=1;

switch handles.algorithm
    
    case 'RLSC'
        data.X=handles.data.x(handles.data.Yss~=0,:);
        data.Y=handles.data.y(handles.data.Yss~=0,:);
        data.K=calckernel(options,data.X);
        handles.classifier=rlsc(options,data);

    case 'SVM (Newton)'
        data.X=handles.data.x(handles.data.Yss~=0,:);
        data.Y=handles.data.y(handles.data.Yss~=0,:);
        data.K=calckernel(options,data.X);
        options.Cg=0;   
        options.MaxIter=30;
        handles.classifier=svmp(options,data);
               
    case 'SVM (PCG)'
        data.X=handles.data.x(handles.data.Yss~=0,:);
        data.Y=handles.data.y(handles.data.Yss~=0,:);
        data.K=calckernel(options,data.X);
        options.Cg=1;        
        options.CgStopType=0;
        options.MaxIter=200;       
        handles.classifier=svmp(options,data);           
                 
    case 'OCSVM (Newton)'          
        data.X=handles.data.x(handles.data.Yss~=0,:);
        data.Y=handles.data.y(handles.data.Yss~=0,:);
        data.K=calckernel(options,data.X);
        options.Cg=0;   
        options.MaxIter=30;        
        handles.classifier=ocsvmp(options,data);       
       
    case 'OCSVM (PCG)'        
        data.X=handles.data.x(handles.data.Yss~=0,:);
        data.Y=handles.data.y(handles.data.Yss~=0,:);
        data.K=calckernel(options,data.X);
        options.Cg=1;        
        options.CgStopType=0;
        options.MaxIter=200;
        handles.classifier=ocsvmp(options,data);   
        
    case 'SVM (Dual)'
        data.X=handles.data.x(handles.data.Yss~=0,:);
        data.Y=handles.data.y(handles.data.Yss~=0,:);
        data.K=calckernel(options,data.X);        
        handles.classifier=svm(options,data);        

    case 'LapRLSC'
        data.X=handles.data.x;
        data.Y=handles.data.Yss;
        data.K=calckernel(options,data.X);
        data.L=laplacian(options,data.X);
        if options.LaplacianNormalize, options.UseBias=0; end
        handles.classifier=laprlsc(options,data);

    case 'LapSVM (Dual)'
        data.X=handles.data.x;
        data.Y=handles.data.Yss;
        data.K=calckernel(options,data.X);
        data.L=laplacian(options,data.X);
        if options.LaplacianNormalize, options.UseBias=0; end
        handles.classifier=lapsvm(options,data);            
        
    case 'LapSVM (Newton)'
        data.X=handles.data.x;
        data.Y=handles.data.Yss;
        data.K=calckernel(options,data.X);
        data.L=laplacian(options,data.X);           
        options.Cg=0;   
        options.MaxIter=30;
        if options.LaplacianNormalize, options.UseBias=0; end
        handles.classifier=lapsvmp(options,data);    
        
    case 'LapSVM (PCG)'     
        data.X=handles.data.x;
        data.Y=handles.data.Yss;
        data.K=calckernel(options,data.X);
        data.L=laplacian(options,data.X);     
        options.Cg=1;        
        options.CgStopType=0;
        options.MaxIter=200;
        if options.LaplacianNormalize, options.UseBias=0; end
        handles.classifier=lapsvmp(options,data); 
        
    case 'OCLapSVM (Newton)'    
        data.X=handles.data.x;
        data.Y=handles.data.Yss;
        data.K=calckernel(options,data.X);
        data.L=laplacian(options,data.X);     
        options.Cg=0;   
        options.MaxIter=30;
        if options.LaplacianNormalize, options.UseBias=0; end
        handles.classifier=oclapsvmp(options,data);         
        
    case 'OCLapSVM (PCG)'  
        data.X=handles.data.x;
        data.Y=handles.data.Yss;
        data.K=calckernel(options,data.X);
        data.L=laplacian(options,data.X);           
        options.Cg=1;        
        options.CgStopType=0;
        options.MaxIter=200;
        if options.LaplacianNormalize, options.UseBias=0; end
        handles.classifier=oclapsvmp(options,data);       
end

do_plotclassifier(handles);

set(handles.hack_button,'Enable','on');
set(handles.sv_checkbox,'Enable','on');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function dataset_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dataset_popup (see GCBO)z
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in dataset_popup.
function dataset_popup_Callback(hObject, eventdata, handles)
% hObject    handle to dataset_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns dataset_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from dataset_popup
contents = get(hObject,'String');
handles.datafile = contents{get(hObject,'Value')};

if isempty(strfind(handles.datafile,'.mat'))
    handles.datafile = '';
    return
end

handles.data = load(handles.datafile, 'x', 'y');
if isfield(handles.data,'x') ~= 0 && isfield(handles.data,'y') ~= 0
    set(handles.generate_button,'Enable','on');
    set(handles.run_button,'Enable','off');
    set(handles.hack_button,'Enable','off');
    set(handles.sv_checkbox,'Enable','off');
    set(handles.bias_text,'String','');
    set(handles.sv_text,'String','');
    handles.classifier=[];
    
    hold on;    
    h = gca; 
    cla(h);    
    set(h,'Xcolor','w'); 
    set(h,'Ycolor','w');
    plot2D(handles.data.x,handles.data.y);
    title([handles.datafile ' ' ...
           num2str(size(handles.data.x,1)) ' data points'],'Color','w');
else
    set(handles.generate_button,'Enable','off');
    set(handles.run_button,'Enable','off');
    set(handles.hack_button,'Enable','off');
    set(handles.sv_checkbox,'Enable','off');
    handles.classifier=[];
    handles.data.x = [];
    handles.data.y = [];
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function labeled_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to labeled_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function labeled_value_Callback(hObject, eventdata, handles)
% hObject    handle to labeled_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of labeled_value as text
%        str2double(get(hObject,'String')) returns contents of labeled_value as a double
handles.labeled=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function kernelparam_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernelparam_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function kernelparam_value_Callback(hObject, eventdata, handles)
% hObject    handle to kernelparam_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of kernelparam_value as text
%        str2double(get(hObject,'String')) returns contents of kernelparam_value as a double

handles.kernelparam=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes on button press in generate_button.
function generate_button_Callback(hObject, eventdata, handles)
% hObject    handle to generate_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.cgstep=-1;

if (handles.labeled>0)    
    ipos=find(handles.data.y>0);
    ineg=find(handles.data.y<0);

    rpos=randperm(length(ipos));
    rneg=randperm(length(ineg));

    if handles.labeled>length(ipos)
        handles.labeled=length(ipos);
    end
    if handles.labeled>length(ineg)
        handles.labeled=length(ineg);
    end

    handles.data.Yss=zeros(length(handles.data.y),1);
    handles.data.Yss(ipos(rpos(1:handles.labeled)))=1;
    handles.data.Yss(ineg(rneg(1:handles.labeled)))=-1;  
else
    handles.data.Yss=handles.data.y;
end

h = gca;
hold on;
cla(h);
set(h,'Xcolor','w');
set(h,'Ycolor','w');
plot2D(handles.data.x,handles.data.Yss);

set(handles.run_button,'Enable','on');
set(handles.hack_button,'Enable','off');
set(handles.sv_checkbox,'Enable','off');
set(handles.bias_text,'String','');
set(handles.sv_text,'String','');
set(handles.minalpha_text,'String','');
set(handles.meanalpha_text,'String','');
set(handles.maxalpha_text,'String','');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function gamma_A_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_A_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function gamma_A_value_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_A_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma_A_value as text
%        str2double(get(hObject,'String')) returns contents of gamma_A_value as a double

handles.gamma_A=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function gamma_I_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gamma_I_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function gamma_I_value_Callback(hObject, eventdata, handles)
% hObject    handle to gamma_I_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gamma_I_value as text
%        str2double(get(hObject,'String')) returns contents of gamma_I_value as a double

handles.gamma_I=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes on selection change in kernel_popup.
function kernel_popup_Callback(hObject, eventdata, handles)
% hObject    handle to kernel_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns kernel_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from kernel_popup
contents = get(hObject,'String');
handles.kernel=contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function kernel_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kernel_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nn_value_Callback(hObject, eventdata, handles)
% hObject    handle to nn_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nn_value as text
%        str2double(get(hObject,'String')) returns contents of nn_value as a double

handles.nn=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function nn_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nn_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function laplacian_degree_value_Callback(hObject, eventdata, handles)
% hObject    handle to laplacian_degree_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of laplacian_degree_value as text
%        str2double(get(hObject,'String')) returns contents of laplacian_degree_value as a double
handles.laplaciandeg=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function laplacian_degree_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to laplacian_degree_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function weights_value_Callback(hObject, eventdata, handles)
% hObject    handle to weights_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of weights_value as text
%        str2double(get(hObject,'String')) returns contents of weights_value as a double
handles.graphweightparam=str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function weights_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to weights_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in distance_popup.
function distance_popup_Callback(hObject, eventdata, handles)
% hObject    handle to distance_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns distance_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from distance_popup
contents = get(hObject,'String');
handles.graphdistance=contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function distance_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distance_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in weights_popup.
function weights_popup_Callback(hObject, eventdata, handles)
% hObject    handle to weights_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns weights_popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from weights_popup
contents = get(hObject,'String');
handles.graphweights=contents{get(hObject,'Value')};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function weights_popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to weights_popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in laplaciannorm_checkbox.
function laplaciannorm_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to laplaciannorm_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of laplaciannorm_checkbox
handles.laplaciannormalized=get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in bias_checkbox.
function bias_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to bias_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bias_checkbox
handles.discardbias = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in hack_button.
function hack_button_Callback(hObject, eventdata, handles)
% hObject    handle to hack_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.discardbias
    b_backup=handles.classifier.b;
    handles.classifier.b=0;
end
if handles.alphathres>0
    valid=abs(handles.classifier.alpha)>handles.alphathres;
    handles.classifier.svs=handles.classifier.svs(valid);
    handles.classifier.alpha=handles.classifier.alpha(valid);
    handles.classifier.xtrain=handles.classifier.xtrain(valid,:);
end

do_plotclassifier(handles);

if handles.discardbias
    handles.classifier.b=b_backup;
end
guidata(hObject,handles);

function alpha_discard_value_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_discard_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_discard_value as text
%        str2double(get(hObject,'String')) returns contents of alpha_discard_value as a double
handles.alphathres = str2double(get(hObject,'String'));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function alpha_discard_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_discard_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in sv_checkbox.
function sv_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to sv_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sv_checkbox
if get(hObject,'Value')
    handles.svplot=plot(gca,handles.classifier.xtrain(:,1), ...
                        handles.classifier.xtrain(:,2), ...
                        'mo','MarkerSize',12);
else
    delete(handles.svplot);
end
guidata(hObject,handles);


function do_plotclassifier(handles)

rmin=min(min(handles.data.x))-0.2;
rmax=max(max(handles.data.x))+0.2;
steps=(rmax-rmin)/100;

% ------------------------------------------
% start of tmp section
% ------------------------------------------
% screen = get(0, 'ScreenSize');
% fig = figure('Position',[(screen(3)/2) (screen(4)/2) 280 180], ...
%              'NumberTitle', 'on');
% axesfig = gca;
% set(axesfig,'Box', 'on', 'XMinorTick', 'off', 'YMinorTick', 'off', ...
%     'Parent', fig, 'FontSize', 10);
% steps=(rmax-rmin)/250;
% ------------------------------------------
% end of tmp section
% ------------------------------------------

xrange=rmin:steps:rmax+0.1;
yrange=rmin:steps:rmax+0.1;

hold on;
h = gca;
cla(h);
set(h,'Xcolor','w');
set(h,'Ycolor','w');

plotclassifier(handles.classifier,xrange,yrange);    

title([handles.datafile ' (' handles.algorithm ') ' ...
       num2str(size(handles.data.x,1)) ' data points'],'Color','w');
plot2D(handles.data.x,handles.data.Yss);

abs_alpha=abs(handles.classifier.alpha);
set(handles.bias_text,'String',num2str(handles.classifier.b));
set(handles.sv_text,'String',num2str(length(handles.classifier.svs)));
set(handles.minalpha_text,'String',num2str(min(abs_alpha)));
set(handles.meanalpha_text,'String',num2str(mean(abs_alpha)));
set(handles.maxalpha_text,'String',num2str(max(abs_alpha)));
set(handles.sv_checkbox,'Value',0);

% ------------------------------------------
% start of tmp section
% ------------------------------------------
% xlabel('');title('');
% ylim([-1,1.5]);
% eval(['print -f' num2str(gcf) ' -depsc test']);
% saveas(gcf,'test.fig');
% ------------------------------------------
% end of tmp section
% ------------------------------------------
