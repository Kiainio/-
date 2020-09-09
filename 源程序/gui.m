function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 25-Jun-2020 15:25:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in origin.
function origin_Callback(hObject, eventdata, handles)
% hObject    handle to origin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    c = handles.Fs;  
    y = handles.y1;
    %[B, A] = butter(5, 700*2/c, 'low');
    %y = filter(B, A, y);
    sound(y, c);
    plot(handles.axes1, handles.y1)
    title(handles.axes1, '时域图');
    ysize = size(handles.y1);
    y1 = fft(handles.y, length(handles.y1));
    ysize = size(y1);
    plot(handles.axes2, abs(y1(1: ysize/2)))
    xlabel(handles.axes2, '频率');
    ylabel(handles.axes2, '幅度');
    title(handles.axes2, '频率特性');
    guidata(hObject, handles);


% --- Executes on button press in boy2girl.
function boy2girl_Callback(hObject, eventdata, handles)
% hObject    handle to boy2girl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    FL = 80;
    WL = 240;
    P = 10;
    x1 = handles.y;
    fs = 8000;
    x1 = x1 / max(x1);
    L = length(x1);
    FN = floor(L / FL) - 2;
    exc = zeros(L, 1);
    zi_pre = zeros(P, 1);
    xl_rec = zeros(L, 1);
    zi_rec = zeros(P, 1);
    exc_syn_t = zeros(L, 1);
    x1_syn_t = zeros(L,1);
    last_syn_t = 0;
    zi_syn_t = zeros(P,1);
    hw = hamming(WL);
    for n = 3 : FN
        x1_w = x1(n * FL - WL + 1 : n * FL).*hw;
        [A, E]=lpc(x1_w,P);
        x1_f = x1((n - 1) * FL + 1 : n * FL);
        [excl, zi_pre] = filter(A, 1, x1_f, zi_pre);
        exc((n-1) * FL + 1 : n * FL) = excl;
        [x1_rec1, zi_rec] = filter(1, A, excl, zi_rec);
        x1_rec((n - 1) * FL + 1 : n * FL) = x1_rec1;
        x1_Pitch = exc(n * FL - 222 : n * FL);
        PT = findpitch(x1_Pitch);
        G = sqrt(E * PT);
        PT1 = floor(PT / 1.5);
        poles = roots(A);
        delta = 150 * 2 * pi / 8000;
        for p = 1 : 10
            if imag(poles(p)) > 0
                poles(p) = poles(p) * exp(1i * delta);
            elseif imag(poles(p)) < 0
                poles(p) = poles(p) * exp(-1i * delta);
            end
        end
        A1 = poly(poles);
        tempn_syn_t = (1 : n * FL - last_syn_t);
        exc_synl_t = zeros(length(tempn_syn_t), 1);
        exc_synl_t(mod(tempn_syn_t, PT1)==0) = G;
        exc_synl_t = exc_synl_t((n-1) * FL - last_syn_t + 1 : n * FL - last_syn_t);
        [xl_synl_t, zi_syn_t] = filter(1, A1, exc_synl_t, zi_syn_t);
        exc_syn_t((n - 1) * FL + 1 : n * FL) = exc_synl_t;
        xl_syn_t((n - 1) * FL + 1 : n * FL) = xl_synl_t;
        last_syn_t = last_syn_t + PT1 * floor((n * FL - last_syn_t)/PT1);
    end
    [B, A] = butter(5, [200/4000 600/4000]);
    xl_syn_t_1 = filter(B, A, xl_syn_t);
    [B, A] = butter(5, [900/4000 1500/4000]);
    xl_syn_t_2 = filter(B, A, xl_syn_t);
    x1_syn_t = smoothdata(x1_syn_t, 'sgolay', WL);
    xl_syn_t = xl_syn_t_1 + xl_syn_t_2;
    sound(xl_syn_t, fs);
    plot(handles.axes1, xl_syn_t)
    title(handles.axes1, '时域图');
    grid on;
    ysize = size(x1_syn_t);
    y = fft(xl_syn_t, length(xl_syn_t));
    ysize = size(y);
    plot(handles.axes2, abs(y))
    xlabel(handles.axes2, '频率');
    ylabel(handles.axes2, '幅度');
    title(handles.axes2,'频率特征');
    guidata(hObject, handles);
    
   


% --- Executes on button press in girl2boy.
function girl2boy_Callback(hObject, eventdata, handles)
% hObject    handle to girl2boy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    FL = 80;
    WL = 240;
    P = 10;
    x1 = handles.y;
    fs = 8000;
    x1 = x1 / max(x1);
    L = length(x1);
    FN = floor(L / FL) - 2;
    exc = zeros(L, 1);
    zi_pre = zeros(P, 1);
    xl_rec = zeros(L, 1);
    zi_rec = zeros(P, 1);
    exc_syn_t = zeros(L, 1);
    x1_syn_t = zeros(L,1);
    last_syn_t = 0;
    zi_syn_t = zeros(P,1);
    hw = hamming(WL);
    for n = 3 : FN
        x1_w = x1(n * FL - WL + 1 : n * FL).*hw;
        [A, E]=lpc(x1_w,P);
        x1_f = x1((n - 1) * FL + 1 : n * FL);
        [excl, zi_pre] = filter(A, 1, x1_f, zi_pre);
        exc((n-1) * FL + 1 : n * FL) = excl;
        [x1_rec1, zi_rec] = filter(1, A, excl, zi_rec);
        x1_rec((n - 1) * FL + 1 : n * FL) = x1_rec1;
        x1_Pitch = exc(n * FL - 222 : n * FL);
        PT = findpitch(x1_Pitch);
        G = sqrt(E * PT);
        PT1 = floor(PT / 0.5);
        poles = roots(A);
        delta = 60 * 2 * pi / 8000;
        for p = 1 : 10
            if imag(poles(p)) > 0
                poles(p) = poles(p) * exp(1i * delta);
            elseif imag(poles(p)) < 0
                poles(p) = poles(p) * exp(-1i * delta);
            end
        end
        A1 = poly(poles);
        tempn_syn_t = (1 : n * FL - last_syn_t);
        exc_synl_t = zeros(length(tempn_syn_t), 1);
        exc_synl_t(mod(tempn_syn_t, PT1)==0) = G;
        exc_synl_t = exc_synl_t((n-1) * FL - last_syn_t + 1 : n * FL - last_syn_t);
        [xl_synl_t, zi_syn_t] = filter(1, A1, exc_synl_t, zi_syn_t);
        exc_syn_t((n - 1) * FL + 1 : n * FL) = exc_synl_t;
        xl_syn_t((n - 1) * FL + 1 : n * FL) = xl_synl_t;
        last_syn_t = last_syn_t + PT1 * floor((n * FL - last_syn_t)/PT1);
    end
    [B, A] = butter(5, [60/4000 200/4000]);
    xl_syn_t_1 = filter(B, A, xl_syn_t);
    %xl_syn_t_1 = 0;
    %xl_syn_t = xl_syn_t - xl_syn_t_1;
    [B, A] = butter(5, [200/4000 600/4000]);
    xl_syn_t_2 = filter(B, A, xl_syn_t);
    %xl_syn_t_2 = 0;
    xl_syn_t = xl_syn_t_2;
    %x1_syn_t = smoothdata(x1_syn_t, 'sgolay', 50);
    x1_syn_t = medfilt1(x1_syn_t, WL * 1000);
    sound(xl_syn_t, fs)
    plot(handles.axes1, xl_syn_t)
    title(handles.axes1, '时域图');
    grid on;
    ysize = size(x1_syn_t);
    y = fft(xl_syn_t, length(xl_syn_t));
    ysize = size(y);
    plot(handles.axes2, abs(y))
    xlabel(handles.axes2, '频率');
    ylabel(handles.axes2, '幅度');
    title(handles.axes2,'频率特征');
    guidata(hObject, handles);
    

% --- Executes on button press in old2young.
function old2young_Callback(hObject, eventdata, handles)
% hObject    handle to old2young (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    FL = 80;
    WL = 240;
    P = 10;
    x1 = handles.y;
    fs = 8000;
    x1 = x1 / max(x1);
    L = length(x1);
    FN = floor(L / FL) - 2;
    exc = zeros(L, 1);
    zi_pre = zeros(P, 1);
    xl_rec = zeros(L, 1);
    zi_rec = zeros(P, 1);
    exc_syn_t = zeros(L, 1);
    x1_syn_t = zeros(L,1);
    last_syn_t = 0;
    zi_syn_t = zeros(P,1);
    hw = hamming(WL);
    for n = 3 : FN
        x1_w = x1(n * FL - WL + 1 : n * FL).*hw;
        [A, E]=lpc(x1_w,P);
        x1_f = x1((n - 1) * FL + 1 : n * FL);
        [excl, zi_pre] = filter(A, 1, x1_f, zi_pre);
        exc((n-1) * FL + 1 : n * FL) = excl;
        [x1_rec1, zi_rec] = filter(1, A, excl, zi_rec);
        x1_rec((n - 1) * FL + 1 : n * FL) = x1_rec1;
        x1_Pitch = exc(n * FL - 222 : n * FL);
        PT = findpitch(x1_Pitch);
        G = sqrt(E * PT);
        PT1 = floor(PT / 1.5);
        poles = roots(A);
        delta = 300 * 2 * pi / 8000;
        for p = 1 : 10
            if imag(poles(p)) > 0
                poles(p) = poles(p) * exp(1i * delta);
            elseif imag(poles(p)) < 0
                poles(p) = poles(p) * exp(-1i * delta);
            end
        end
        A1 = poly(poles);
        tempn_syn_t = (1 : n * FL - last_syn_t);
        exc_synl_t = zeros(length(tempn_syn_t), 1);
        exc_synl_t(mod(tempn_syn_t, PT1)==0) = G;
        exc_synl_t = exc_synl_t((n-1) * FL - last_syn_t + 1 : n * FL - last_syn_t);
        [xl_synl_t, zi_syn_t] = filter(1, A1, exc_synl_t, zi_syn_t);
        exc_syn_t((n - 1) * FL + 1 : n * FL) = exc_synl_t;
        xl_syn_t((n - 1) * FL + 1 : n * FL) = xl_synl_t;
        last_syn_t = last_syn_t + PT1 * floor((n * FL - last_syn_t)/PT1);
    end
    [B, A] = butter(5, [1000/4000 1500/4000]);
    xl_syn_t_1 = filter(B, A, xl_syn_t);
    [B, A] = butter(5, [2200/4000 2500/4000]);
    xl_syn_t_2 = filter(B, A, xl_syn_t);
    %xl_syn_t_2 = 0;
    [B, A] = butter(5, [3300/4000 3500/4000]);
    xl_syn_t_3 = filter(B, A, xl_syn_t);
    %xl_syn_t_3 = 0;
    [B, A] = butter(5, [100/4000 500/4000]);
    xl_syn_t_4 = filter(B, A, xl_syn_t);
    %xl_syn_t_4 = 0;
    xl_syn_t = xl_syn_t_1 + xl_syn_t_2 + xl_syn_t_3 + xl_syn_t_4;
    x1_syn_t = smoothdata(x1_syn_t, 'movmedian', WL);
    sound(xl_syn_t, fs)
    plot(handles.axes1, xl_syn_t)
    title(handles.axes1, '时域图');
    grid on;
    ysize = size(x1_syn_t);
    y = fft(xl_syn_t, length(xl_syn_t));
    ysize = size(y);
    plot(handles.axes2, abs(y))
    xlabel(handles.axes2, '频率');
    ylabel(handles.axes2, '幅度');
    title(handles.axes2,'频率特征');
    guidata(hObject, handles);

% --- Executes on button press in open.
function open_Callback(hObject, eventdata, handles)
% hObject    handle to open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    [filename, pathname] = uigetfile({'*.wav','ALL FILES'}, '选择音频文件');
    str = [pathname filename];
    [temp, Fs] = audioread(str);
    temp1 = resample(temp, 80, 441);
    handles.y1 = temp;
    handles.y = temp1;
    handles.Fs = Fs;
    guidata(hObject, handles);
    
    
function PT = findpitch(s) %低通滤波器、求基音
    [B, A] = butter(5, 700/4000);
    s = filter(B, A, s);
    R = zeros(143, 1);
    for k = 1 : 143
        R(k) = s(144 : 223)' * s(144 - k : 223 - k);
    end
    [R1, T1] = max(R(80: 143));
    T1 = T1 + 79;
    R1 = R1 / (norm(s(144 - T1 : 223 - T1)) + 1);
    [R2, T2] = max(R(40 : 79));
    T2 = T2 + 39;
    R2 = R2 / (norm(s(144 - T2 : 223 - T2)) + 1);
    [R3, T3] = max(R(20: 39));
    T3 = T3 + 19;
    R3 = R3 / (norm(s(144 - T3 : 223 - T3)) + 1);
    Top = T1;
    Rop = R1;
    if R2 >= 0.85 * Rop
        Rop = R2;
        Top = T2;
    end
    if R3 > 0.85 * Rop
        Rop = R3;
        Top = T3;
    end
    PT = Top;
    return


% --- Executes on button press in trump.
function trump_Callback(hObject, eventdata, handles)
% hObject    handle to trump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    FL = 80;
    WL = 240;
    P = 10;
    x1 = handles.y;
    fs = 8000;
    x1 = x1 / max(x1);
    L = length(x1);
    FN = floor(L / FL) - 2;
    exc = zeros(L, 1);
    zi_pre = zeros(P, 1);
    xl_rec = zeros(L, 1);
    zi_rec = zeros(P, 1);
    exc_syn_t = zeros(L, 1);
    x1_syn_t = zeros(L,1);
    last_syn_t = 0;
    zi_syn_t = zeros(P,1);
    hw = hamming(WL);
    for n = 3 : FN
        x1_w = x1(n * FL - WL + 1 : n * FL).*hw;
        [A, E]=lpc(x1_w,P);
        x1_f = x1((n - 1) * FL + 1 : n * FL);
        [excl, zi_pre] = filter(A, 1, x1_f, zi_pre);
        exc((n-1) * FL + 1 : n * FL) = excl;
        [x1_rec1, zi_rec] = filter(1, A, excl, zi_rec);
        x1_rec((n - 1) * FL + 1 : n * FL) = x1_rec1;
        x1_Pitch = exc(n * FL - 222 : n * FL);
        PT = findpitch(x1_Pitch);
        G = sqrt(E * PT);
        PT1 = floor(PT / 0.85);
        poles = roots(A);
        delta = 30 * 2 * pi / 8000;
        for p = 1 : 10
            if imag(poles(p)) > 0
                poles(p) = poles(p) * exp(1i * delta);
            elseif imag(poles(p)) < 0
                poles(p) = poles(p) * exp(-1i * delta);
            end
        end
        A1 = poly(poles);
        tempn_syn_t = (1 : n * FL - last_syn_t);
        exc_synl_t = zeros(length(tempn_syn_t), 1);
        exc_synl_t(mod(tempn_syn_t, PT1)==0) = G;
        exc_synl_t = exc_synl_t((n-1) * FL - last_syn_t + 1 : n * FL - last_syn_t);
        [xl_synl_t, zi_syn_t] = filter(1, A1, exc_synl_t, zi_syn_t);
        exc_syn_t((n - 1) * FL + 1 : n * FL) = exc_synl_t;
        xl_syn_t((n - 1) * FL + 1 : n * FL) = xl_synl_t;
        last_syn_t = last_syn_t + PT1 * floor((n * FL - last_syn_t)/PT1);
    end
    [B, A] = butter(5, [60/4000 200/4000]);
    xl_syn_t_1 = filter(B, A, xl_syn_t);
    %xl_syn_t_1 = 0;
    %xl_syn_t = xl_syn_t - xl_syn_t_1;
    [B, A] = butter(5, [200/4000 600/4000]);
    xl_syn_t_2 = filter(B, A, xl_syn_t);
    %xl_syn_t_2 = 0;
    %xl_syn_t = xl_syn_t_2;
    %x1_syn_t = smoothdata(x1_syn_t, 'sgolay', 50);
    x1_syn_t = medfilt1(x1_syn_t, WL * 1000);
    sound(xl_syn_t, 8500)
    plot(handles.axes1, xl_syn_t)
    title(handles.axes1, '时域图');
    grid on;
    ysize = size(x1_syn_t);
    y = fft(xl_syn_t, length(xl_syn_t));
    ysize = size(y);
    plot(handles.axes2, abs(y))
    xlabel(handles.axes2, '频率');
    ylabel(handles.axes2, '幅度');
    title(handles.axes2,'频率特征');
    guidata(hObject, handles);
    
    % --- Executes on button press in trump.
function xx_Callback(hObject, eventdata, handles)
% hObject    handle to trump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    FL = 80;
    WL = 240;
    P = 10;
    x1 = handles.y;
    fs = 8000;
    x1 = x1 / max(x1);
    L = length(x1);
    FN = floor(L / FL) - 2;
    exc = zeros(L, 1);
    zi_pre = zeros(P, 1);
    xl_rec = zeros(L, 1);
    zi_rec = zeros(P, 1);
    exc_syn_t = zeros(L, 1);
    x1_syn_t = zeros(L,1);
    last_syn_t = 0;
    zi_syn_t = zeros(P,1);
    hw = hamming(WL);
    for n = 3 : FN
        x1_w = x1(n * FL - WL + 1 : n * FL).*hw;
        [A, E]=lpc(x1_w,P);
        x1_f = x1((n - 1) * FL + 1 : n * FL);
        [excl, zi_pre] = filter(A, 1, x1_f, zi_pre);
        exc((n-1) * FL + 1 : n * FL) = excl;
        [x1_rec1, zi_rec] = filter(1, A, excl, zi_rec);
        x1_rec((n - 1) * FL + 1 : n * FL) = x1_rec1;
        x1_Pitch = exc(n * FL - 222 : n * FL);
        PT = findpitch(x1_Pitch);
        G = sqrt(E * PT);
        PT1 = floor(PT / 0.5);
        poles = roots(A);
        delta = 60 * 2 * pi / 8000;
        for p = 1 : 10
            if imag(poles(p)) > 0
                poles(p) = poles(p) * exp(1i * delta);
            elseif imag(poles(p)) < 0
                poles(p) = poles(p) * exp(-1i * delta);
            end
        end
        A1 = poly(poles);
        tempn_syn_t = (1 : n * FL - last_syn_t);
        exc_synl_t = zeros(length(tempn_syn_t), 1);
        exc_synl_t(mod(tempn_syn_t, PT1)==0) = G;
        exc_synl_t = exc_synl_t((n-1) * FL - last_syn_t + 1 : n * FL - last_syn_t);
        [xl_synl_t, zi_syn_t] = filter(1, A1, exc_synl_t, zi_syn_t);
        exc_syn_t((n - 1) * FL + 1 : n * FL) = exc_synl_t;
        xl_syn_t((n - 1) * FL + 1 : n * FL) = xl_synl_t;
        last_syn_t = last_syn_t + PT1 * floor((n * FL - last_syn_t)/PT1);
    end
    [B, A] = butter(5, [60/4000 200/4000]);
    xl_syn_t_1 = filter(B, A, xl_syn_t);
    %xl_syn_t_1 = 0;
    %xl_syn_t = xl_syn_t - xl_syn_t_1;
    [B, A] = butter(5, [200/4000 600/4000]);
    xl_syn_t_2 = filter(B, A, xl_syn_t);
    %xl_syn_t_2 = 0;
    xl_syn_t = xl_syn_t_2;
    x1_syn_t = medfilt1(x1_syn_t, WL * 1000);
    sound(xl_syn_t, 8500)
    plot(handles.axes1, xl_syn_t)
    title(handles.axes1, '时域图');
    grid on;
    ysize = size(x1_syn_t);
    y = fft(xl_syn_t, length(xl_syn_t));
    ysize = size(y);
    plot(handles.axes2, abs(y))
    xlabel(handles.axes2, '频率');
    ylabel(handles.axes2, '幅度');
    title(handles.axes2,'频率特征');
    guidata(hObject, handles);
