function varargout = GUI_beam(varargin)
% GUI_BEAM MATLAB code for GUI_beam.fig
%      GUI_BEAM, by itself, creates a new GUI_BEAM or raises the existing
%      singleton*.
%
%      H = GUI_BEAM returns the handle to a new GUI_BEAM or the handle to
%      the existing singleton*.
%
%      GUI_BEAM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_BEAM.M with the given input arguments.
%
%      GUI_BEAM('Property','Value',...) creates a new GUI_BEAM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_beam_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_beam_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_beam

% Last Modified by GUIDE v2.5 15-Apr-2022 23:21:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_beam_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_beam_OutputFcn, ...
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


% --- Executes just before GUI_beam is made visible.
function GUI_beam_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_beam (see VARARGIN)

% Choose default command line output for GUI_beam
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes GUI_beam wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = GUI_beam_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function E_Callback(hObject, eventdata, handles)
% hObject    handle to E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of E as text
%        str2double(get(hObject,'String')) returns contents of E as a double

global E

handles.E = get(hObject,'String');
guidata(hObject,handles)

E = str2double(handles.E);


% --- Executes during object creation, after setting all properties.
function E_CreateFcn(hObject, eventdata, handles)
% hObject    handle to E (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function A_Callback(hObject, eventdata, handles)
% hObject    handle to A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of A as text
%        str2double(get(hObject,'String')) returns contents of A as a double
global A
A = str2double(get(hObject, 'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function A_CreateFcn(hObject, eventdata, handles)
% hObject    handle to A (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Nut_Callback(hObject, eventdata, handles)
% hObject    handle to Nut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nut as text
%        str2double(get(hObject,'String')) returns contents of Nut as a double

global nut
handles.Nut = get(hObject,'String');
guidata(hObject,handles)

handles.Nut(1) = [];
handles.Nut(length(handles.Nut)) = [];

k = str2double(regexp(handles.Nut,'[+-]?\d+\.?\d*','match'));



nut = zeros(length(k)/2,2);

for i = 1:size(nut,1)
        nut(i,1) = nut(i,1)+k(2*i-1);
        nut(i,2) = nut(i,2)+k(2*i);
end



% --- Executes during object creation, after setting all properties.
function Nut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Phan_tu_Callback(hObject, eventdata, handles)
% hObject    handle to Phan_tu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Phan_tu as text
%        str2double(get(hObject,'String')) returns contents of Phan_tu as a double
global ptu
handles.Phan_tu = get(hObject,'String');
guidata(hObject,handles)

handles.Phan_tu(1) = [];
handles.Phan_tu(length(handles.Phan_tu)) = [];

k = str2double(regexp(handles.Phan_tu,'[+-]?\d+\.?\d*','match'));



ptu = zeros(length(k)/2,2);

for i = 1:size(ptu,1)
        ptu(i,1) = ptu(i,1)+k(2*i-1);
        ptu(i,2) = ptu(i,2)+k(2*i);
end

% --- Executes during object creation, after setting all properties.
function Phan_tu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Phan_tu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rangbuoc_Callback(hObject, eventdata, handles)
% hObject    handle to rangbuoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rangbuoc as text
%        str2double(get(hObject,'String')) returns contents of rangbuoc as a double
global rb_u
handles.rbuoc_u = get(hObject,'String');
guidata(hObject,handles)

handles.rbuoc_u(1) = [];
handles.rbuoc_u(length(handles.rbuoc_u)) = [];

rb_u = str2double(regexp(handles.rbuoc_u,'[+-]?\d+\.?\d*','match'));


% --- Executes during object creation, after setting all properties.
function rangbuoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rangbuoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VOU_Callback(hObject, eventdata, handles)
% hObject    handle to VOU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VOU as text
%        str2double(get(hObject,'String')) returns contents of VOU as a double
global valueOfU
handles.vou = get(hObject,'String');
guidata(hObject,handles)

handles.vou(1) = [];
handles.vou(length(handles.vou)) = [];

valueOfU = str2double(regexp(handles.vou,'[+-]?\d+\.?\d*','match'));


% --- Executes during object creation, after setting all properties.
function VOU_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VOU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function rangbuocF_Callback(hObject, eventdata, handles)
% hObject    handle to rangbuocF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rangbuocF as text
%        str2double(get(hObject,'String')) returns contents of rangbuocF as a double
global rb_F
handles.rbuoc_F = get(hObject,'String');
guidata(hObject,handles)

handles.rbuoc_F(1) = [];
handles.rbuoc_F(length(handles.rbuoc_F)) = [];

rb_F = str2double(regexp(handles.rbuoc_F,'[+-]?\d+\.?\d*','match'));

% --- Executes during object creation, after setting all properties.
function rangbuocF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rangbuocF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function VOF_Callback(hObject, eventdata, handles)
% hObject    handle to VOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of VOF as text
%        str2double(get(hObject,'String')) returns contents of VOF as a double
global valueOfF
handles.vof = get(hObject,'String');
guidata(hObject,handles)

handles.vof(1) = [];
handles.vof(length(handles.vof)) = [];

valueOfF = str2double(regexp(handles.vof,'[+-]?\d+\.?\d*','match'));

% --- Executes during object creation, after setting all properties.
function VOF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VOF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Nut Phan_tu mt_bool_phantu Utt

global E A J
global nut
global ptu
global rb_u
global valueOfU
global rb_F
global valueOfF
global phantu_ux taitrong_ux
global phantu_kn taitrong_kn

BTD_nut = 3;
    
Nut = nut;
Phan_tu = ptu;
rangbuoc_u = rb_u;
VOU = valueOfU;
rangbuoc_F = rb_F;
VOF = valueOfF;

sophantu = size(Phan_tu,1);
sonut_pt = size(Phan_tu,2);
so_nut = size(Nut,1);
so_BTD_he = so_nut*BTD_nut;

% VE MO HINH BAI TOAN

axes(handles.mohinh1)
hold on
% axis([min(Nut(:,1))-0.2 max(Nut(:,1))+0.2 min(Nut(:,2))-0.2 max(Nut(:,2))+0.2])
axis padded
grid on
for i = 1:length(rangbuoc_u)
    m = mod(rangbuoc_u(i),3);
    if m == 0
        x = Nut(rangbuoc_u(i)/3,1);
        y = Nut(rangbuoc_u(i)/3,2);
        plot(x,y,'r+','markersize',20,'linewidth',3);
        plot(x,y,'ro','markersize',15,'linewidth',2);
    end
    m = mod(rangbuoc_u(i)+1,3);
    if m == 0
        x = Nut((rangbuoc_u(i)+1)/3,1);
        y = Nut((rangbuoc_u(i)+1)/3,2);
        plot(x,y-0.05,'^','markersize',17,'markerfacecolor','g');
    end
    m = mod(rangbuoc_u(i)+2,3);
    if m == 0
        x = Nut((rangbuoc_u(i)+2)/3,1);
        y = Nut((rangbuoc_u(i)+2)/3,2); 
        plot(x-0.05,y,'>','markersize',17,'markerfacecolor','b');
    end
end

for i = 1 : sophantu
    nutdau = Phan_tu(i,1);
    nutcuoi = Phan_tu(i,2);
    mt_bool_phantu(i,:) = tinh_bool([nutdau nutcuoi], BTD_nut);
    xi = Nut(nutdau,1);
    xj = Nut(nutcuoi,1);
    yi = Nut(nutdau,2);
    yj = Nut(nutcuoi,2);
    L(i) = sqrt((xj-xi)^2 + (yj-yi)^2);%chieu dai phan tu
    C(i) = (xj-xi)/L(i);%cos(alpha) cua phan tu
    S(i) = (yj-yi)/L(i);%sin(alpha) cua phan tu
    plot([xi xj],[yi yj],'-','linewidth',3);
    plot([xi xj],[yi yj],'o','markersize',3);
    text(xi+0.03,yi,num2str(nutdau));
    text(xj+0.03,yj,num2str(nutcuoi));
    if xj>xi
        text((xj-xi)/2+xi,(yj-yi)/2 + yi,num2str(i));
    elseif xj == xi
        text(xj,(yj-yi)/2 + yi,num2str(i));
    elseif xj<xi
        text((xi-xj)/2+xj,(yj-yi)/2 + yi,num2str(i));
    end
end


% Dieu kien bien rang buoc ve tai trong (boundary condition) (gom ca tai trong)
Ftt = zeros(so_BTD_he,1);

for i=1:size(rangbuoc_F,2)
    Ftt(rangbuoc_F(:,i),1) = Ftt(rangbuoc_F(:,i),1)+VOF(1,i);
end

F_quydoi = zeros(sophantu, sonut_pt*BTD_nut);

% quy doi luc uon
for i=1:size(phantu_ux,2)
    index = phantu_ux(i);
    qd = quydoi_Fuon_PBD(taitrong_ux(i),C(index),S(index),L(index));
    
    for j=1:size(qd,1)
        F_quydoi(index,j) = F_quydoi(index,j)+qd(j,1);
    end
end

% quy doi luc keo nen
for i=1:size(phantu_kn,2)
    index = phantu_kn(i);
    qd = quydoi_Fdoc_PBD(taitrong_kn(i),C(index),S(index),L(index));
    
    for j=1:size(qd,1)
        F_quydoi(index,j) = F_quydoi(index,j)+qd(j,1);
    end
end

% combine vao vector tai tong the
for i=1:sophantu
    for j=1:(sonut_pt*BTD_nut)
        Ftt(mt_bool_phantu(i,j),1) = Ftt(mt_bool_phantu(i,j),1)+F_quydoi(i,j);
    end 
end

% XAY DUNG MA TRAN CUNG TONG THE
Ktt = zeros(so_BTD_he,so_BTD_he);
for i=1:sophantu
    Kpt = tinh_Kpt(E,A,J,L(i),S(i),C(i));
    Ktt = tinh_Ktt(Ktt,Kpt,mt_bool_phantu(i,:));
end

% KHU DIEU KIEN BIEN GIAI HE PHUONG TRINH
[Ktt_bd, Ftt_bd] = khu_dkb(Ktt,Ftt,rangbuoc_u);
Utt = Ktt_bd\Ftt_bd;

% TINH PHAN LUC LIEN KET
phanluc = zeros(so_BTD_he,1);
F_pl = Ktt*Utt-Ftt;
for i = 1:length(rangbuoc_u)
    phanluc(rangbuoc_u(i),1) = phanluc(rangbuoc_u(i),1)+F_pl(rangbuoc_u(i),1);
end

% TINH MOMENT UON
for i = 1:sophantu
    Se = mt_tinh_noiluc(E,J,L(i),C(i),S(i));
    M(:,i) = Se*Utt(mt_bool_phantu(i,:),1);
end


% VE BIEU DO MOMENT
Gmax = max(abs(Nut(:)))*0.2;

Mmax = max(abs(M(:)));
Me = M/Mmax*Gmax;

axes(handles.mohinh3)
axis padded
hold on; 
h1=trisurf(Phan_tu, Nut(:,1),Nut(:,2), zeros(length(Nut),1));
set(h1, 'facecolor', 'none', 'edgecolor', 'k', 'linewidth' ,3);

for i = 1:sophantu 
    nd = Phan_tu(i,:); % cac nut cua phan tu thu i 
    % xac dinh toa do cua cac nut cua phan tu thu i
    x = Nut(nd,1); 
    y = Nut(nd,2);
    Mve = Me(:,i);
    if x(1) == x(2)
        X = x(1) + [0 Mve(1) Mve(2) 0];
        Y = [y(1) y(1) y(2) y(2)];
        C = [M(1,i) M(1,i) M(2,i) M(2,i)];
        fill(X,Y,C, 'facealpha', 0.8)
        hold on
    end
    
    if y(1) == y(2)
        X = [x(1) x(1) x(2) x(2)];
        Y = y(1) + [0 Mve(1) Mve(2) 0];
        C = [M(1,i) M(1,i) M(2,i) M(2,i)];
        fill(X,Y,C, 'facealpha', 0.8)
        hold on
    end
    
    colorbar
    colormap(jet(100))
end

% XUAT DU LIEU RA BANG
set(handles.chuyenvi, 'Data', Utt)
set(handles.noiluc, 'Data', transpose(M))
set(handles.phanluc, 'Data', phanluc)
% set(handles.ungsuat, 'Data', ungSuat)

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in chuyenvi.
function chuyenvi_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to chuyenvi (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles)


% --- Executes when entered data in editable cell(s) in noiluc.
function noiluc_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to noiluc (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles)


% --- Executes when entered data in editable cell(s) in phanluc.
function phanluc_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to phanluc (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles)


% --- Executes when entered data in editable cell(s) in ungsuat.
function ungsuat_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to ungsuat (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject,handles)



function J_Callback(hObject, eventdata, handles)
% hObject    handle to J (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of J as text
%        str2double(get(hObject,'String')) returns contents of J as a double
global J
J = str2double(get(hObject, 'String'));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function J_CreateFcn(hObject, eventdata, handles)
% hObject    handle to J (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% HAM TINH MA TRAN BOOL PHAN TU
function [mt_bool_pt] = tinh_bool(mt_nut_pt,btd_nut)
    k = 0;
    sonut_pt = size(mt_nut_pt,2);
    for i = 1:sonut_pt
        s = (mt_nut_pt(i) - 1)*btd_nut;
        for j = 1:btd_nut
            k=k+1;
            mt_bool_pt(k) = s + j;
        end
    end


% HAM QUY DOI LUC
function luc_uon_phanbodeu = quydoi_Fuon_PBD(q0,C,S,L)
    luc_uon_phanbodeu(1,1) = (L/2)*-q0*S;
    luc_uon_phanbodeu(2,1) = (L/2)*q0*C;
    luc_uon_phanbodeu(3,1) = (L/2)*q0*L/6;
    luc_uon_phanbodeu(4,1) = (L/2)*-q0*S;
    luc_uon_phanbodeu(5,1) = (L/2)*q0*C;
    luc_uon_phanbodeu(6,1) = (L/2)*-q0*L/6;

   
function luc_doc_phanbodeu = quydoi_Fdoc_PBD(p0,C,S,L)
    luc_doc_phanbodeu(1,1) = (L/2)*p0*C;
    luc_doc_phanbodeu(2,1) = (L/2)*p0*S;
    luc_doc_phanbodeu(3,1) = (L/2)*0;
    luc_doc_phanbodeu(4,1) = (L/2)*p0*C;
    luc_doc_phanbodeu(5,1) = (L/2)*p0*S;
    luc_doc_phanbodeu(6,1) = (L/2)*0;


function luc_uon_taptrung = quydoi_Fuon_TT(Q,C,S,L)
    luc_uon_taptrung(1,1) = -Q*S/2;
    luc_uon_taptrung(2,1) = Q*C/2;
    luc_uon_taptrung(3,1) = Q*L/8;
    luc_uon_taptrung(4,1) = -Q*S/2;
    luc_uon_taptrung(5,1) = Q*C/2;
    luc_uon_taptrung(6,1) = -Q*L/8;


function moment_taptrung = quydoi_moment_TT(M,C,S,L)
    moment_taptrung(1,1) = 3*M*S/(2*L);
    moment_taptrung(2,1) = -3*M*C/(2*L);
    moment_taptrung(3,1) = -M/4;
    moment_taptrung(4,1) = -3*M*S/(2*L);
    moment_taptrung(5,1) = 3*M*C/(2*L);
    moment_taptrung(6,1) = -M/4;


% HAM TINH MA TRAN CUNG PHAN TU
function [Kpt] = tinh_Kpt(E,A,J,L,S,C)
    B=12*J/(L^2);
    Kpt = (E/L)*[A*C^2+B*S^2    (A-B)*C*S       -B*L*S/2    -(A*C^2+B*S^2)  -(A-B)*C*S      -B*L*S/2;
                (A-B)*C*S        A*S^2+B*C^2     B*L*C/2    -(A-B)*C*S      -(A*S^2+B*C^2)   B*L*C/2;
                -B*L*S/2         B*L*C/2         4*J         B*L*S/2        -B*L*C/2         2*J;
                -(A*C^2+B*S^2) -(A-B)*C*S        B*L*S/2     A*C^2+B*S^2    (A-B)*C*S        B*L*S/2;
                -(A-B)*C*S     -(A*S^2+B*C^2)   -B*L*C/2    (A-B)*C*S        A*S^2+B*C^2    -B*L*C/2;
                -B*L*S/2         B*L*C/2         2*J         B*L*S/2        -B*L*C/2         4*J];


% HAM TINH MA TRAN CUNG TONG THE
function [Ktt] = tinh_Ktt(Ktt,Kpt,mt_bool_phantu)
    Ktt(mt_bool_phantu,mt_bool_phantu) = Ktt(mt_bool_phantu,mt_bool_phantu)+Kpt;


function [Ktt, Ftt] = khu_dkb(Ktt,Ftt,rangbuoc)
    n = length(rangbuoc);
    btdtt = length(Ktt);
    for i = 1:n
        c = rangbuoc(i);
        for j = 1:btdtt
            Ktt(c,j) = 0; %khu hang
            Ktt(j,c) = 0; %khu cot
        end
        Ktt(c,c) = 1;
        Ftt(c) = 0;
    end


% HAM TINH MA TRAN TINH NOI LUC
function [Se] = mt_tinh_noiluc(E,J,L,C,S)
    Se = E*(J/(L^3))*[6*L*S  -6*L*C  -4*L^2   -6*L*S     6*L*C   -2*L^2;
                     -6*L*S  6*L*C   2*L^2    6*L*S     -6*L*C   4*L^2];



function ptu_ux_Callback(hObject, eventdata, handles)
% hObject    handle to ptu_ux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptu_ux as text
%        str2double(get(hObject,'String')) returns contents of ptu_ux as a double
global phantu_ux
handles.ptu_ux = get(hObject,'String');
guidata(hObject,handles)

handles.ptu_ux(1) = [];
handles.ptu_ux(length(handles.ptu_ux)) = [];

phantu_ux = str2double(regexp(handles.ptu_ux,'[+-]?\d+\.?\d*','match'));



% --- Executes during object creation, after setting all properties.
function ptu_ux_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptu_ux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tt_ux_Callback(hObject, eventdata, handles)
% hObject    handle to tt_ux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tt_ux as text
%        str2double(get(hObject,'String')) returns contents of tt_ux as a double
global taitrong_ux
handles.tt_ux = get(hObject,'String');
guidata(hObject,handles)

handles.tt_ux(1) = [];
handles.tt_ux(length(handles.tt_ux)) = [];

taitrong_ux = str2double(regexp(handles.tt_ux,'[+-]?\d+\.?\d*([eE][+-]?\d+)?','match'));

% --- Executes during object creation, after setting all properties.
function tt_ux_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tt_ux (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ptu_kn_Callback(hObject, eventdata, handles)
% hObject    handle to ptu_kn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ptu_kn as text
%        str2double(get(hObject,'String')) returns contents of ptu_kn as a double
global phantu_kn
handles.ptu_kn = get(hObject,'String');
guidata(hObject,handles)

handles.ptu_kn(1) = [];
handles.ptu_kn(length(handles.ptu_kn)) = [];

phantu_kn = str2double(regexp(handles.ptu_kn,'[+-]?\d+\.?\d*','match'));

% --- Executes during object creation, after setting all properties.
function ptu_kn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ptu_kn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tt_kn_Callback(hObject, eventdata, handles)
% hObject    handle to tt_kn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tt_kn as text
%        str2double(get(hObject,'String')) returns contents of tt_kn as a double
global taitrong_kn
handles.tt_kn = get(hObject,'String');
guidata(hObject,handles)

handles.tt_kn(1) = [];
handles.tt_kn(length(handles.tt_kn)) = [];

taitrong_kn = str2double(regexp(handles.tt_kn,'[+-]?\d+\.?\d*([eE][+-]?\d+)?','match'));

% --- Executes during object creation, after setting all properties.
function tt_kn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tt_kn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function scale_Callback(hObject, eventdata, handles)
% hObject    handle to scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scale as text
%        str2double(get(hObject,'String')) returns contents of scale as a double
handles.scale = str2double(get(hObject,'String'));
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function scale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in animate.
function animate_Callback(hObject, eventdata, handles)
% hObject    handle to animate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% MO PHONG QUA TRINH BIEN DANG

global Nut Phan_tu mt_bool_phantu Utt

axes(handles.mohinh2)

scale = handles.scale;
N = 20;

for i = 1:N
    hold off
    for j = 1:size(Phan_tu,1)
        xi = Nut(Phan_tu(j,1),1);% nut i (nut dau)
        yi = Nut(Phan_tu(j,1),2); 

        xj = Nut(Phan_tu(j,2),1);%nut j (nut cuoi)
        yj = Nut(Phan_tu(j,2),2);
        
        plot([xi xj],[yi yj],'-','linewidth',3);
        plot([xi xj],[yi yj],'o','markersize',5);
        
        hold on
        
        upt = Utt(mt_bool_phantu(j,:));
        plot([xi+scale*(i/N)*upt(1) xj+scale*(i/N)*upt(4)], [yi+scale*(i/N)*upt(2) yj+scale*(i/N)*upt(5)], '-ob', 'markerfacecolor', [1 0 0], 'markersize', 3, 'linewidth', 2)% Co he sau khi bi dien dang
%         axis([min(Nut(:,1))-0.2 max(Nut(:,1))+0.2 min(Nut(:,2))-0.2 max(Nut(:,2))+0.2])
        axis padded
    end
    getframe
end
