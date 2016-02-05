function varargout = histoGUI(varargin)
% HISTOGUI MATLAB code for histoGUI.fig
%      HISTOGUI, by itself, creates a new HISTOGUI or raises the existing
%      singleton*.
%
%      H = HISTOGUI returns the handle to a new HISTOGUI or the handle to
%      the existing singleton*.
%
%      HISTOGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HISTOGUI.M with the given input arguments.
%
%      HISTOGUI('Property','Value',...) creates a new HISTOGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before histoGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to histoGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help histoGUI

% Last Modified by GUIDE v2.5 11-Jul-2012 16:28:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @histoGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @histoGUI_OutputFcn, ...
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


% --- Executes just before histoGUI is made visible.
function histoGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to histoGUI (see VARARGIN)

set(handles.bildVorschau,'xtick',[]);
set(handles.bildVorschau,'ytick',[]);
set(handles.extraktionVorschau,'xtick',[]);
set(handles.extraktionVorschau,'ytick',[]);
set(handles.vermessungZeigeVermessung,'xtick',[]);
set(handles.vermessungZeigeVermessung,'ytick',[]);
set(handles.trabekelZeigeBild,'xtick',[]);
set(handles.trabekelZeigeBild,'ytick',[]);
set(handles.einzeltrabekelZeige,'xtick',[]);
set(handles.einzeltrabekelZeige,'ytick',[]);

%set(handles.axesInvis1,'xtick',[]);
%set(handles.axesInvis1,'ytick',[]);

handles.wirbel1 = 0;
handles.wirbel2 = 0;
handles.wirbel3 = 0;
handles.wirbel1Seg = 0;
handles.wirbel2Seg = 0;

handles.pixelgroesse = 0.25;

%Daten
handles.ergebniszaehler = 1;
handles.einstellungenZaehler = 1;
handles.umfang = 0;
handles.volumen = 0;
handles.wirbelHoehe = 0;
handles.wirbelHoeheStart = 0;
handles.wirbelHoeheEnde = 0;
handles.wirbelBreiteOben = 0;
handles.wirbelBreiteUnten = 0;
handles.calcifizierteFlaeche = 0;
handles.cortikalisflaeche = 0;
handles.cortikalisdicke = 0;
handles.trabekelflaeche = 0;
handles.trabekelListe = {};
handles.trabekelDurchschnittsabstand = 0;
handles.trabekelDurchschnittdicke = 0;
handles.pfad = 0;
handles.trabekelProPixel = 0;

handles.trabekelStep = 0;

%Einstellungen fuer Extraktion
handles.einstellungExtraktionNU = 0.3;
handles.einstellungExtraktionMU = 1;
handles.einstellungExtraktionDT = 1;
handles.einstellungExtraktionIterations = 50;


% Choose default command line output for histoGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes histoGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = histoGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in schritteBilderLaden.
function schritteBilderLaden_Callback(hObject, eventdata, handles)
    
hidePanels(hObject, eventdata, handles);
set(handles.panelBilderLaden,'Visible','On');
guidata(hObject, handles);

% --- Executes on button press in bildLaden.
function bildLaden_Callback(hObject, eventdata, handles)
    
    if (handles.pfad)
        [FileName,PathName] = uigetfile([handles.pfad,'*.*'],'Select an image file');
    else
        [FileName,PathName] = uigetfile('*.*','Select an image file');
    end
    
    handles.pfad = PathName;
    
    if (~FileName)
        return;
    end
    
    bild = im2double(rgb2gray(imread([PathName,FileName])));
    
    if (get(handles.bildTyp,'Value') == 1)
        handles.wirbel1 = bild;
        
    elseif (get(handles.bildTyp,'Value') == 2)
        handles.wirbel2 = bild;
    else
        handles.wirbel3 = bild;
    end

    axes(handles.bildVorschau);
    imagesc(bild);colormap(gray);set(handles.bildVorschau,'xtick',[]);set(handles.bildVorschau,'ytick',[]);axis equal;
    set(handles.wirbel1VorschauText,'String',['Groesse: ',num2str(size(bild,2)),'x',num2str(size(bild,1))]);

guidata(hObject, handles);



function hidePanels(hObject, eventdata, handles)
    set(handles.panelBilderLaden,'Visible','Off');
    set(handles.panelExtraktion,'Visible','Off');
    set(handles.panelVermessung,'Visible','Off');
    set(handles.panelTrabekel,'Visible','Off');
    set(handles.panelEinzeltrabekel,'Visible','Off');
    set(handles.panelTrabekelanalyse,'Visible','Off');
    guidata(hObject, handles);


function schritteExtraktion_Callback(hObject, eventdata, handles)
    hidePanels(hObject, eventdata, handles);
    set(handles.panelExtraktion,'Visible','On');
    handles = loeschePanelEinstellungen(hObject, eventdata, handles);
    handles=extraktionTyp_Callback(hObject, eventdata, handles);
    
    guidata(hObject, handles);

function schritteVermessung_Callback(hObject, eventdata, handles)
    hidePanels(hObject, eventdata, handles);
    set(handles.panelVermessung,'Visible','On');
    handles = loeschePanelEinstellungen(hObject, eventdata, handles);
    
    %Data for histogram
    [n,xout] = hist(handles.wirbel1(:),100);
    
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Angle start','-20','vermessungWinkelbereich1'},{n,xout});
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Angle end','20','vermessungWinkelbereich2'},{n,xout});
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Distance top','0.1','vermessungAbstandOben'},{n,xout});
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Distance bottom','0.1','vermessungAbstandUnten'},{n,xout});
    

    
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Threshold calc. area','0.25','knochendatenSchwellwertCalc'},{n,xout});
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Threshold cort. area','0.25','knochendatenSchwellwertCort'},{n,xout});
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Shift for inner segmentation','0','knochendatenInnerSegOffset'},{n,xout});
    
    guidata(hObject, handles);
    
function schritteTrabekel_Callback(hObject, eventdata, handles)
    hidePanels(hObject, eventdata, handles);
    set(handles.panelTrabekel,'Visible','On');
    handles = loeschePanelEinstellungen(hObject, eventdata, handles);
    
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Variance Gaussian','1','trabekelVarianzGaussglaettung'});
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Threshold size','20','trabekelSchwellwertTrabekelgroesse'});
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Threshold bone','0.25','trabekelSchwellwertKnochen'});
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Shift for inner segmentation','0','knochendatenInnerSegOffset'});
    
    handles.trabekelStep = 1;
    
    guidata(hObject, handles);
    
function schritteEinzeltrabekel_Callback(hObject, eventdata, handles)
    if (numel(handles.trabekelListe) == 0)
        warndlg('Keine Trabekelliste gefunden!')
        return;
    end


    hidePanels(hObject, eventdata, handles);
    set(handles.panelEinzeltrabekel,'Visible','On');
    handles = loeschePanelEinstellungen(hObject, eventdata, handles);
    
    %Slider initialisieren
    set(handles.einzeltrabekelSlider,'Visible','On');
    
    set(handles.einzeltrabekelSlider,'Min',1);
    set(handles.einzeltrabekelSlider,'Max',numel(handles.trabekelListe));
    set(handles.einzeltrabekelSlider,'Value',1);
    set(handles.einzeltrabekelSlider,'SliderStep',[1/(numel(handles.trabekelListe)-1) 1]);
    
    handles=zeigeEinzeltrabekel(hObject, eventdata, handles,1,handles.einzeltrabekelZeige,[]);
    
    handles=einzeltrabekelSlider_Callback(hObject, eventdata, handles);
    
    handles.trabekelStep = 2;
    
    guidata(hObject, handles);
    
function schritteTrabekelanalyse_Callback(hObject, eventdata, handles)
    if (numel(handles.trabekelListe) == 0)
        warndlg('Keine Trabekelliste gefunden!')
        return;
    end
    
    hidePanels(hObject, eventdata, handles);
    set(handles.panelTrabekelanalyse,'Visible','On');
    handles = loeschePanelEinstellungen(hObject, eventdata, handles);
    
    guidata(hObject, handles);

    
% --- Executes on selection change in extraktionTyp.
function handles=extraktionTyp_Callback(hObject, eventdata, handles)
handles = loeschePanelEinstellungen(hObject, eventdata, handles);
handles = erzeugeEinstellung(hObject, eventdata, handles,{'nu','0.3','einstellungNU'});
handles = erzeugeEinstellung(hObject, eventdata, handles,{'mu','1','einstellungMU'});
handles = erzeugeEinstellung(hObject, eventdata, handles,{'dt','1','einstellungDT'});
handles = erzeugeEinstellung(hObject, eventdata, handles,{'Iterations','50','einstellungIterationen'});

if (get(handles.extraktionTyp,'Value') == 1)
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Scale start','8','einstellungScaleStart'});
else
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'Scale start','1','einstellungScaleStart'});
end



handles = erzeugeEinstellung(hObject, eventdata, handles,{'Scale end','1','einstellungScaleEnd'});
handles = erzeugeEinstellung(hObject, eventdata, handles,{'Tolerance','0.1','einstellungTolerance'});

if (get(handles.extraktionTyp,'Value') == 2)
    handles = erzeugeEinstellung(hObject, eventdata, handles,{'InnerShift','5','einstellungInnerShift'});
end

% --- Executes on button press in extraktionReset.
function extraktionReset_Callback(hObject, eventdata, handles)
    if (get(handles.extraktionTyp,'Value') == 1)
        handles.wirbel1Seg = 0;
    elseif (get(handles.extraktionTyp,'Value') == 2)
        handles.wirbel2Seg = 0;
    end
    
    guidata(hObject, handles);

% --- Executes on button press in extraktion1StartenButton.
function extraktion1StartenButton_Callback(hObject, eventdata, handles)
    
if (handles.wirbel1 == 0)
    warndlg('Bitte zunächst ein Bild laden!')
    return;
end

set(handles.extraktionStop,'UserData',0);

if (get(handles.extraktionTyp,'Value') == 1)
    if (handles.wirbel1Seg ~= 0)
        seg = handles.wirbel1Seg;
    else
        seg = zeros(size(handles.wirbel1));
        seg(2:end-2,2:end-2) = -1;
        seg = double(bwdist(seg < 0) - bwdist(1-(seg < 0)));
    end
    scale = 1/str2double(leseEinstellung(handles,'einstellungScaleStart'));
    
elseif (get(handles.extraktionTyp,'Value') == 2)
    if (handles.wirbel1Seg == 0)
        warndlg('Bitte zunächst die erste Extraktion durchführen!')
        return;
    end
    
    if (handles.wirbel2Seg ~= 0)
        seg = handles.wirbel2Seg;
    else
        seg = handles.wirbel1Seg;
        seg = seg + str2double(leseEinstellung(handles,'einstellungInnerShift'));
    end
    scale = 1/str2double(leseEinstellung(handles,'einstellungScaleStart'));
end



%profile on;
for i=1:999
    if (get(handles.extraktionTyp,'Value') == 1)
        handles.wirbel1Seg = seg;
    elseif (get(handles.extraktionTyp,'Value') == 2)
        handles.wirbel2Seg = seg;
    end
    
    guidata(hObject, handles);
    if (get(handles.extraktionStop,'UserData') == 0)

        showImageWithOverlay(handles.extraktionVorschau,handles.wirbel1,double(seg < 0)*0.5);drawnow;
        
        einstellungIterationen = str2double(leseEinstellung(handles,'einstellungIterationen'));
        einstellungNU = str2double(leseEinstellung(handles,'einstellungNU'));
        einstellungMU = str2double(leseEinstellung(handles,'einstellungMU'));
        einstellungDT = str2double(leseEinstellung(handles,'einstellungDT'));
        
        oldSegSize = sum(sum(seg<0));
        seg = performSegmentation(einstellungIterationen,einstellungNU,einstellungMU,einstellungDT,seg,handles.wirbel1,scale);
        newSegSize = sum(sum(seg<0));
        
        %oldSegSize
        %newSegSize
        ( 1-str2double(leseEinstellung(handles,'einstellungTolerance') )*0.01 )*oldSegSize
        
        if (newSegSize >= ( 1-str2double(leseEinstellung(handles,'einstellungTolerance') )*0.01 )*oldSegSize)
            scale = scale * 2;
            
            if (scale > 1/str2double(leseEinstellung(handles,'einstellungScaleEnd')))
                showImageWithOverlay(handles.extraktionVorschau,handles.wirbel1,double(seg < 0)*0.5);drawnow;
                break;
            end
            
            seg = seg - 2;
        end
    else
        break;
    end
end
%profile off;profile viewer;

function phi = performSegmentation(numIteration,nu,mu,dt,phi,image,scale)
    dxc = 0.5 * [-1 0 1];
    dyc =dxc';
    epsilon = 10e-5;
    
    Dxpk = [0; -1; 1]; % D^x_+
    Dypk = Dxpk'; % D^x_+
    Dxmk = [-1; 1; 0]; % D^x_-
    Dymk = Dxmk'; % D^y_-
    
    Dxck = 0.5 * (Dxpk + Dxmk);
    Dyck = 0.5 * (Dypk + Dymk);
    
    h = waitbar(0,'','Name','Please wait...'); 
    
    imageSizeOriginal = size(image);
    
    %First resize to a maximum of 1000 pixels width
    scaleFirst = 1000 / imageSizeOriginal(2);
    phi =   imresize(phi,'OutputSize',  [imageSizeOriginal(1)*scaleFirst imageSizeOriginal(2)*scaleFirst]);
    image = imresize(image,'OutputSize',[imageSizeOriginal(1)*scaleFirst imageSizeOriginal(2)*scaleFirst]);
    
    imageSizeSmaller = size(image);
    
    phi =   imresize(phi,'OutputSize',  [imageSizeSmaller(1)*scale imageSizeSmaller(2)*scale]);
    image = imresize(image,'OutputSize',[imageSizeSmaller(1)*scale imageSizeSmaller(2)*scale]);
    
    phi(1,:) = 1;
    phi(end,:) = 1;
    phi(:,1) = 1;
    phi(:,end) = 1;
    
    [gradientAX,gradientAY]=gradient(image,0.1);
    edgemap = (gradientAX).^2 + (gradientAY).^2;
    edgemap = (1 ./ (1+edgemap)).^2;
    
    phi = double(bwdist(phi < 0) - bwdist(1-(phi <0)));

    for i=1:numIteration
        phix = imfilter(phi,dxc);
        phiy = imfilter(phi,dyc);
        
        gradPhi = (phix.^2 + phiy.^2 + epsilon).^0.5;
        
        change = dt .* edgemap .* gradPhi .* ( mu*(imfilter(phix ./ gradPhi ,dxc) + imfilter(phiy ./ gradPhi ,dyc)) + nu );
        %epsi = 0.1;
        %adeltaphi = 1/pi * epsi ./ (epsi.^2 + phi.^2);
        
        %Dxp = imfilter(phi,Dxpk,'same','replicate','conv');
        %Dyp = imfilter(phi,Dypk,'same','replicate','conv');
        
        %Dxc = imfilter(phi,Dxck,'same','replicate','conv');
        %Dyc = imfilter(phi,Dyck,'same','replicate','conv');
        
        %normGradPhi = sqrt(epsilon+Dxc.^2+Dyc.^2);
        
        %Dx = Dxp ./ normGradPhi;
        %Dx(isnan(Dx)) = 0;

        %Dy = Dyp ./ normGradPhi;
        %Dy(isnan(Dy)) = 0;
        
        %DxM = imfilter(Dxp./sqrt(Dxp.^2 + Dyc.^2 + epsilon),Dxmk,'same','replicate','conv');
        %DyM = imfilter(Dyp./sqrt(Dxc.^2 + Dyp.^2 + epsilon),Dymk,'same','replicate','conv');
        
        %change = dt .* edgemap .* adeltaphi .* normGradPhi .* ( (DxM + DyM)  * mu + nu);

        change(change<0) = 0;
        
        
%        figure(125);clf;imagesc(phi);
        phi = phi + change;
%         figure(126);clf;imagesc(phi);
%         figure(127);clf;imagesc(change);
%         pause

        if (mod(i,10) == 0)
            phi = double(bwdist(phi < 0) - bwdist(1-(phi <0)));
        end
        
        waitbar(i/numIteration,h,[num2str((i*100)/numIteration) '%']);
    end
    save
    
    %Delete isolated small regions
    CC = bwconncomp(phi < 0);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    
    for i=1:numel(numPixels)
        if (i ~= idx)
            phi(CC.PixelIdxList{i}) = 1;
        end
    end
    
    close(h);
    
    phi = imresize(phi,'OutputSize',[imageSizeOriginal(1) imageSizeOriginal(2)]);
    phi = double(bwdist(phi < 0) - bwdist(1-(phi <0)));

% --- Executes on button press in extraktionStop.
function extraktionStop_Callback(hObject, eventdata, handles)
set(handles.extraktionStop,'UserData',1);
guidata(hObject, handles);


% --------------------------------------------------------------------
function menueAllgemeinLaden_Callback(hObject, eventdata, handles)
[filename, pathname] = uigetfile('*.mat', 'Choose the GUI settings file to load');

if (~pathname)
    return;
end

%construct the path name of the file to be loaded
loadDataName = fullfile(pathname,filename);

load(loadDataName);

if (exist('pixelgroesse','var'))
    handles.pixelgroesse=pixelgroesse;
else
    handles.pixelgroesse = 0.25;
end

if (exist('wirbelHoeheStart','var'))
    handles.wirbelHoeheStart = wirbelHoeheStart;
    handles.wirbelHoeheEnde = wirbelHoeheEnde;
else
    handles.wirbelHoeheStart = 0;
    handles.wirbelHoeheEnde = 0;
end

if (exist('trabekelProPixel','var'))
    handles.trabekelProPixel = trabekelProPixel;
else
    handles.trabekelProPixel = 0;
end

handles.wirbel1 = wirbel1;
handles.wirbel2 = wirbel2;
handles.wirbel3 = wirbel3;
handles.wirbel1Seg = wirbel1Seg;
handles.wirbel2Seg = wirbel2Seg;
handles.trabekelListe = trabekelListe;
handles.umfang = umfang;
handles.volumen = volumen;
handles.wirbelHoehe = wirbelHoehe;
handles.wirbelBreiteOben = wirbelBreiteOben;
handles.wirbelBreiteUnten = wirbelBreiteUnten;
handles.calcifizierteFlaeche = calcifizierteFlaeche;
handles.cortikalisflaeche = cortikalisflaeche;
handles.cortikalisdicke = cortikalisdicke;
handles.trabekelflaeche = trabekelflaeche;
handles.trabekelDurchschnittsabstand = trabekelDurchschnittsabstand;
handles.trabekelDurchschnittdicke = trabekelDurchschnittdicke;

%Update GUI elements
set(handles.editPixelgroesse,'String',handles.pixelgroesse);


handles=zeigeErgebnisse(hObject, eventdata, handles);

guidata(hObject, handles);

% --------------------------------------------------------------------
function menueAllgemeinSpeichern_Callback(hObject, eventdata, handles)
[filename,pathname] = uiputfile('*.mat','Save your GUI settings');

if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
end
%construct the path name of the save location
saveDataName = fullfile(pathname,filename); 

%saves the gui data
%hgsave(saveDataName);

pixelgroesse = handles.pixelgroesse;

wirbel1 = handles.wirbel1;
wirbel2 = handles.wirbel2;
wirbel3 = handles.wirbel3;
wirbel1Seg = handles.wirbel1Seg;
wirbel2Seg = handles.wirbel2Seg;
trabekelListe = handles.trabekelListe;

umfang = handles.umfang;
volumen = handles.volumen;
wirbelHoehe = handles.wirbelHoehe;
wirbelHoeheStart = handles.wirbelHoeheStart;
wirbelHoeheEnde = handles.wirbelHoeheEnde;
wirbelBreiteOben = handles.wirbelBreiteOben;
wirbelBreiteUnten = handles.wirbelBreiteUnten;
calcifizierteFlaeche = handles.calcifizierteFlaeche;
cortikalisflaeche = handles.cortikalisflaeche;
cortikalisdicke = handles.cortikalisdicke;
trabekelflaeche = handles.trabekelflaeche;
trabekelListe = handles.trabekelListe;
trabekelDurchschnittsabstand = handles.trabekelDurchschnittsabstand;
trabekelDurchschnittdicke = handles.trabekelDurchschnittdicke;
trabekelProPixel = handles.trabekelProPixel;

einstellungExtraktionNU = handles.einstellungExtraktionNU;
einstellungExtraktionMU = handles.einstellungExtraktionMU;
einstellungExtraktionDT = handles.einstellungExtraktionDT;
einstellungExtraktionIterations = handles.einstellungExtraktionIterations;

save(saveDataName, ...
    'wirbel1',...
    'wirbel2',...
    'wirbel3',...
    'wirbel1Seg',...
    'wirbel2Seg',...
    'trabekelListe',...
    ...%Berechnete Daten
    'umfang',...
    'volumen',...
    'wirbelHoehe',...
    'wirbelHoeheStart',...
    'wirbelHoeheEnde',...
    'wirbelBreiteOben',...
    'wirbelBreiteUnten',...
    'calcifizierteFlaeche',...
    'cortikalisflaeche',...
    'cortikalisdicke',...
    'trabekelflaeche',...
    'trabekelDurchschnittsabstand',...
    'trabekelDurchschnittdicke',...
    'trabekelProPixel',...
    'pixelgroesse',...
    'einstellungExtraktionNU',...
    'einstellungExtraktionMU',...
    'einstellungExtraktionDT',...
    'einstellungExtraktionIterations'...
);

function handles=erzeugeErgebnisfeld(hObject, eventdata, handles,ergebnisText,ergebnisWert,einheit)
format long;
position = 1.02 - handles.ergebniszaehler * 0.055;

handles.ergebnisText{handles.ergebniszaehler} = uicontrol(handles.panelErgebnisse,'Style','text','Units','normalized','Position',[0.2,position,0.6,0.03],'String',ergebnisText);
handles.ergebnisDaten{handles.ergebniszaehler} = uicontrol(handles.panelErgebnisse,'Style','text','Units','normalized','Position',[0.2,position-0.025,0.6,0.03],'String',[num2str(ergebnisWert, '%14.2f'),einheit]);


handles.ergebniszaehler = handles.ergebniszaehler + 1;
format short;
%
function handles=loescheErgebnisfeld(hObject, eventdata, handles)
    for i=1:handles.ergebniszaehler-1
        if (handles.ergebnisDaten{i} ~= 0)
            delete(handles.ergebnisText{i});
            delete(handles.ergebnisDaten{i});
        end
    end
    
    handles.ergebniszaehler = 1;
%

function handles=zeigeErgebnisse(hObject, eventdata, handles)

handles = loescheErgebnisfeld(hObject, eventdata, handles);

if (handles.volumen ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Area:',handles.volumen*handles.pixelgroesse^2,'µm²');
end

if (handles.umfang ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Perimeter:',handles.umfang*handles.pixelgroesse,'µm');
end

if (handles.wirbelHoehe ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Height:',handles.wirbelHoehe*handles.pixelgroesse,'µm');
end

if (handles.wirbelBreiteOben ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Width top:',handles.wirbelBreiteOben*handles.pixelgroesse,'µm');
end

if (handles.wirbelBreiteUnten ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Width bottom:',handles.wirbelBreiteUnten*handles.pixelgroesse,'µm');
end

if (handles.calcifizierteFlaeche ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Calcified area:',handles.calcifizierteFlaeche*handles.pixelgroesse^2,'µm²');
end

if (handles.calcifizierteFlaeche ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Amount calcified area:',handles.calcifizierteFlaeche / handles.volumen,'');
end


if (handles.trabekelflaeche ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Trabecular area:',handles.trabekelflaeche*handles.pixelgroesse^2,'µm²');
end

if (handles.trabekelflaeche ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Amount trabecular area:',handles.trabekelflaeche/handles.volumen,'');
end

if (numel(handles.trabekelListe) > 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Number of trabeculae:',numel(handles.trabekelListe),'');
end

if (handles.trabekelProPixel > 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Trabeculae per mm:',1000*handles.trabekelProPixel/handles.pixelgroesse,'');
end

if (handles.trabekelDurchschnittsabstand > 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Avg. dist between trabeculae:',handles.trabekelDurchschnittsabstand*handles.pixelgroesse,'µm');
end

if (handles.trabekelDurchschnittdicke > 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Avg. trabecular width:',handles.trabekelDurchschnittdicke*handles.pixelgroesse,'µm');
end


if (handles.cortikalisflaeche ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Cortical area:',handles.cortikalisflaeche*handles.pixelgroesse^2,'µm²');
end

if (handles.cortikalisflaeche ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Amount cortical area:',handles.cortikalisflaeche/handles.volumen,'');
end

if (handles.cortikalisdicke ~= 0)
    handles = erzeugeErgebnisfeld(hObject, eventdata, handles,'Avg. cortical with:',handles.cortikalisdicke*handles.pixelgroesse,'µm');
end

guidata(hObject, handles);

function x = schwerpunkt( A )
%SCHWERPUNKT 
% Erwartet Matrix und berechnet ihren Schwerpunkt
%
    x = [0 0];
    anzahlMassepunkte = 0;
    for i=1:size(A,1)
        for j=1:size(A,2)
            if A(i,j) > 0
               x = x + [i j]*A(i,j);
               anzahlMassepunkte = anzahlMassepunkte + 1;
            end
        end
    end
    
    x = x / anzahlMassepunkte;
    
function [ lastPosition,image ] = goForward( start , direction, image,colorImage )

stop = 0;
position = start;
while (stop == 0)
    position = position + direction;
    if (round(position(1)) > 0 && round(position(1)) < size(image,1) && round(position(2)) > 0 && round(position(2)) < size(image,2) && image(round(position(1)),round(position(2))) > 0)
        %if (colorImage > 0)
        %    image(round(position(1)),round(position(2))) = colorImage;
        %end
    else
        stop = 1;
        lastPosition = position - direction;
    end
end



% --- Executes on button press in startVermessung.
function startVermessung_Callback(hObject, eventdata, handles)

if (handles.wirbel1Seg == 0)
    warndlg('Please do the extraction step first')
    return;
end

axes(handles.vermessungZeigeVermessung);
cla

waitbar_handle = waitbar(0,'Running measurements');

%Zunächst Umfang und Volumen berechnen
phi = handles.wirbel1Seg;
largerPhi = zeros(size(phi,1)+100,size(phi,2)+100);
largerPhi(50:size(phi,1)+49,50:size(phi,2)+49) = phi;

largerPhi = double(bwdist(largerPhi < 0) - bwdist(1-(largerPhi <0)));

epsilon = 3;deltaPhi = 1/pi * epsilon ./ (epsilon^2 + largerPhi.^2);

Dxck = [1 0 -1];
Dyck = [1;0;-1];

Dxc = imfilter(largerPhi,Dxck,'conv');
Dyc = imfilter(largerPhi,Dyck,'conv');

waitbar(0.1,waitbar_handle);

%Umfang: Integral über Betrag von Gradient des Levelsets 
%Volumen : Integral über phi < 0
handles.umfang = sum(sum(sqrt(Dxc.^2 + Dyc.^2) / 2 .* deltaPhi));
handles.volumen = sum(sum(phi < 0));


%Nun Höhe und Breite vermessen
imshow(handles.wirbel1,'Parent',handles.vermessungZeigeVermessung);

waitbar(0.2);
wirbelBinaer = (phi < 0);
schwerpunktWirbel = schwerpunkt(wirbelBinaer);
schwerpunktWirbelBild = double(wirbelBinaer);

%Bestimme höhe, indem Linien von Schwerpunkt ausgehen nach oben und unten
%gezeichnet werden
minWinkel = str2double(leseEinstellung(handles,'vermessungWinkelbereich1'));
maxWinkel = str2double(leseEinstellung(handles,'vermessungWinkelbereich2'));

minHoehe = Inf;

%sSchwerpunkt verschieben
radius = 10;

l=1;
for k=schwerpunktWirbel(1)-radius:1:schwerpunktWirbel(1)+radius
    waitbar(0.2+0.6*l/(2*radius),waitbar_handle);
    for j=schwerpunktWirbel(2)-radius:1:schwerpunktWirbel(2)+radius
        
        schwerpunktTmp = [k,j];
        %Verschiedene Winkel probieren
        for i=minWinkel:1:maxWinkel
            %Linien nach unten zeichnen
            pointBottom = goForward(schwerpunktTmp,[1 i/90],schwerpunktWirbelBild,0);
            %Linien nach oben zeichnen
            pointTop = goForward(schwerpunktTmp,-[1 i/90],schwerpunktWirbelBild,0);

            if (norm(pointTop-pointBottom) < minHoehe)
                lastPointBottom = pointBottom;
                lastPointTop = pointTop;
                minHoehe = norm(lastPointTop-lastPointBottom,2);
                minWinkel = i;
            end
        end
    end
    l=l+1;
end



%Hoehe einzeichnen nach oben
[lastPosUnten,schwerpunktWirbelBild ] = goForward( schwerpunktWirbel , [1 minWinkel/90], schwerpunktWirbelBild,2 );

%Hoehe einzeichnen nach unten
[lastPosOben,schwerpunktWirbelBild ] = goForward( schwerpunktWirbel ,- [1 minWinkel/90], schwerpunktWirbelBild,2 );

handles.wirbelHoeheStart = lastPosOben;
handles.wirbelHoeheEnde = lastPosUnten;
handles.wirbelHoehe = norm(lastPosOben-lastPosUnten);

hold on;
line([lastPosOben(2),lastPosUnten(2)],[lastPosOben(1),lastPosUnten(1)],'LineWidth',2,'Color',[1 0 0],'Parent',handles.vermessungZeigeVermessung);
waitbar(0.9);

%Obere Höhe feststellen

%Linie nach links zeichnen
winkelOben = [1 minWinkel/90];
if (minWinkel == 0)
    winkelLinks = [0 1];
else
    winkelLinks = [1 -90/minWinkel];
    winkelLinks = winkelLinks / norm(winkelLinks);
end
winkelRechts = -winkelLinks;


abstandLinieOben = str2double(leseEinstellung(handles,'vermessungAbstandOben')) * handles.wirbelHoehe;
[lastPosLinks,~ ] = goForward (lastPointTop+abstandLinieOben*winkelOben ,winkelLinks, schwerpunktWirbelBild,2 );
[lastPosRechts,~ ] = goForward (lastPointTop+abstandLinieOben*winkelOben ,winkelRechts, schwerpunktWirbelBild,2 );


handles.wirbelBreiteOben = norm(lastPosLinks-lastPosRechts);

hold on;
line([lastPosLinks(2),lastPosRechts(2)],[lastPosLinks(1),lastPosRechts(1)],'LineWidth',2,'Color',[1 0 0],'Parent',handles.vermessungZeigeVermessung);

%Wirbel Breite unten berechnen
abstandLinieUnten =   str2double(leseEinstellung(handles,'vermessungAbstandUnten')) * handles.wirbelHoehe;
[lastPosLinks,~ ] = goForward (lastPointBottom-abstandLinieUnten*winkelOben ,winkelLinks, schwerpunktWirbelBild,2 );
[lastPosRechts,~ ] = goForward (lastPointBottom-abstandLinieUnten*winkelOben ,winkelRechts, schwerpunktWirbelBild,2 );
waitbar(0.95);
handles.wirbelBreiteUnten = norm(lastPosLinks-lastPosRechts);

hold on;
line([lastPosLinks(2),lastPosRechts(2)],[lastPosLinks(1),lastPosRechts(1)],'LineWidth',2,'Color',[1 0 0],'Parent',handles.vermessungZeigeVermessung);

handles = zeigeErgebnisse(hObject, eventdata, handles);

delete( waitbar_handle );
clear waitbar_handle;

guidata(hObject, handles);

function einstellung = leseEinstellung(handles,wert)
    for i=1:handles.einstellungenZaehler-1
        if (ishandle(handles.einstellungenDaten{i}) && strcmp(get(handles.einstellungenDaten{i},'UserData'),wert)  == 1)
            einstellung = get(handles.einstellungenDaten{i},'String');
        end
    end


function handles=erzeugeEinstellung(hObject, eventdata, handles,einstellungDaten,histData)

sizePerField = 0.13;

position = (handles.einstellungenZaehlerVerschiebung-1) * sizePerField;

handles.einstellungenDaten{handles.einstellungenZaehler} = uipanel(handles.panelEinstellungen,'BorderType','none','Units','normalized','Position',[position,0,sizePerField,1],'Title','');
zaehlerParent = handles.einstellungenZaehler;

handles.einstellungenZaehler = handles.einstellungenZaehler + 1;

handles.einstellungenDaten{handles.einstellungenZaehler} = uicontrol(handles.einstellungenDaten{zaehlerParent},'Style','text','Units','normalized','Position',[0,0.75,1,0.25],'String',einstellungDaten{1});

handles.einstellungenZaehler = handles.einstellungenZaehler + 1;

handles.einstellungenDaten{handles.einstellungenZaehler} = uicontrol(handles.einstellungenDaten{zaehlerParent},'Style','edit','Units','normalized','Position',[0,0.5,1,0.25],'String',einstellungDaten{2},'UserData',einstellungDaten{3},'BackgroundColor','w');

if (nargin == 5)
    handles.einstellungenZaehler = handles.einstellungenZaehler + 1;

    handles.einstellungenDaten{handles.einstellungenZaehler} = axes('Parent',handles.einstellungenDaten{zaehlerParent},'Units','normalized','Position',[position,0,1,0.5]);
    plot(handles.einstellungenDaten{handles.einstellungenZaehler},histData{2},smooth(histData{1}))
end


handles.einstellungenZaehler = handles.einstellungenZaehler + 1;
handles.einstellungenZaehlerVerschiebung = handles.einstellungenZaehlerVerschiebung + 1;

guidata(hObject, handles);

%
function handles=loeschePanelEinstellungen(hObject, eventdata, handles)
    for i=1:handles.einstellungenZaehler-1
        if (ishandle(handles.einstellungenDaten{i}))
                delete(handles.einstellungenDaten{i});
            handles.einstellungenDaten{i} = 0;
        end
    end
    
    handles.einstellungenZaehler = 1;
    handles.einstellungenZaehlerVerschiebung = 1;
    
    guidata(hObject, handles);
%

% --- Executes on button press in knochendatenCalcifizierteFlaecheButton.
function knochendatenCalcifizierteFlaecheButton_Callback(hObject, eventdata, handles)

if (handles.wirbel1Seg == 0)
    warndlg('Please do the extraction steps first!')
    return;
end

% Calcifizierte Fläche
wirbel = handles.wirbel1;
wirbel(handles.wirbel1Seg > 0) = 1;

schwellwert = str2double(leseEinstellung(handles,'knochendatenSchwellwertCalc'));

handles.calcifizierteFlaeche = sum(sum(wirbel < schwellwert));

showImageWithOverlay(handles.vermessungZeigeVermessung,handles.wirbel1,double(wirbel < schwellwert)*0.6);

handles = zeigeErgebnisse(hObject, eventdata, handles);

guidata(hObject, handles);

% --- Executes on button press in knochendatenCortikalisflaecheButton.
function knochendatenCortikalisflaecheButton_Callback(hObject, eventdata, handles)

if (handles.wirbel2Seg == 0)
    warndlg('Please do the extraction steps first!')
    return;
end

if (handles.wirbel2 == 0)
    warndlg('Can´t find an image without endplate!')
    return;
end

schwellwert = str2double(leseEinstellung(handles,'knochendatenSchwellwertCort'));

segAussen = handles.wirbel1Seg;
segInnen = handles.wirbel2Seg;
bildOhneCorticalis = handles.wirbel2;

%check size
if (size(bildOhneCorticalis,1) ~= size(segInnen,1) || size(bildOhneCorticalis,2) ~= size(segInnen,2))
    errordlg('Image sizes do not match!') 
end

segInnerPhi = segInnen + str2double(leseEinstellung(handles,'knochendatenInnerSegOffset'));

cortikalisflaeche = segAussen;
cortikalisflaeche(segInnerPhi <= 0) = 0;
cortikalisflaeche(cortikalisflaeche>0) = 0;


cortikalisflaecheBild = bildOhneCorticalis;
cortikalisflaecheBild(cortikalisflaeche==0) = 1;
cortikalisflaecheBild(cortikalisflaecheBild < schwellwert) = -1;
cortikalisflaecheBild(cortikalisflaecheBild > 0) = 0;
cortikalisflaecheBild(cortikalisflaecheBild < 0) = 1;

handles.cortikalisflaeche = sum(sum(cortikalisflaecheBild > 0));


%Da Cortikalisflaeche bekannt ist, kann auch Dicke berechnet werden
phi = double(bwdist(cortikalisflaecheBild) - bwdist(1-(cortikalisflaecheBild)));

Dxck = [1 0 -1];
Dyck = [1;0;-1];

normalsU = imfilter(phi,Dxck,'conv') / 2;
normalsV = imfilter(phi,Dyck,'conv') / 2;


h = fspecial('gaussian', 5, 0.5);

normalsU = imfilter(normalsU,h);
normalsV = imfilter(normalsV,h);

betrag = sqrt(normalsU.^2 + normalsV.^2);

normalsU = normalsU ./ betrag;
normalsV = normalsV ./ betrag;

for i=1:size(phi,1)
   for j=1:size(phi,2)
       if (phi(i,j) >= 0 || phi(i,j) < -1)
          normalsU(i,j) = 0; 
          normalsV(i,j) = 0; 
       end
   end
end

image = cortikalisflaecheBild;

width = [];
for i=1:size(phi,1)
   for j=1:size(phi,2)
       if (abs(normalsU(i,j)) > 0 || abs(normalsV(i,j)) > 0 )
           [ lastPosition,image ] = goForward( [i j] ,-[normalsV(i,j) normalsU(i,j)], image,[] );
           width(end+1) = norm([i j]-lastPosition);
       end
   end
end

avg = sum(width) / numel(width);

handles.cortikalisdicke = avg;

showImageWithOverlay(handles.vermessungZeigeVermessung,handles.wirbel1,double(cortikalisflaecheBild)*0.6);

handles = zeigeErgebnisse(hObject, eventdata, handles);
guidata(hObject, handles);

function showImageWithOverlay(imageHandle,imageToShow,imageOverlay)
    axes(imageHandle);cla;
    imshow(imageToShow);
    hold on;
    h=imshow(cat(3, ones(size(imageToShow)), zeros(size(imageToShow)), zeros(size(imageToShow))));
    set(h,'AlphaData', imageOverlay,'AlphaDataMapping','none');
    

function trabekelTrabekelFinden_Callback(hObject, eventdata, handles)

    if (handles.wirbel2Seg == 0)
        warndlg('Please do the extraction steps first!')
        return;
    end

    schwellwertTrabekelgroesse = str2double(leseEinstellung(handles,'trabekelSchwellwertTrabekelgroesse'));
    schwellwertKnochen = str2double(leseEinstellung(handles,'trabekelSchwellwertKnochen'));
    trabekelVarianzGaussglaettung = str2double(leseEinstellung(handles,'trabekelVarianzGaussglaettung'));
    
    bildOhneDeckplatten = handles.wirbel2;
    
    trabekelbild = bildOhneDeckplatten;
    
    trabekelbild = imfilter(trabekelbild,fspecial('gaussian', 7, trabekelVarianzGaussglaettung));
    
    segInnerPhi = handles.wirbel2Seg + str2double(leseEinstellung(handles,'knochendatenInnerSegOffset'));
    
    trabekelbild(segInnerPhi >= 0) = 1;
    trabekelbild(trabekelbild > schwellwertKnochen) = 1;
    trabekelbild(trabekelbild <= schwellwertKnochen) = 0;
    
    trabekelbild = abs(1-trabekelbild);
    
    CC = bwconncomp(trabekelbild,8);
    %numPixels = cellfun(@numel,CC.PixelIdxList);

    %Filtere zu kleine Objekte
    i = 1;
    while i <= CC.NumObjects
        if (numel(CC.PixelIdxList{i}) < schwellwertTrabekelgroesse)
            CC.PixelIdxList(i) = [];
            CC.NumObjects = CC.NumObjects - 1;
        else
            i = i + 1;
        end
    end

    %Trabekelliste generieren
    handles.trabekelListe = {};
    for i=1:CC.NumObjects
       trab = struct('Name',['Trabecula ',num2str(i)],'Region',CC.PixelIdxList{i}); 
       handles.trabekelListe{i} = trab;
       handles.trabekelListe{i}.Neighbours = {};
    end
    
    if (numel(handles.trabekelListe) > 0)
        set(handles.trabekelListeBox,'Visible','On');
        set(handles.trabekelLoeschen,'Visible','On');
        set(handles.trabekelVerbinden,'Visible','On');
        
        %Liste der Trabekel fuellen
         listeTrabekelNamen = {};
         for i=1:numel(handles.trabekelListe)
             listeTrabekelNamen{end+1} = ['Trabecula ',num2str(i)];
         end
         set(handles.trabekelListeBox,'String',listeTrabekelNamen,'Value',1);
     
    end   
    
    
    handles=zeigeEinzeltrabekel(hObject, eventdata, handles,1,handles.trabekelZeigeBild,0,0);
    
    handles = zeigeErgebnisse(hObject, eventdata, handles);
    
    guidata(hObject, handles);
    

function handles=zeigeEinzeltrabekel(hObject, eventdata, handles,trabekelNummer,displayFigure,asd,zeigeNachbarn)
    axes(displayFigure);
    cla;
    segAussen = handles.wirbel1Seg;
    

    trabekelGesamtbild = zeros(size(handles.wirbel1));

    
    backgroundColoredImage = cat(3, zeros(size(handles.wirbel1)), zeros(size(handles.wirbel1)), zeros(size(handles.wirbel1)));
    for i=1:numel(handles.trabekelListe)
        [rows,cols] = ind2sub(size(handles.wirbel1),handles.trabekelListe{i}.Region);
        trabekelGesamtbild(handles.trabekelListe{i}.Region) = 1;
        
        
        if (sum(find(trabekelNummer==i)) > 0)
            for k=1:numel(rows)
                backgroundColoredImage(rows(k),cols(k),:) = [0 1 0];
            end
        else
            for k=1:numel(rows)
                backgroundColoredImage(rows(k),cols(k),:) = [1 0 0];
            end
        end
    end
    

    h2 = imshow(handles.wirbel1);
    hold on;
     
    if (exist('zeigeNachbarn') && zeigeNachbarn == 1 && numel(handles.trabekelListe{trabekelNummer}.Neighbours) > 0)
        %Alle Nachbarn in blau
        for i=1:numel(handles.trabekelListe{trabekelNummer}.Neighbours)
            [rows,cols] = ind2sub(size(handles.wirbel1),handles.trabekelListe{handles.trabekelListe{trabekelNummer}.Neighbours{i}}.Region);
            for k=1:numel(rows)
                backgroundColoredImage(rows(k),cols(k),:) = [0 0 1];
            end
        end
        
        %Aktuellen Nachbarn in hellblau
        trabekelNachbarNummer = get(handles.listboxNachbarn,'Value');
        [rows,cols] = ind2sub(size(handles.wirbel1),handles.trabekelListe{handles.trabekelListe{trabekelNummer}.Neighbours{trabekelNachbarNummer}}.Region);
        for k=1:numel(rows)
            backgroundColoredImage(rows(k),cols(k),:) = [0 1 1];
        end
    end
    
    h=imshow(backgroundColoredImage,'Parent',displayFigure);
    set(h,'AlphaData', 1*double(trabekelGesamtbild));
      
    
    hold off;
    
    set(h,'ButtonDownFcn',{@dynamicViewerFunction,displayFigure}); 
    set(h2,'ButtonDownFcn',{@dynamicViewerFunction,displayFigure}); 
    
    guidata(hObject, handles);


% --- Executes on button press in trabekelLoeschen.
function trabekelLoeschen_Callback(hObject, eventdata, handles)

    if (numel(handles.trabekelListe) == 0)
        warndlg('Can´t find any trabeculae!')
        return;
    end
    
    trabekelNummer = get(handles.trabekelListeBox,'Value');
    
    for i=1:numel(trabekelNummer)
        handles.trabekelListe(trabekelNummer(i)-(i-1)) = [];
    end
    
    %Liste der Trabekel fuellen
    listeTrabekelNamen = {};
    for i=1:numel(handles.trabekelListe)
        listeTrabekelNamen{end+1} = ['Trabecula ',num2str(i)];
    end
    if (trabekelNummer > numel(listeTrabekelNamen) )
        set(handles.trabekelListeBox,'String',listeTrabekelNamen,'Value',trabekelNummer-1);
        handles = zeigeEinzeltrabekel(hObject, eventdata, handles,trabekelNummer-1,handles.trabekelZeigeBild);
    else
        set(handles.trabekelListeBox,'String',listeTrabekelNamen,'Value',trabekelNummer);
        handles = zeigeEinzeltrabekel(hObject, eventdata, handles,trabekelNummer,handles.trabekelZeigeBild);
    end
    
    uicontrol(handles.trabekelListeBox);
    handles = zeigeErgebnisse(hObject, eventdata, handles);
    
    guidata(hObject, handles);


% --- Executes on slider movement.
function handles=einzeltrabekelSlider_Callback(hObject, eventdata, handles)
    if (numel(handles.trabekelListe) == 0)
        warndlg('Can´t find any trabeculae!')
        return;
    end
    
    trabNummer = round(get(hObject,'Value'));
    
    %Liste der Nachbarn fuellen
    listeTrabekelNamen = {};
    for i=1:numel(handles.trabekelListe{trabNummer}.Neighbours)
        listeTrabekelNamen{end+1} = ['Trabecula ',num2str(handles.trabekelListe{trabNummer}.Neighbours{i})];
    end
    set(handles.listboxNachbarn,'String',listeTrabekelNamen,'Value',1);

    handles=zeigeEinzeltrabekel(hObject, eventdata, handles,trabNummer,handles.einzeltrabekelZeige,[],1);
    
    guidata(hObject, handles);

function einzeltrabekelLoescheNachbarn_Callback(hObject, eventdata, handles)
    if (numel(handles.trabekelListe) == 0)
        warndlg('Can´t find any trabeculae!')
        return;
    end
    
    trabNummer = get(handles.einzeltrabekelSlider,'Value');
    trabekelNachbarNummer = get(handles.listboxNachbarn,'Value');
    
    neighbours=handles.trabekelListe{trabNummer}.Neighbours;
    if numel(neighbours) > 0
        neighbours(trabekelNachbarNummer) = [];
        handles.trabekelListe{trabNummer}.Neighbours = neighbours;
    end
    handles=einzeltrabekelSlider_Callback(handles.einzeltrabekelSlider, eventdata, handles);
    guidata(hObject, handles);
 
function trabekelanalyseStart_Callback(hObject, eventdata, handles)
    
    waitbar_handle = waitbar(0,'Analyzing trabeculae');
    %Trabekel pro MM hoehe
    
    %Erzeuge Trabekelbild
    trabekelProMMBild = zeros(size(handles.wirbel1));
    for i=1:numel(handles.trabekelListe)
        trabekelProMMBild(handles.trabekelListe{i}.Region) = 1;
    end
    
    if (~handles.wirbel3)
        warndlg('Can´t find an image without endplate and transvese processes')
        return
    end
    trabekelProMMBild(handles.wirbel3 == 1) = 0;
    
    sum1 = 0;
    for i=floor(handles.wirbelHoeheStart(1)):ceil(handles.wirbelHoeheEnde(1))
        sum1 = sum1 + numel(unique(trabekelProMMBild(i,:))) - 1;
    end
    mean = sum1 / (ceil(handles.wirbelHoeheEnde(1)) - floor(handles.wirbelHoeheStart(1)));
    handles.trabekelProPixel = mean;
    


    dist = [];
    width = [];
    for i=1:numel(handles.trabekelListe)
        waitbar(i / numel(handles.trabekelListe))
        %Durchschnittsabstand zu Nachbarn berechnen
        trabekelEinzelbild = zeros(size(handles.wirbel1));
        trabekelEinzelbild(handles.trabekelListe{i}.Region) = 1;

        trabekelSchwerpunkt = schwerpunkt(trabekelEinzelbild);

        distNachbar = [];
        for j=1:numel(handles.trabekelListe{i}.Neighbours)
            trabekelEinzelbildNachbarn = zeros(size(handles.wirbel1));
            trabekelEinzelbildNachbarn(handles.trabekelListe{handles.trabekelListe{i}.Neighbours{j}}.Region) = 1;

            nachbarSchwerpunkt = schwerpunkt(trabekelEinzelbildNachbarn);
            distNachbar(end+1) = norm(trabekelSchwerpunkt-nachbarSchwerpunkt);
        end
        if (numel(distNachbar) > 0)
            dist(end+1) = sum(distNachbar) / numel(distNachbar);
        end
        
        %Trabekeldicke berechnen
        phi = double(bwdist(trabekelEinzelbild) - bwdist(1-(trabekelEinzelbild)));

        Dxck = [1 0 -1];
        Dyck = [1;0;-1];

        normalsU = imfilter(phi,Dxck,'conv') / 2;
        normalsV = imfilter(phi,Dyck,'conv') / 2;


        h = fspecial('gaussian', 5, 0.5);

        normalsU = imfilter(normalsU,h);
        normalsV = imfilter(normalsV,h);

        betrag = sqrt(normalsU.^2 + normalsV.^2);

        normalsU = normalsU ./ betrag;
        normalsV = normalsV ./ betrag;

        for m=1:size(phi,1)
           for n=1:size(phi,2)
               if (phi(m,n) >= 0 || phi(m,n) < -1)
                  normalsU(m,n) = 0; 
                  normalsV(m,n) = 0; 
               end
           end
        end

        image = trabekelEinzelbild;
        image = abs(image-1);
        tmpWidth = [];
        for m=1:size(phi,1)
           for n=1:size(phi,2)
               if (abs(normalsU(m,n)) > 0 || abs(normalsV(m,n)) > 0 )
                   [ lastPosition,image ] = goForward( [m n] ,-[normalsV(m,n) normalsU(m,n)], image,[] );
                   %line([n,lastPosition(2)],[m,lastPosition(1)],'LineWidth',2,'Color',[1 0 0]);
                   tmpWidth(end+1) = norm([m n]-lastPosition);
               end
           end
        end
        width(end+1) = sum(tmpWidth) / numel(tmpWidth);
    end
    handles.trabekelDurchschnittsabstand = sum(dist) / numel(dist);
    handles.trabekelDurchschnittdicke = sum(width) / numel(width);
    
    %Trabekelfläche:
    trabekelBild = zeros(size(handles.wirbel1));
    for i=1:numel(handles.trabekelListe)
        trabekelBild(handles.trabekelListe{i}.Region) = 1;
    end
    handles.trabekelflaeche = sum(sum(trabekelBild > 0));
    
    handles = zeigeErgebnisse(hObject, eventdata, handles);
    
    delete( waitbar_handle );
    clear waitbar_handle;
    
    guidata(hObject, handles);

function trabekelWeiterZuEinzeltrabekel_Callback(hObject, eventdata, handles)
schritteEinzeltrabekel_Callback(hObject, eventdata, handles);

function einzeltrabekelWeiterAuswertung_Callback(hObject, eventdata, handles)
schritteTrabekelanalyse_Callback(hObject, eventdata, handles);

function einzeltrabekelZurueck_Callback(hObject, eventdata, handles)
schritteTrabekel_Callback(hObject, eventdata, handles);

function trabekelanalyseZurueck_Callback(hObject, eventdata, handles)
schritteEinzeltrabekel_Callback(hObject, eventdata, handles);


function handles=trabekelListeBox_Callback(hObject, eventdata, handles)
handles = zeigeEinzeltrabekel(hObject, eventdata, handles,get(handles.trabekelListeBox,'Value'),handles.trabekelZeigeBild,[]);

guidata(hObject, handles);

function trabekelVerbinden_Callback(hObject, eventdata, handles)

    if (numel(handles.trabekelListe) == 0)
        errordlg('Can´t find any trabeculae!')
        return;
    end
    
    trabekelNummer = sort(get(handles.trabekelListeBox,'Value'));
    
    if (numel(trabekelNummer) < 2)
        return;
    end
    
    reg = [];
    for i=1:numel(trabekelNummer)
        reg = [reg;handles.trabekelListe{trabekelNummer(i)}.Region];
    end
    
    handles.trabekelListe{trabekelNummer(1)}.Region = reg;

    for i=2:numel(trabekelNummer)
        handles.trabekelListe(trabekelNummer(i)-(i-2)) = [];
    end
    
    %Liste der Trabekel fuellen
    listeTrabekelNamen = {};
    for i=1:numel(handles.trabekelListe)
        listeTrabekelNamen{end+1} = ['Trabecula ',num2str(i)];
    end
    set(handles.trabekelListeBox,'String',listeTrabekelNamen,'Value',min(trabekelNummer));
    
    handles = zeigeEinzeltrabekel(hObject, eventdata, handles,trabekelNummer(1),handles.trabekelZeigeBild);
    handles = zeigeErgebnisse(hObject, eventdata, handles);
    
    guidata(hObject, handles);

function editPixelgroesse_Callback(hObject, eventdata, handles)
    handles.pixelgroesse = str2double(get(hObject,'String'));
    handles = zeigeErgebnisse(hObject, eventdata, handles);
    guidata(hObject, handles);

function menueInfo_Callback(hObject, eventdata, handles)
msgbox({'histoGUI version 1.1 beta (11.07.2012)','','Created by:','','Hendrik Dirks','Institute for Computational and Applied Mathematics','Faculty of Mathematics and Computational Science','Westfälische Wilhelms-Universität Münster','','Email: hendrik.dirks@math.uni-muenster.de'},'About') 

function menueAllgemeinExport_Callback(hObject, eventdata, handles)
    
[filename,pathname] = uiputfile('*.xls','Export your results');

if pathname == 0 %if the user pressed cancelled, then we exit this callback
    return
end
%construct the path name of the save location
saveDataName = fullfile(pathname,filename); 

d = {...
    'Area:',...
    'Perimeter:',...
    'Height:',...
    'Width top:',...
    'Width bottom:',...
    'Calcified area:',...
    'Amount calcified area:',...
    'Trabecular area:',...
    'Amount trabecular area:',...
    'Number of trabeculae:',...
    'Trabeculae per mm:',...
    'Avg. dist between trabeculae:',...
    'Avg. trabecular width:',...
    'Cortical area:',...
    'Amount cortical area:',...
    'Avg. cortical with:'...
    ;...
    [num2str(handles.volumen*handles.pixelgroesse^2),'µm²'],...
    [num2str(handles.umfang*handles.pixelgroesse),'µm'],...
    [num2str(handles.wirbelHoehe*handles.pixelgroesse),'µm'],...
    [num2str(handles.wirbelBreiteOben*handles.pixelgroesse),'µm'],...
    [num2str(handles.wirbelBreiteUnten*handles.pixelgroesse),'µm'],...
    [num2str(handles.calcifizierteFlaeche*handles.pixelgroesse^2),'µm²'],...
    num2str(handles.calcifizierteFlaeche / handles.volumen),...
    [num2str(handles.trabekelflaeche*handles.pixelgroesse^2),'µm²'],...
    num2str(handles.trabekelflaeche/handles.volumen),...
    num2str(numel(handles.trabekelListe)),...
    num2str(handles.trabekelProPixel/handles.pixelgroesse),...
    [num2str(handles.trabekelDurchschnittsabstand*handles.pixelgroesse),'µm'] ...
    [num2str(handles.trabekelDurchschnittdicke*handles.pixelgroesse),'µm'] ...
    [num2str(handles.cortikalisflaeche*handles.pixelgroesse^2),'µm²'] ...
    num2str(handles.cortikalisflaeche/handles.volumen) ...
    [num2str(handles.cortikalisdicke*handles.pixelgroesse),'µm']...
    ;};

xlswrite(saveDataName, d);

function dynamicViewerFunction(hObject, eventdata,displayFigure)
    handles = guidata(hObject);
    trabekelZeigeBild_ButtonDownFcn(hObject, eventdata, handles,displayFigure);


% --- Executes on mouse press over axes background.
function trabekelZeigeBild_ButtonDownFcn(hObject, eventdata, handles,displayFigure)
% hObject    handle to trabekelZeigeBild (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pt = get(displayFigure, 'CurrentPoint');

trabekelGesamtbild = zeros(size(handles.wirbel1));
%Identify Trabekel
for i=1:numel(handles.trabekelListe)
    trabekelGesamtbild(handles.trabekelListe{i}.Region) = i;
end

trabekelNumber = trabekelGesamtbild(round(pt(1,2)),round(pt(1,1)));

if (handles.trabekelStep == 1)

    if (trabekelNumber > 0)
        oldValue = get(handles.trabekelListeBox,'Value');
        set(handles.trabekelListeBox,'Value',[oldValue,trabekelNumber]);
        handles = zeigeEinzeltrabekel(handles.trabekelListeBox, eventdata, handles,get(handles.trabekelListeBox,'Value'),handles.trabekelZeigeBild,[]);
    end
elseif (handles.trabekelStep == 2)
    if (trabekelNumber > 0)
        
        %Add neighbour to this trabecula
        neighbours = handles.trabekelListe{get(handles.einzeltrabekelSlider,'Value')}.Neighbours;
        
        noInsert = 0;
        for i=1:numel(neighbours)
            if neighbours{i} == trabekelNumber
                noInsert = 1;
            end
        end
        
        if (noInsert == 0)
            neighbours{end+1} = trabekelNumber;
        end
        
        handles.trabekelListe{get(handles.einzeltrabekelSlider,'Value')}.Neighbours = neighbours;
        
        if (noInsert == 0)
            %Add this trabecular to neighbour
            neighbours = handles.trabekelListe{trabekelNumber}.Neighbours;

            for i=1:numel(neighbours)
                if neighbours{i} == get(handles.einzeltrabekelSlider,'Value')
                    noInsert = 1;
                end
            end

            if (noInsert == 0)
                neighbours{end+1} = get(handles.einzeltrabekelSlider,'Value');
            end

            handles.trabekelListe{trabekelNumber}.Neighbours = neighbours;
        end
        handles=einzeltrabekelSlider_Callback(handles.einzeltrabekelSlider, eventdata, handles);

    end
end

guidata(handles.trabekelListeBox, handles);


% --- Executes on selection change in listboxNachbarn.
function handles=listboxNachbarn_Callback(hObject, eventdata, handles)

trabNummer = round(get(handles.einzeltrabekelSlider,'Value'));

handles=zeigeEinzeltrabekel(hObject, eventdata, handles,trabNummer,handles.einzeltrabekelZeige,[],1);

guidata(hObject, handles);


function einzeltrabekelSlider_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function bildTyp_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editPixelgroesse_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function trabekelListeBox_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Max',2);

function listboxNachbarn_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
