function FLIMZ(varargin)
%% FLIMZ.m USAGE NOTES
%{

Syntax
-----------------------------------------------------
    FLIMZ()
    FLIMZ(datfilefdir)


Description
-----------------------------------------------------
    FLIMZ() can be run with no arguments passed in. In this case user
    will be prompted to select a directory which contains the FLIM dat 
    file along with the corresponding CCD images. Optionally this function can 
    be called using FLIMZ(datfilefdir) where the full path to the data directory
    is explicitly provided.
    

Useage Definitions
-----------------------------------------------------


    FLIMZ()
        launches a GUI that will first ask whether you want to compile
        a dataset output from Bh SPC-Image. Specifically it requires
        that...
            - Color Coded VAlue
            - Chi
            - Pixel Intensities
            - Color Coded Image
        ...are exported from the FLIM analysis software. This GUI will also
        ask the user if it wants to load one of these compiled .dat files.
        If this 'Load data file' option is clicked, the user is prompted to
        select a .dat file. After this the main FLIMCCD analysis GUI is
        launched.
 


Example
-----------------------------------------------------

% Create 2D triangulated mesh
    XY = randn(10,2);
    TR2D = delaunayTriangulation(XY);
    vrts = TR2D.Points;
    tets = TR2D.ConnectivityList;

    xmlmesh(vrts,tets,'xmlmesh_2D.xml')



See Also
-----------------------------------------------------
http://bradleymonk.com/xmlmesh
http://fenicsproject.org
>> web(fullfile(docroot, 'matlab/math/triangulation-representations.html'))


Attribution
-----------------------------------------------------
% Created by: Bradley Monk
% email: brad.monk@gmail.com
% website: bradleymonk.com
% 2016.04.19

%}


%% ESTABLISH STARTING PATHS


clc; close all; clear all;
disp('clearing matlab workspace');

thisfile = 'FLIMZ.m';
thisfilepath = fileparts(which(thisfile));
cd(thisfilepath);


global datfilefdir ccdfilegdir ccdfilerdir
global datfilef ccdfileg ccdfiler 
datfilefdir = '';
ccdfilegdir = '';
ccdfilerdir = '';
datfilef    = '';
ccdfileg    = '';
ccdfiler    = '';


%% MANUALLY SET PER-SESSION PATH PARAMETERS IF WANTED

global datadir datafile datadate
datadir = '';
datafile = '';
datadate = '';

global imgpath
imgpath = '';


%% ESTABLISH GLOBALS AND SET STARTING VALUES


global ROI

ROI.POS  = {};
ROI.PH   = {};
ROI.POSd = {};
ROI.PHd  = {};
ROI.TYPE  = {};


global haxROIS haxTABH haxAPI hREC hLINE
global LifeImageFile FLIMcmap gROI
global intenseThreshMIN intenseThreshMAX intenseThreshPMIN intenseThreshPMAX
global lifeThreshMIN lifeThreshMAX chiThreshMIN chiThreshMAX magnification maglevel
global flimdata flimdat flimtable flimd ROInames Datafilename 
global hROI hROIs ROImask ROIpos ROIarea dendritesize dpos
global ChiGood IntensityGood LifeGood AllGood
global ROI_LIFETIME ROI_INTENSITY ROI_CHI
global ROI_LIFETIME_MEAN ROI_INTENSITY_MEAN ROI_CHI_MEAN
global ROI_imgG ROI_imgR ROI_imgG_MEAN ROI_imgR_MEAN
global imXlim imYlim VxD dVOL
global imgG imgR imgF haxes haxnum stampSize
global phFLIM phCCDR phCCDG
global sROI sROIpos sROIarea sROImask flimdats
global tempV1 tempV2 tempV3 tempV4 tempV5 tempV6 tempV7 tempV8
global ROIcsv CLIMslider APIgroup ROIgroup TRACEgroup
global radioh doLineTrace pairedRadio singleRadio delROImodeOFF delROImodeON

tempV1 = [];
tempV2 = [];
tempV3 = [];
tempV4 = [];
tempV5 = [];
tempV6 = [];
tempV7 = [];
tempV8 = [];


flimdats = {};
sROI = [];
sROIpos = [];
sROIarea = [];
sROImask = [];
LifeImageFile = 0;
FLIMcmap = FLIMcolormap;
intenseThreshMIN = 85.000;
intenseThreshMAX = 99.999;
intenseThreshPMIN = 2;
intenseThreshPMAX = 10;
lifeThreshMIN = 500;
lifeThreshMAX = 2900;
chiThreshMIN = 0.7;
chiThreshMAX = 2.0;
magnification = 6;
maglevel = 6;
dendritesize = maglevel*5;
dpos = [];
flimdata = {};
flimdat = [];
flimtab = [];
flimd = [];
ROInames = '';
Datafilename = '';
hROI = [];
ROImask = [];
ROIpos = [];
ROIarea = [];
ChiGood = [];
IntensityGood = [];
LifeGood = [];
AllGood = [];
ROI_LIFETIME = [];
ROI_INTENSITY = [];
ROI_CHI = [];
ROI_LIFETIME_MEAN = [];
ROI_INTENSITY_MEAN = [];
ROI_CHI_MEAN = [];
ROI_imgG = [];
ROI_imgR = [];
ROI_imgG_MEAN = [];
ROI_imgR_MEAN = [];
VxD = 1;
dVOL = 1;
imgG = [];
imgR = [];
imgF = [];
haxnum = 1:3;
stampSize = 11;
hROIs = {};

global boxtype
boxtype = 'freehand'; % freehand:1  rectangle:2  elipse:3




%----------------------------------------------------
%           INITIALIZE CONTAINERS
%----------------------------------------------------


datastack = zeros(1,1,3,'double');
lifetime = zeros(1, 1);
intensity = zeros(1, 1);
chi = zeros(1, 1);
lifetimeimage = zeros(1, 1);
intensityimage = zeros(1, 1);
xdim = 0;
ydim = 0;
saveROI = zeros(200, 17);
saveData = zeros(200, 9);





%########################################################################
%%              MAIN ANALYSIS GUI WINDOW SETUP 
%########################################################################

MAINGUIFIG = figure('Units', 'normalized','OuterPosition', [.02 .05 .85 .87], 'BusyAction',...
    'cancel', 'Name', 'Lifetime image', 'Tag', 'lifetime image','Visible', 'Off', ...
    'KeyPressFcn', {@keypresszoom,1},'Color',[.99 .99 .99],'Resize','off','MenuBar','none'); % 


haxCCDG = axes('Parent', MAINGUIFIG, 'NextPlot', 'Add',...
    'Position', [0.05 0.14 0.8 0.8], 'OuterPosition', [-.2 0 1 1],...
    'PlotBoxAspectRatio', [1 1 1],'XColor','none','YColor','none'); 
    % ,'XDir','reverse',...

haxCCDR = axes('Parent', MAINGUIFIG, 'NextPlot', 'Add',...
    'Position', [0.05 0.14 0.8 0.8], 'OuterPosition', [-.2 0 1 1],...
    'PlotBoxAspectRatio', [1 1 1],'XColor','none','YColor','none'); 
    % ,'XDir','reverse',...

haxFLIM = axes('Parent', MAINGUIFIG, 'NextPlot', 'Add',...
    'Position', [0.05 0.14 0.8 0.8], 'OuterPosition', [-.2 0 1 1],...
    'PlotBoxAspectRatio', [1 1 1],'XColor','none','YColor','none');

haxAPI = axes('Parent', MAINGUIFIG, 'NextPlot', 'Add','Color','none',...
    'Position', [0.05 0.14 0.8 0.8], 'OuterPosition', [-.2 0 1 1],...
    'PlotBoxAspectRatio', [1 1 1],'XColor','none','YColor','none');

axes(haxAPI)
TRACEgroup = hggroup;
APIgroup = hggroup;

haxROIS = axes('Parent', MAINGUIFIG, 'NextPlot', 'Add','Color','none',...
    'Position', [0.05 0.15 0.8 0.8], 'OuterPosition', [-.2 0 1 1],...
    'PlotBoxAspectRatio', [1 1 1],'XColor','none','YColor','none');

axes(haxROIS)
ROIgroup = hggroup;

linkaxes([haxFLIM,haxCCDG,haxCCDR,haxAPI,haxROIS],'xy')
haxes = {haxFLIM, haxCCDG, haxCCDR, haxAPI, haxROIS};


CLIMsliderTxt = uicontrol('Parent', MAINGUIFIG, 'Style', 'Text', 'Units', 'normalized', ...
    'Position', [.01 .97 .06 .02], 'FontSize', 10, 'String', 'FLIM Level','BackgroundColor',[1 1 1]);
CLIMslider = uicontrol('Parent', MAINGUIFIG, 'Units', 'normalized','Style','slider',...
	'Max',150,'Min',1,'Value',50,'SliderStep',[0.01 0.10],...
	'Position', [.07 .97 .23 .02], 'Callback', @setCLimFLIM);

CLIMsliderCCDGTxt = uicontrol('Parent', MAINGUIFIG, 'Style', 'Text', 'Units', 'normalized', ...
    'Position', [.01 .95 .06 .02], 'FontSize', 10, 'String', 'GFP Level','BackgroundColor',[1 1 1]);
CLIMsliderCCDG = uicontrol('Parent', MAINGUIFIG, 'Units', 'normalized','Style','slider',...
	'Max',150,'Min',1,'Value',50,'SliderStep',[0.01 0.10],...
	'Position', [.07 .95 .23 .02], 'Callback', @setCLimGFP);

CLIMsliderCCDRTxt = uicontrol('Parent', MAINGUIFIG, 'Style', 'Text', 'Units', 'normalized', ...
    'Position', [.01 .93 .06 .02], 'FontSize', 10, 'String', 'RFP Level','BackgroundColor',[1 1 1]);
CLIMsliderCCDR = uicontrol('Parent', MAINGUIFIG, 'Units', 'normalized','Style','slider',...
	'Max',150,'Min',1,'Value',50,'SliderStep',[0.01 0.10],...
	'Position', [.07 .93 .23 .02], 'Callback', @setCLimRFP);


CONsliderTxt = uicontrol('Parent', MAINGUIFIG, 'Style', 'Text', 'Units', 'normalized', ...
    'Position', [.30 .97 .06 .02], 'FontSize', 10, 'String', 'F Contrast','BackgroundColor',[1 1 1]);
CONslider = uicontrol('Parent', MAINGUIFIG, 'Units', 'normalized','Style','slider',...
	'Max',100,'Min',50,'Value',80,'SliderStep',[0.01 0.10],...
	'Position', [.36 .97 .21 .02], 'Callback', @setCLimFLIM);

CONsliderCCDGTxt = uicontrol('Parent', MAINGUIFIG, 'Style', 'Text', 'Units', 'normalized', ...
    'Position', [.30 .95 .06 .02], 'FontSize', 10, 'String', 'G Contrast','BackgroundColor',[1 1 1]);
CONsliderCCDG = uicontrol('Parent', MAINGUIFIG, 'Units', 'normalized','Style','slider',...
	'Max',100,'Min',50,'Value',80,'SliderStep',[0.01 0.10],...
	'Position', [.36 .95 .21 .02], 'Callback', @setCLimGFP);

CONsliderCCDRTxt = uicontrol('Parent', MAINGUIFIG, 'Style', 'Text', 'Units', 'normalized', ...
    'Position', [.30 .93 .06 .02], 'FontSize', 10, 'String', 'R Contrast','BackgroundColor',[1 1 1]);
CONsliderCCDR = uicontrol('Parent', MAINGUIFIG, 'Units', 'normalized','Style','slider',...
	'Max',100,'Min',50,'Value',80,'SliderStep',[0.01 0.10],...
	'Position', [.36 .93 .21 .02], 'Callback', @setCLimRFP);




%----------------------------------------------------
%           SWAP IMAGE PANNEL
%----------------------------------------------------
ShowIMGpanelH = uibuttongroup('Parent', MAINGUIFIG, 'Visible','on',...
                  'Units', 'normalized','BackgroundColor',[1 1 1],'FontSize', 11,...
                  'Position',[.001 .60 .06 .20],'Title','Show IMG',...
                  'SelectionChangedFcn',@imTOG,'Tag','imgTOG');

showFIMh = uicontrol('Parent', ShowIMGpanelH, 'Units', 'normalized','Style','togglebutton', ...
    'Position', [.001 .67 .99 .32], 'String', 'FLIM', 'FontSize', 11);

showGIMh = uicontrol('Parent', ShowIMGpanelH, 'Units', 'normalized','Style','togglebutton', ...
    'Position', [.001 .34 .99 .32], 'String', 'GFP', 'FontSize', 11);

showRIMh = uicontrol('Parent', ShowIMGpanelH, 'Units', 'normalized','Style','togglebutton', ...
    'Position', [.001 .001 .99 .32], 'String', 'RFP', 'FontSize', 11);



%----------------------------------------------------
%           CREATE CUSTOM TOOLBAR
%----------------------------------------------------

MAINGUItoolbar = uitoolbar(MAINGUIFIG);


UITB_chgIM = uipushtool(MAINGUItoolbar,'TooltipString','Reset Workspace',...
               'ClickedCallback',@resetwsFun,'Tag','restwsFUN');


UITB_delROI = uitoggletool(MAINGUItoolbar,'TooltipString','Delete ROI Mode',...
               'ClickedCallback',@delROImode,'Tag','delModeTB');


% SET TOOLBAR ICONS
% fullfile(matlabroot,'toolbox','matlab','icons')
[img,map] = imread(fullfile(matlabroot,'toolbox','matlab','icons','greenarrowicon.gif'));             
UITB_chgIM.CData = ind2rgb(img,map);

[img,map] = imread(fullfile(matlabroot,'toolbox','matlab','icons','tool_data_brush.png'));
UITB_delROI.CData = im2double(img);

%---------------------
delROITOGh = uibuttongroup('Parent', MAINGUIFIG, 'Visible','on',...
                  'Units', 'normalized','BackgroundColor',[1 1 1],...
                  'Position',[.30 .02 .15 .06],'Title','DELETE ROI MODE',...
                  'SelectionChangedFcn',@delROImode,'Tag','delModeTOG');


delROImodeOFF  = uicontrol('Parent', delROITOGh,'Style','togglebutton', 'Units', 'normalized', ...
    'Position', [.01 .01 .45 .98], 'String', 'OFF', 'FontSize', 11);

delROImodeON = uicontrol('Parent', delROITOGh,'Style','togglebutton', 'Units', 'normalized', ...
    'Position', [.51 .01 .45 .98], 'String', 'ON', 'FontSize', 11);










%----------------------------------------------------
%           ROI ACQUISITION PANEL
%----------------------------------------------------
IPpanelH = uipanel('Title','Image Processing','FontSize',10,...
    'BackgroundColor',[1 1 1],...
    'Position', [0.62 0.71 0.37 0.28]); % 'Visible', 'Off',
%--------

getROIh = uicontrol('Parent', IPpanelH, 'Units', 'normalized', ...
    'Position', [.02 .78 .42 .18], 'FontSize', 12, 'String', 'MEASURE ROI',...
    'Callback', @getROI); 


ROICounth = uicontrol('Parent', IPpanelH, 'Units', 'normalized', ...
    'Position', [.50 .87 .20 .12], 'FontSize', 10, 'String', 'ROI Count',...
    'Callback', @recountROIs);

ROI_IDh = uicontrol('Parent', IPpanelH, 'Style', 'Text', 'Units', 'normalized', ...
    'Position', [.50 .80 .20 .10], 'FontSize', 14, 'String', '0','BackgroundColor',[1 .95 .95]);
% ROI_IDh.String = int2str(0);


stampSizeTxt = uicontrol('Parent', IPpanelH, 'Style', 'Text', 'Units', 'normalized', ...
    'Position', [.78 .87 .20 .11], 'FontSize', 10, 'String', 'ROI Size','BackgroundColor',[1 1 1]);


stampvals = {'9','10','11','12','14','16','18','20'};
stampSizeH = uicontrol('Parent', IPpanelH,'Style', 'popup',...
    'Units', 'normalized', 'String',stampvals,...
    'Position', [.78 .78 .20 .12],'BackgroundColor',[.9 .9 .9],...
    'Callback', @getStampSize, 'Enable','on');
stampSizeH.Value = 3;
stampSizeH.ListboxTop = 3;


%---------------------
radioh = uibuttongroup('Parent', IPpanelH, 'Visible','off',...
                  'Units', 'normalized','BackgroundColor',[1 1 1],...
                  'Position',[.01 .50 .95 .24],...
                  'SelectionChangedFcn',@ROImode);

% Create two radio buttons in the button group.
pairedRadio = uicontrol(radioh,'Style','radiobutton','FontSize', 11,...
                  'String','Paired ROI Mode','BackgroundColor',[1 1 1],...
                  'Units', 'normalized',...
                  'Position',[0.05 0.05 0.32 0.9],...
                  'HandleVisibility','off');

singleRadio = uicontrol(radioh,'Style','radiobutton','FontSize', 11,...
                  'String','Single ROI Mode','BackgroundColor',[1 1 1],...
                  'Units', 'normalized',...
                  'Position',[0.40 0.05 0.32 0.9],...
                  'HandleVisibility','off');

doLineTrace = uicontrol(radioh,'Style','checkbox','FontSize', 11,...
                  'String','Line Trace','BackgroundColor',[1 1 1],...
                  'Units', 'normalized',...
                  'Position',[0.75 0.05 0.23 0.9],...
                  'HandleVisibility','off');              
radioh.Visible = 'on';
doLineTrace.Value = 1;
%------------------


loadROISh = uicontrol('Parent', IPpanelH, 'Units', 'normalized', ...
    'Position', [0.05 0.25 0.29 0.20], 'String', 'Load ROIs', 'FontSize', 11,...
    'Callback', @loadROIfun);

saveROISh = uicontrol('Parent', IPpanelH, 'Units', 'normalized', ...
    'Position', [0.36 0.25 0.29 0.20], 'String', 'Save ROIs', 'FontSize', 11,...
    'Callback', @saveROIfun);


exploreFLIMh = uicontrol('Parent', IPpanelH, 'Units', 'normalized', ...
    'Position', [0.67 0.25 0.29 0.20], 'String', 'Explore FLIM', 'FontSize', 11,...
    'Callback', @exploreFLIM);


savefileh = uicontrol('Parent', IPpanelH, 'Units', 'normalized', ...
    'Position', [0.05 0.02 0.45 0.20], 'String', 'Save Data', 'FontSize', 11,...
    'Callback', @saveDataset);


loadIMGh = uicontrol('Parent', IPpanelH, 'Units', 'normalized', ...
    'Position', [0.52 0.02 0.45 0.20], 'String', 'Import Imaging Data', 'FontSize', 11,...
    'Callback', @loadIMG);





%----------------------------------------------------
%%     RIGHT PANE FIGURE PANELS
%----------------------------------------------------

tabgp = uitabgroup(MAINGUIFIG,'Position',[.62 .27 .37 .42],'TabLocation','top');
FLIMtab = uitab(tabgp,'Title','FLIM','BackgroundColor',[1 1 1]);
PLOTtab = uitab(tabgp,'Title','PLOT','BackgroundColor',[1 1 1]);


haxTABH = axes('Parent', PLOTtab, 'NextPlot', 'Add','Color','none',...
    'Position', [.10 .10 0.8 0.8]); % 'XColor','none','YColor','none'


%----------------------------------------------------
%           FLIM MIN MAX PANEL
%----------------------------------------------------

FLIMminmaxPanel = uipanel('Parent', FLIMtab,'Title','FLIM Min Max Panel','FontSize',10,...
    'BackgroundColor',[1 1 1],...
    'Position', [0 0 1 1]); % 'Visible', 'Off',


setintenh = uicontrol('Parent', FLIMminmaxPanel, 'Units', 'normalized', ...
    'Position', [.02 .89 .45 .10], 'FontSize', 11, 'String', 'Set intensity',...
    'Callback', @setinten);


dftintenh = uicontrol('Parent', FLIMminmaxPanel, 'Units', 'normalized', ...
    'Position', [.52 .89 .45 .10], 'FontSize', 11,...
    'String', 'Default intensities','Callback', @defaultinten);


intThreshMinh = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Text', 'Units', 'normalized', ...
    'Position', [.02 .80 .45 .06], 'FontSize', 11, 'String', 'Min Intensity');
intThreshMin = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Edit', 'FontSize', 11, 'Units', 'normalized', ...
    'Position', [.02 .72 .45 .08]);

intThreshMaxh = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Text', 'Units', 'normalized', ...
    'Position', [.52 .80 .45 .06], 'FontSize', 11, 'String', 'Max Intensity');
intThreshMax = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Edit', 'FontSize', 11, 'Units', 'normalized', ...
    'Position', [.52 .72 .45 .08]);


lifetimethresholdh = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Text',  'Units', 'normalized',...
    'Position', [.02 .60 .45 .06], 'FontSize', 11, 'String', 'Lifetime Min');
lftthresholdMINh = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Edit',  'Units', 'normalized',...
    'Position', [.02 .52 .45 .08], 'FontSize', 11);


lifetimethreshMAXh = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Text',  'Units', 'normalized',...
    'Position', [.52 .60 .45 .06], 'FontSize', 11, 'String', 'Lifetime Max');
lftthresholdMAXh = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Edit',  'Units', 'normalized',...
    'Position', [.52 .52 .45 .08], 'FontSize', 11);



chithresholdminh = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Text',  'Units', 'normalized',...
    'Position', [.02 .40 .45 .06], 'FontSize', 11, 'String', 'Chi Min');
chiminh = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Edit',  'Units', 'normalized', ...
    'Position', [.02 .32 .45 .08], 'FontSize', 11);


chithresholdmaxh = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Text',  'Units', 'normalized', ...
    'Position', [.52 .40 .45 .06], 'FontSize', 11, 'String', 'Chi Max');
chimaxh = uicontrol('Parent', FLIMminmaxPanel, 'Style', 'Edit',  'Units', 'normalized', ...
    'Position', [.52 .32 .45 .08], 'FontSize', 11);



%----------------------------------------------------
%           MEMO CONSOLE GUI WINDOW
%----------------------------------------------------

memopanelH = uipanel('Parent', MAINGUIFIG,'Title','Memo Log ','FontSize',10,...
    'BackgroundColor',[1 1 1],...
    'Position', [0.62 0.01 0.37 0.25]); % 'Visible', 'Off',


memes = {' ',' ',' ', ' ',' ',' ',' ', ...
         'Welcome to FLIMIM!', 'Getting things ready for you...'};

conboxH = uicontrol('Parent',memopanelH,'Style','listbox','Units','normalized',...
        'Max',9,'Min',0,'Value',9,'FontSize', 13,'FontName', 'FixedWidth',...
        'String',memes,'FontWeight', 'bold',...
        'Position',[.0 .0 1 1]);  
    


%----------------------------------------------------
%     IMPORT IMAGE & LOAD DEFAULT TOOLBOX PARAMETERS
%----------------------------------------------------

% loadfile()

disableButtons()
loadIMGh.Enable = 'on';

set(MAINGUIFIG, 'Visible', 'On');

memocon(' ');
memocon('All set.')

axes(haxFLIM)
axes(haxROIS)

memocon(' '); pause(.2);
memocon('Click "Import Imaging Data" button to begin.')

pause(.5);

loadIMGh.BackgroundColor = [.05 .95 .50]; pause(.15)
loadIMGh.BackgroundColor = [0.94 0.94 0.94];  pause(.10)
loadIMGh.BackgroundColor = [.05 .95 .50]; pause(.10)
loadIMGh.BackgroundColor = [0.94 0.94 0.94]; pause(.1);






% -----------------------------------------------------------------------------
%%                     GUI TOOLBOX FUNCTIONS
% -----------------------------------------------------------------------------


%----------------------------------------------------
%        LOAD IMAGING DATA
%----------------------------------------------------
function loadIMG(loadIMGh, eventData)
    
    
    
    memocon('Select a .dat file'); pause(.5);
    [datfilef, datfilefdir] = uigetfile('*.dat', 'load a .dat file');
    
    cd(datfilefdir);
    imgpath = [datfilefdir '/' datfilef];
    
    
    memocon('Select red channel image'); pause(.5);
    [ccdfiler, ccdfilerdir] = uigetfile({'*.tif*';'*.bmp';'*.jpg';'*.png'}, 'load green-channel .bmp file');
    
    
    memocon('Select green channel image'); pause(.5);
    [ccdfileg, ccdfilegdir] = uigetfile({'*.tif*';'*.bmp';'*.jpg';'*.png'}, 'load red-channel .dat file');
    
    

    
    
    tempdata = load([datfilefdir datfilef]);
    tempdatadim = size(tempdata);
    totxdim = tempdatadim(1);
    
    
    ydim = tempdatadim(2);
    if mod(totxdim,3)~=0
        disp('This does not appear to be a properly compiled file.');
        return
    end
    xdim = totxdim/3;
    datastack = zeros(xdim,ydim,3,'double');

    
    set(MAINGUIFIG, 'Visible', 'On');

    datastack(1:xdim,1:ydim,1) = tempdata(1:xdim,1:ydim);
    datastack(1:xdim,1:ydim,2) = tempdata(xdim+1:2*xdim,1:ydim);
    datastack(1:xdim,1:ydim,3) = tempdata(2*xdim+1:3*xdim,1:ydim);

    lifetime = datastack(:,:,1);
    intensity = datastack(:,:,2);
    chi = datastack(:,:,3);
    
    
    imgG = imformat([ccdfilegdir ccdfileg]);
    imgR = imformat([ccdfilerdir ccdfiler]);
    
    
    
    if size(lifetime,1) > size(imgG,1)
    
        lifetime = lifetime(1:size(imgG,1),:);
        intensity = intensity(1:size(imgG,1),:);
        chi = chi(1:size(imgG,1),:);
        
    elseif size(lifetime,1) < size(imgG,1)
        
        imgG = imgG(1:size(lifetime,1),:);
        imgR = imgR(1:size(lifetime,1),:);
        
    end
        
    if size(lifetime,2) > size(imgG,2)
    
        lifetime = lifetime(:,1:size(imgG,2));
        intensity = intensity(:,1:size(imgG,2));
        chi = chi(:,1:size(imgG,2));
        
    elseif size(lifetime,2) < size(imgG,2)
        
        imgG = imgG(:,1:size(lifetime,2));
        imgR = imgR(:,1:size(lifetime,2));
        
    end
        
        imgF = intensity;
        
      

    %----------------------------------------------------
    %           SET USER-EDITABLE GUI VALUES
    %----------------------------------------------------
    set(intThreshMin, 'String', num2str(intenseThreshMIN));
    set(intThreshMax, 'String', num2str(intenseThreshMAX));

    set(intThreshMin, 'String', num2str(intenseThreshMIN));
    set(intThreshMax, 'String', num2str(intenseThreshMAX));

    set(lftthresholdMINh, 'String', num2str(lifeThreshMIN));
    set(lftthresholdMAXh, 'String', num2str(lifeThreshMAX));

    set(chiminh, 'String', num2str(chiThreshMIN));
    set(chimaxh, 'String', num2str(chiThreshMAX));

    
    ROI_IDh.String = int2str(0);
    
    radioh.SelectedObject = pairedRadio; % Set radiobutton to stamp
    
    set(MAINGUIFIG, 'Name', datfilef);
    
    % set(magh, 'String', num2str(magnification));
    % boxtype = boxtypeh.SelectedObject.String;
    
    %----------------------------------------------------
    %           AXES LIMITS
    %----------------------------------------------------
    
    set(haxCCDG, 'XLim', [1 size(imgG,2)], 'YLim', [1 size(imgG,1)]);
    set(haxCCDR, 'XLim', [1 size(imgR,2)], 'YLim', [1 size(imgR,1)]);
    set(haxFLIM, 'XLim', [1 xdim]);
    set(haxFLIM, 'YLim', [1 ydim]);
    
    haxROIS.XLim = [1 xdim];
    haxROIS.YLim = [1 ydim];
    
    %----------------------------------------------------
    
    
    %----------------------------------------------------
    %                   DRAW IMAGE
    %----------------------------------------------------
    
    axes(haxCCDG)
    colormap(haxCCDG,hot)
    phCCDG = imagesc(imgG , 'Parent', haxCCDG);
    memocon('haxCCDG: Green Channel')
              pause(.8)
              
    axes(haxCCDR)
    colormap(haxCCDR,hot)
    phCCDR = imagesc(imgR, 'Parent', haxCCDR);
    memocon('haxCCDR: Red Channel')
        pause(.8)
        
    axes(haxFLIM)    
    colormap(haxFLIM,hot)
    phFLIM = imagesc(intensity, 'Parent', haxFLIM,...
                  [prctile(intensity(:),intenseThreshMIN) prctile(intensity(:),intenseThreshMAX)]);
        memocon('haxFLIM: FLIM Channel')
        pause(.8)
    
    imXlim = haxFLIM.XLim;
    imYlim = haxFLIM.YLim;
    
    
    axes(haxROIS)
    

%{

iminfo = imfinfo(imgpath);
[im, map] = imread(imgpath);


im_size = size(im);
im_nmap = numel(map);
im_ctype = iminfo.ColorType;


if strcmp(im_ctype, 'truecolor') || numel(im_size) > 2

IMG = rgb2gray(im);
IMG = im2double(IMG);

elseif strcmp(im_ctype, 'indexed')

IMG = ind2gray(im,map);
IMG = im2double(IMG);

elseif strcmp(im_ctype, 'grayscale')

IMG = im2double(im);

else

IMG = im;

end


axes(haxCCD)
colormap(haxCCD,bone); % parula
phCCD = imagesc(IMG , 'Parent', haxCCD);
      pause(1)

ccmap = bone;
cmmap = [zeros(10,3); ccmap(end-40:end,:)];
colormap(haxCCD,cmmap)
MAINGUIFIG.Colormap = cmmap;



pause(.2)
imXlim = haxCCD.XLim;
imYlim = haxCCD.YLim;


xdim = size(IMG,2); 
ydim = size(IMG,1);



%----------------------------------------------------
%           SET USER-EDITABLE GUI VALUES
%----------------------------------------------------
set(MAINGUIFIG, 'Name', datafile);
set(ROIIDh, 'String', int2str(1));
set(haxCCD, 'XLim', [1 xdim]);
set(haxCCD, 'YLim', [1 ydim]);
%----------------------------------------------------
% axes(haxCCD)
%}    
    
    enableButtons()
    
    memocon(' ')
    memocon('Notes:')
    memocon('Click on the image to make it active.')
    memocon('Use the =/- keys to zoom in and out')
    memocon('Use s d f e keys to move left down right up')
    memocon(' ')
    
    
    pause(.2)
    memocon('Import previously saved ROI data.')
    pause(.2)
    loadROISh.BackgroundColor = [.05 .95 .50]; pause(.15)
    loadROISh.BackgroundColor = [0.94 0.94 0.94]; pause(.05)
    
    memocon('OR...'); memocon(' ')
    
    pause(.5)
    memocon('Click MEASURE ROI to begin data collection.')
    pause(.5)
    getROIh.BackgroundColor = [.05 .95 .50]; pause(.15)
    getROIh.BackgroundColor = [0.94 0.94 0.94]; pause(.1)
    getROIh.BackgroundColor = [.05 .95 .50]; pause(.1)
    getROIh.BackgroundColor = [0.94 0.94 0.94]; pause(.1);
    
end








%----------------------------------------------------
%           GET ROI
%----------------------------------------------------
function getROI(boxidselecth, eventdata)

    ROI_IDh.String = int2str(str2num(ROI_IDh.String) + 1);
    Rmode = radioh.SelectedObject.String;
    doTrace = doLineTrace.Value;
        
    
    if strcmp(Rmode,'Paired ROI Mode')
    %---------------------------
    % GET SPINE ROI
    %---------------------------
        memocon('Click center of ROI #1 (spine)')
        %axes(haxROIS); pause(0.05);
        axes(haxAPI); APIgroup.Visible='on'; pause(0.05);
    
        hROI = impoint;
        ROIpos = hROI.getPosition;
        delete(hROI)
        hROI = imellipse(haxAPI, [ROIpos-round(stampSize/2) stampSize stampSize]);
        ROIpos = hROI.getPosition;
        ROIarea = pi * (stampSize/2)^2;
        setColor(hROI,[0 0 1]);
        % verts = getVertices(hROI);

        haxAPI.Children(1).DisplayName = ['S' ROI_IDh.String];

        %nROI = numel(ROI)+1;
        %ROI(nROI).POS   = ROIpos;
        %ROI(nROI).PH    = hROI;
        %ROI(nROI).TYPE  = Rmode;

        
        
        set(hROI,'Parent',APIgroup)
        APIgroup.Visible = 'off';
        axes(haxROIS); pause(0.05);
        %axes(haxAPI); APIgroup.Visible='on'; pause(0.05);
        
        hREC = rectangle('Parent',haxROIS,'Position',ROIpos,'Curvature',1,'LineWidth',2,'EdgeColor',[.4 .4 1]);
        set(hREC,'Parent',ROIgroup)
        
    %---------------------------
    % GET DENDRITE ROI
    %---------------------------
        memocon('Click center of ROI #2 (dendrite)')
        axes(haxAPI); APIgroup.Visible='on'; pause(0.05);
        
        hROI = impoint;
        ROIpos = hROI.getPosition;
        delete(hROI)
        hROI = imellipse(haxAPI, [ROIpos-round(stampSize/2) stampSize stampSize]);
        ROIpos = hROI.getPosition;
        ROIarea = pi * (stampSize/2)^2;
        setColor(hROI,[1 1 0]);

        haxAPI.Children(1).DisplayName = ['D' ROI_IDh.String];
        
        set(hROI,'Parent',APIgroup)
        APIgroup.Visible = 'off';
        axes(haxROIS); pause(0.05);
        %axes(haxAPI); APIgroup.Visible='on'; pause(0.05);
        
        hREC = rectangle('Parent',haxROIS,'Position',ROIpos,'Curvature',1,'LineWidth',2,'EdgeColor',[1 1 0]);
        set(hREC,'Parent',ROIgroup)

    else
    %---------------------------
    % GET SPINE ROI
    %--------------------------- 
        memocon('Click center of ROI')
        %axes(haxROIS); pause(0.05);
        axes(haxAPI); APIgroup.Visible='on'; pause(0.05);
    
        hROI = impoint;
        ROIpos = hROI.getPosition;
        delete(hROI)
        hROI = imellipse(haxAPI, [ROIpos-round(stampSize/2) stampSize stampSize]);
        ROIpos = hROI.getPosition;
        ROIarea = pi * (stampSize/2)^2;
        setColor(hROI,[0 0 1]);
        % verts = getVertices(hROI);

        haxAPI.Children(1).DisplayName = ['S' ROI_IDh.String];

        nROI = numel(ROI)+1;
        ROI(nROI).POS   = ROIpos;
        ROI(nROI).PH    = hROI;
        ROI(nROI).TYPE  = Rmode;

        
        
        set(hROI,'Parent',APIgroup)
        APIgroup.Visible = 'off';
        axes(haxROIS); pause(0.05);
        %axes(haxAPI); APIgroup.Visible='on'; pause(0.05);
        
        hREC = rectangle('Parent',haxROIS,'Position',ROIpos,'Curvature',1,'LineWidth',2,'EdgeColor',[.4 .4 1]);
        set(hREC,'Parent',ROIgroup)
    
    end
    
    
    
    
    %---------------------------
    % GET LINE-TRACE ROI
    %---------------------------
    if doTrace
        memocon('Click drag a line - release to finish')
        axes(haxAPI); TRACEgroup.Visible='on'; pause(0.05);

        hROI = imline(haxAPI);
        dpos = hROI.getPosition;
        setColor(hROI,[0 1 0]);
        haxAPI.Children(1).DisplayName = ['L' ROI_IDh.String];
        
        spineextent = sqrt((dpos(1,1)-dpos(2,1))^2 + (dpos(1,2)-dpos(2,2))^2);
        
        [SPIx,SPIy,SPIf] = improfile(imgF, dpos(:,1), dpos(:,2), round(spineextent));
        [SPIx,SPIy,SPIg] = improfile(imgG, dpos(:,1), dpos(:,2), round(spineextent));
        [SPIx,SPIy,SPIr] = improfile(imgR, dpos(:,1), dpos(:,2), round(spineextent));
        
        fmax = max(max(SPIf)); 
        fmin = min(min(SPIf));
        SPIft = lintrans(SPIf,fmin,fmax,0,1);

        axes(haxTABH); hold off;
        plot(haxTABH, SPIft); hold on;
        plot(haxTABH, SPIg); hold on;
        plot(haxTABH, SPIr)
        axes(haxROIS)
        
        
        set(hROI,'Parent',TRACEgroup)
        TRACEgroup.Visible = 'off';
        axes(haxROIS); pause(0.05);
        
        hLINE = line(dpos(:,1), dpos(:,2),'Parent',haxROIS,'Color',[0 1 0]);
        set(hLINE,'Parent',ROIgroup)
        
    end

    axes(haxROIS);


    %---------------------------
    % GET ANOTHER ROI
    %---------------------------
    doagainROI = questdlg('Select next ROI?', 'Select next ROI?', 'Yes', 'No', 'No');
    switch doagainROI
       case 'Yes'
            getROI
       case 'No'
           
    end

    set(gcf,'Pointer','arrow')

end






%----------------------------------------------------
%           GET MOUSE LOCATION
%----------------------------------------------------
function GetMouseLoc(boxidselecth, eventdata)

% set(gcf,'Pointer','hand')

        if(saveROI(str2double(get(ROI_IDh, 'String')),1)==0)
            %[x, y] = ginput(2);
            [x,y] = FLIMginput(2,'custom');
            x1=x(1);
            y1=y(1);
            x2=x(2);
            y2=y(2);
            calcROIcoor(x1, y1, x2, y2, str2double(get(ROI_IDh, 'String')));
        else 
            duplicateROI = questdlg('Box already exists. Overwrite?', 'Duplicate ROI', 'Yes', 'No', 'No');
            switch duplicateROI
                case 'Yes'
                    [x, y] = ginput(2);
                    x1=x(1);
                    y1=y(1);
                    x2=x(2);
                    y2=y(2);
                    calcROIcoor(x1, y1, x2, y2, str2double(get(ROI_IDh, 'String')));
                case 'No'
            end
        end

        doagainROI = questdlg('Select next ROI?', 'Select next ROI?', 'Yes', 'No', 'No');
        switch doagainROI
           case 'Yes'
                set(ROI_IDh,'String',num2str((str2num(ROI_IDh.String)+1)) );
                GetMouseLoc
           case 'No'
        end

set(gcf,'Pointer','arrow')

end








%----------------------------------------------------
%           LIFETIME PREVIEW
%----------------------------------------------------
function exploreFLIM(lifetimeviewerh, eventData)

    set(MAINGUIFIG, 'Visible', 'Off');
    % set(initmenuh, 'Visible', 'Off');

    lftthresholdMIN = str2double(get(lftthresholdMINh, 'String'));
    lftthresholdMAX = str2double(get(lftthresholdMAXh, 'String'));
    
    intenthresholdMIN = str2double(get(intThreshMin, 'String'));
    intenthresholdMAX = str2double(get(intThreshMax, 'String'));
    
    intPminmax = prctile(intensity(:),[intenthresholdMIN intenthresholdMAX]);
    
    chimin = str2double(get(chiminh, 'String'));
    chimax = str2double(get(chimaxh, 'String'));


    ChiG = (chi >= chimin & chi <= chimax);
    IntG = (intensity >= intPminmax(1) & intensity <= intPminmax(2));
    LifG = (lifetime >= lftthresholdMIN & lifetime <= lftthresholdMAX);
    AllG = (ChiG .* IntG .* LifG) > 0;

        
    % close all
    fh3=figure('Units','normalized','OuterPosition',[.05 .27 .9 .7],'Color','w');
    ah1 = axes('Position',[.05 .55 .2 .4],'Color','none','XTick',[],'YTick',[]);
    ah2 = axes('Position',[.30 .55 .2 .4],'Color','none','XTick',[],'YTick',[]);
    ah3 = axes('Position',[.05 .05 .2 .4],'Color','none','XTick',[],'YTick',[]);
    ah4 = axes('Position',[.30 .05 .2 .4],'Color','none','XTick',[],'YTick',[]);
    ah5 = axes('Position',[.55 .15 .40 .7],'Color','none','XTick',[],'YTick',[]);
    
        axes(ah1)
    imagesc(ChiG); title('Pixels within CHI thresholds')
        axes(ah2)
    imagesc(IntG); title('Pixels within INTENSITY thresholds')
        axes(ah3)
    imagesc(LifG); title('Pixels within LIFETIME thresholds')
        axes(ah4)
    imagesc(AllG); title('Pixels within ALL thresholds')
        axes(ah5)
    imagesc(lifetime .* AllG); title('Fluorescent Lifetime of pixels above ALL thresholds')
        colormap(ah5,[0 0 0; flipud(jet(15))])
        caxis([1600 2800])
        colorbar
        set(ah1,'YDir','normal')
        set(ah2,'YDir','normal')
        set(ah3,'YDir','normal')
        set(ah4,'YDir','normal')
        set(ah5,'YDir','normal')

    disp('Close figure to continue')
    uiwait(fh3)
    
    
    set(MAINGUIFIG, 'Visible', 'On');

end






%----------------------------------------------------
%           SET INTENSITY
%----------------------------------------------------
function setinten(hObject, eventdata)
    
       lowerinten = str2num(intThreshMin.String);
       upperinten = str2num(intThreshMax.String);
       
       lowerintenPCT = prctile(intensity(:),lowerinten);
       upperintenPCT = prctile(intensity(:),upperinten);
       
       
       set(haxFLIM,'CLim',[lowerintenPCT upperintenPCT])

end




%----------------------------------------------------
%        SET CLIM SLIDER CALLBACK
%----------------------------------------------------
function setCLimFLIM(hObject, eventdata)
% Hints: hObject.Value returns position of slider
%        hObject.Min and hObject.Max determine range of slider
% sunel = get(handles.CLIMslider,'value'); % Get current light elev.
% sunaz = get(hObject,'value');   % Varies from -180 -> 0 deg


slideCON = CONslider.Value;

lowerinten = str2num(intThreshMin.String);
upperinten = str2num(intThreshMax.String);
lowerintenPCT = prctile(intensity(:),lowerinten);
upperintenPCT = prctile(intensity(:),upperinten);


slideVal = ceil(CLIMslider.Value);

haxFLIM.CLim = [(lowerintenPCT / ((slideVal/50)*(slideVal/50))),...
    (upperintenPCT / ((slideVal/50)*(slideVal/50)) )];

% memocon(['image' num2str(slideVal)])

end

function setCLimRFP(hObject, eventdata)
    
    slideCON = CONsliderCCDR.Value;

    lowerintenPCT = prctile(imgR(:),slideCON);
    upperintenPCT = prctile(imgR(:),99.999);

    slideVal = ceil(CLIMsliderCCDR.Value);

    haxCCDR.CLim = [(lowerintenPCT / ((slideVal/50)*(slideVal/50))),...
                    (upperintenPCT / ((slideVal/50)*(slideVal/50)) )];

    % memocon(['image' num2str(slideVal)])

end

function setCLimGFP(hObject, eventdata)
    
    slideCON = CONsliderCCDG.Value;

    lowerintenPCT = prctile(imgG(:),slideCON);
    upperintenPCT = prctile(imgG(:),99.999);

    slideVal = ceil(CLIMsliderCCDG.Value);

    haxCCDG.CLim = [(lowerintenPCT / ((slideVal/50)*(slideVal/50))),...
                    (upperintenPCT / ((slideVal/50)*(slideVal/50)) )];

    % memocon(['image' num2str(slideVal)])

end











%----------------------------------------------------
%           SET COLORMAP
%----------------------------------------------------
function setcolormap
              
    % set(MAINGUIFIG, 'Colormap', gray);
    set(MAINGUIFIG, 'Colormap', hot);
    % colormap([0 0 0; jet(20)])
        
end






%----------------------------------------------------
%           DEFAULT INTENSITY
%----------------------------------------------------
function defaultinten(hObject, eventdata)
        
       set(intThreshMin, 'String', num2str(intenseThreshMIN));
       set(intThreshMax, 'String', num2str(intenseThreshMAX));

end






%----------------------------------------------------
%           CLOSE INTENSITY WINDOW
%----------------------------------------------------
function closelftintenw(hObject, eventdata)
%Closelftintenw sets both lifetime image and intensity image windows to
%invisible. The initial menu becomes visible again for further selection. 
    
       set(MAINGUIFIG, 'Visible', 'Off');
       set(initmenuh, 'Visible', 'On');
       saveROI = zeros(200, 17);
       saveData = zeros(200, 9);
       datastack = zeros(1,1,3,'double');
       lifetime = zeros(1, 1);
       intensity = zeros(1, 1);
       chi = zeros(1, 1);
       lifetimeimage = zeros(1, 1);
       intensityimage = zeros(1, 1);
       xdim = 0;
       ydim = 0;
end













%----------------------------------------------------
%           KEY PRESS ZOOM
%----------------------------------------------------
function keypresszoom(hObject, eventData, key)
    
    

    
        % --- ZOOM ---
        
        if strcmp(MAINGUIFIG.CurrentCharacter,'=')
            
            % IN THE FUTURE USE MOUSE LOCATION TO ZOOM
            % INTO A SPECIFIC POINT. TO QUERY MOUSE LOCATION
            % USE THE METHOD: MAINGUIFIG.CurrentPoint
            
            zoom(1.5)
            pause(.05)
        end
        
        if strcmp(MAINGUIFIG.CurrentCharacter,'-')
            zoom(.5)
            pause(.05)
        end
                
        
        % --- PAN ---
        
        if strcmp(MAINGUIFIG.CurrentCharacter,'p')

            pan('on')        
            % h = pan;
            % h.ActionPreCallback = @myprecallback;
            % h.ActionPostCallback = @mypostcallback;
            % h.Enable = 'on';
            drawnow; % pause(.1)
        end
        if strcmp(MAINGUIFIG.CurrentCharacter,'o')
            pan('off')        
            pause(.05)
        end
        
        if strcmp(MAINGUIFIG.CurrentCharacter,'f')
            %haxCCDG.XLim = haxCCDG.XLim+20;
            %haxCCDR.XLim = haxCCDR.XLim+20;
            haxFLIM.XLim = haxFLIM.XLim+20;
            %haxROIS.XLim = haxROIS.XLim+20;
            pause(.05)
        end
        
        if strcmp(MAINGUIFIG.CurrentCharacter,'s')
            %haxCCDG.XLim = haxCCDG.XLim-20;
            %haxCCDR.XLim = haxCCDR.XLim-20;
            haxFLIM.XLim = haxFLIM.XLim-20;
            %haxROIS.XLim = haxROIS.XLim-20;
            pause(.05)
        end
        
        if strcmp(MAINGUIFIG.CurrentCharacter,'e')
            %haxCCDG.YLim = haxCCDG.YLim+20;
            %haxCCDR.YLim = haxCCDR.YLim+20;
            haxFLIM.YLim = haxFLIM.YLim+20;
            %haxROIS.YLim = haxROIS.YLim+20;
            pause(.05)
        end
        
        if strcmp(MAINGUIFIG.CurrentCharacter,'d')
            %haxCCDG.YLim = haxCCDG.YLim-20;
            %haxCCDR.YLim = haxCCDR.YLim-20;
            haxFLIM.YLim = haxFLIM.YLim-20;
            %haxROIS.YLim = haxROIS.YLim-20;
            pause(.05)
        end
        
        
        % --- RESET ZOOM & PAN ---
        
        if strcmp(MAINGUIFIG.CurrentCharacter,'0')
            zoom out
            zoom reset
            haxCCDG.XLim = imXlim;
            haxCCDG.YLim = imYlim;
            haxCCDR.XLim = imXlim;
            haxCCDR.YLim = imYlim;
            haxFLIM.XLim = imXlim;
            haxFLIM.YLim = imYlim;
            haxROIS.XLim = imXlim;
            haxROIS.YLim = imYlim;
            pause(.05)
        end
        
        
end






%----------------------------------------------------
%           RADIO BUTTON BOX SELECTION
%----------------------------------------------------
function ROImode(source,callbackdata)
    
    % callbackdata.OldValue.String
    % radioh.SelectedObject.String
    % boxtype = callbackdata.NewValue.String;
    
    memocon(['ROI selection mode: ' callbackdata.NewValue.String])

end






%----------------------------------------------------
%           COMPILE DATA FILE
%----------------------------------------------------
function compilefile(hObject, eventdata)
%Compile function triggers the user input for three datafiles. (Will return
%error if files are not of the right size, etc.) Compile then writes a new 
%file that is the stacked version of all three files. 


    home = cd;
    
    memocon('Select a folder with files to compile.')
    compiledir = uigetdir;
    cd(compiledir);

    filelist = dir(compiledir);
    filelistsize = size(filelist);
        

        
    for ii=1:filelistsize

        currentfile = filelist(ii,1).name;

        if(length(currentfile) >=7)

            if (strcmp('Chi of ', currentfile(1:7))==1) ||...
               (strcmp('chi.asc', currentfile(end-6:end))==1)

                chifile = currentfile;
                if (strcmp('Chi of ', chifile(1:7))==1)
                    chifilename = chifile(8:end);
                else
                    chifilename = chifile(1:end-8);
                end


                for jj=1:filelistsize

                    searchcolorfile = filelist(jj,1).name;

                    if(length(searchcolorfile) >=21)


                        if (strcmp('Color coded value of ', searchcolorfile(1:21))==1) ||...
                           (strcmp('color coded value.asc', searchcolorfile(end-20:end))==1)

                            if (strcmp(chifilename, searchcolorfile(22:end))==1) ||...
                               (strcmp(chifilename, searchcolorfile(1:end-22))==1)

                                colorfile = searchcolorfile;
                                if (strcmp(chifilename, colorfile(22:end))==1)
                                    colorfilename = colorfile(22:end);
                                else
                                    colorfilename = colorfile(1:end-22);
                                end


                                for kk=1:filelistsize
                                    searchintenfile = filelist(kk,1).name;

                                    if(length(searchintenfile) >=12)

                                        if (strcmp('Photons of ', searchintenfile(1:11))==1) ||...
                                           (strcmp('photons.asc', searchintenfile(end-10:end))==1)


                                            if (strcmp(chifilename, searchintenfile(12:end))==1) ||...
                                               (strcmp(chifilename, searchintenfile(1:end-12))==1)

                                                intenfile = searchintenfile;
                                                if (strcmp(chifilename, colorfile(22:end))==1)
                                                    intenfilename = intenfile(12:end);
                                                else
                                                    intenfilename = intenfile(1:end-12);
                                                end


                                                if(strcmp(chifilename, colorfilename)==1 && strcmp(chifilename, intenfilename)==1)

                                                    lifetime = load(colorfile);
                                                    lifetimedim = size(lifetime);
                                                    intensity = load(intenfile);
                                                    intensitydim = size(intensity);
                                                    chitemp = load(chifile);
                                                    chidim = size(chitemp);
                                                    chi = zeros(chidim(1), chidim(2));
                                                    chitemp(chitemp>100) = 0;
                                                    chi = chitemp;
                                                    
                                                    if isequal(lifetimedim, intensitydim, chidim)==1    
                                                        savefilename = mat2str(strcat(chifilename, '.dat'));
                                                        savefilename = savefilename(2:end-1);
                                                        save(savefilename, 'lifetime', '-ascii');
                                                        save(savefilename, 'intensity', '-ascii', '-append');
                                                        save(savefilename, 'chi', '-ascii', '-append');

                                                        memocon(['Successfully compiled: ' savefilename])
                                                        
                                                    else
                                                        savefilename = mat2str(strcat(chifilename, '.dat'));
                                                        savefilename = savefilename(2:end-1);
                                                        memocon(['Error compiling ' savefilename])
                                                        memocon('One or more component files may be incorrect.')

                                                    end

                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end                 
        end
    end
    memocon('ALL FILES SUCCESSFULLY COMPILED!')
    cd(home);

end
 






%----------------------------------------------------
%        CHANGE IMAGE
%----------------------------------------------------
function resetwsFun(hObject, eventdata)

    resetFZ = questdlg('Reset FLIMZ Toolbox?', 'Reset FLIMZ Toolbox?', 'Yes', 'No', 'Yes');
    switch resetFZ
       case 'Yes'
        memocon('Spinning-up a fresh FLIMZ workspace.');
        clc; close all; clear;
        FLIMZ
       case 'No'
        memocon(' ');
        memocon(' ');
        memocon('KEEP'); pause(.7)
        memocon('      IT'); pause(.7)
        memocon('          ROLLIN'' !!!'); pause(.7)
        memocon('                            ');
        memocon('                            ');
        memocon('     _________              ');
        memocon('      O     O               ');
        memocon('----------------------------'); pause(.05)
        memocon('                            ');
        memocon('                            ');
        memocon('        _________           ');
        memocon('         O     O            ');
        memocon('----------------------------'); pause(.05)
        memocon('                            ');
        memocon('                            ');
        memocon('           _________        ');
        memocon('            O     O         ');
        memocon('----------------------------'); pause(.05)
        memocon('                            ');
        memocon('                            ');
        memocon('               _________    ');
        memocon('                O     O     ');
        memocon('----------------------------'); pause(.05)
        memocon('                            ');
        memocon('                            ');
        memocon('                   _________');
        memocon('                    O     O ');
        memocon('----------------------------'); pause(.05)
        memocon('                            ');
        memocon('                          / ');
        memocon('                         /  ');
        memocon('                        /O  ');
        memocon('----------------------------'); pause(.05)
        memocon('                            ');
        memocon('                            ');
        memocon('                          __');
        memocon('                           O');
        memocon('----------------------------'); pause(.05)
        memocon('                            ');
        memocon('                            ');
        memocon('                            ');
        memocon('                            ');
        memocon('                            ');
        memocon('                            ');
        memocon('                            ');
        memocon('                            ');
        memocon('                            ');
        memocon('Keep it rollin''. ');
    end
    

end




%----------------------------------------------------
%        CHANGE IMAGE
%----------------------------------------------------
function changeimg(hObject, eventdata)

    haxnum = circshift(haxnum,[0 -1]);
    
    if haxnum(1) == 1
    
        axes(haxFLIM)
        axes(haxROIS)
        memocon('haxFLIM: FLIM Channel')
        
    
    elseif haxnum(1) == 2
        
        axes(haxCCDG)
        axes(haxROIS)
        memocon('haxCCDG: Green Channel')
        
    
    else
    
        axes(haxCCDR)
        axes(haxROIS)
        memocon('haxCCDR: Red Channel')
    
    end
    



end




%----------------------------------------------------
%        TOGGLE DISPLAYED IMAGE
%----------------------------------------------------
function imTOG(hObject, eventdata)

        
    if showFIMh.Value == 1;

        axes(haxFLIM)
        axes(haxROIS)
        memocon('FLIM channel IMG')
        pause(0.05);

    elseif showGIMh.Value == 1;

        axes(haxCCDG)
        axes(haxROIS)
        memocon('GREEN channel IMG')

    elseif showRIMh.Value == 1;

        axes(haxCCDR)
        axes(haxROIS)
        memocon('RED channel IMG')

    end
    
end










%----------------------------------------------------
%        DELETE ROI MODE
%----------------------------------------------------
function delROImode(hObject, eventdata)


    if strcmp(hObject.Tag,'delModeTOG')
        
        if delROImodeON.Value == 1;
            
            UITB_delROI.State = 'on';

            memocon('ROI DELETION MODE ~~ON~~')
            disableButtons
            axes(haxAPI); 
            APIgroup.Visible='on';
            TRACEgroup.Visible='on';
            pause(0.05);

        else
            
            UITB_delROI.State = 'off';
            
            redrawROIs()

            memocon('ROI DELETION MODE ~~OFF~~')
            enableButtons
            APIgroup.Visible = 'off';
            TRACEgroup.Visible='off';
            axes(haxROIS); 
            pause(0.05);

        end
    
        
    else
    
        if strcmp(hObject.State,'on')
            
            delROImodeON.Value = 1;
            delROImodeOFF.Value = 0;

            memocon('ROI DELETION MODE ~~ON~~')
            disableButtons
            axes(haxAPI); 
            APIgroup.Visible='on';
            TRACEgroup.Visible='on';
            pause(0.05);

        else
            
            delROImodeON.Value = 0;
            delROImodeOFF.Value = 1;
            
            redrawROIs()

            memocon('ROI DELETION MODE ~~OFF~~')
            enableButtons
            APIgroup.Visible = 'off';
            TRACEgroup.Visible='off';
            axes(haxROIS); 
            pause(0.05);

        end
    
    end

end



%----------------------------------------------------
%        REDRAW ROIS
%----------------------------------------------------
function redrawROIs(hObject, eventdata)

    memocon('Redrawing ROIs...')
    
    delete(ROIgroup.Children)
    
    %-------------
    % REDRAW ROI CIRCLES
    %-------------
    sROI = findobj(APIgroup,'Type','patch');
    gROI = findobj(APIgroup,'Type','hggroup');
    axes(haxROIS); pause(0.05);
    
    for nn = 2:length(gROI)
    %for nn = 1:length(sROI)
        %hROI = sROI(nn);
        %ROIpos = hROI.Vertices;
        
        DispName = gROI(nn).DisplayName;
        ROIpos = gROI(nn).Children(9).Vertices;
        hROI = impoint(haxROIS, mean(ROIpos));
        ROIpos = hROI.getPosition;
        delete(hROI);
        hROI = imellipse(haxROIS, [ROIpos-round(stampSize/2) stampSize stampSize]);
        ROIpos = hROI.getPosition;
        delete(hROI);
        
        if strcmp(DispName(1),'S')
        hREC = rectangle('Parent',haxROIS,'Position',ROIpos,'Curvature',1,'LineWidth',2,'EdgeColor',[.4 .4 1]);
        else
        hREC = rectangle('Parent',haxROIS,'Position',ROIpos,'Curvature',1,'LineWidth',2,'EdgeColor',[1 1 0]);
        end
        
        set(hREC,'Parent',ROIgroup)
        
    end
    
    
    %-------------
    % REDRAW LINE TRACES
    %-------------
    gROI = findobj(TRACEgroup,'Type','hggroup');
    
    for nn = 2:length(gROI)
        
        ROIposE2 = [gROI(nn).Children(1).XData  gROI(nn).Children(1).YData];
        ROIposE1 = [gROI(nn).Children(2).XData  gROI(nn).Children(2).YData];
        hROI = imline(haxROIS,[ROIposE1; ROIposE2]);
        dpos = hROI.getPosition;
        delete(hROI);
        
        hLINE = line(dpos(:,1), dpos(:,2),'Parent',haxROIS,'Color',[0 1 0]);
        set(hLINE,'Parent',ROIgroup)
        
    end
    
	APIgroup.Visible = 'off';
	TRACEgroup.Visible='off';
    
    recountROI
        
    memocon('Done redrawing ROIs!')

end



%----------------------------------------------------
%        RECOUNT ROIS
%----------------------------------------------------
function recountROIs(hObject, eventdata)
    
    if pairedRadio.Value == 1
        ROI_IDh.String = int2str(size(APIgroup.Children,1) ./ 2);
    else
        ROI_IDh.String = int2str(size(APIgroup.Children,1));
    end
    
    memocon(sprintf('ROI Count: % s ', ROI_IDh.String));
    
end

function recountROI()
    
    if pairedRadio.Value == 1
        ROI_IDh.String = int2str(size(APIgroup.Children,1) ./ 2);
    else
        ROI_IDh.String = int2str(size(APIgroup.Children,1));
    end
    
    memocon(sprintf('ROI Count: % s ', ROI_IDh.String));
    
end




%----------------------------------------------------
%        SAVE ROI ANALYSIS DATA
%----------------------------------------------------
function saveDataset(hObject, eventdata)
    
    memocon('Preparing ROI analysis data for export...')
    
    %------
    lftthresholdMIN = str2double(lftthresholdMINh.String);
    lftthresholdMAX = str2double(lftthresholdMAXh.String);   
    intPminmax = prctile(intensity(:),...
        [str2double(intThreshMin.String) str2double(intThreshMax.String)]);
    chimin = str2double(chiminh.String);
    chimax = str2double(chimaxh.String);
    ChiGood       = (chi >= chimin & chi <= chimax);
    IntensityGood = (intensity >= intPminmax(1) & intensity <= intPminmax(2));
    LifeGood      = (lifetime >= lftthresholdMIN & lifetime <= lftthresholdMAX);
    AllGood       = (ChiGood .* IntensityGood .* LifeGood) > 0;
    %------
    
    
    
    recountROI
    Rmode = radioh.SelectedObject.String;
    doTrace = doLineTrace.Value;
    axes(haxAPI); 
    APIgroup.Visible='on'; 
    pause(0.05);
    
    
    %------
    sROI = findobj(APIgroup,'Type','patch');
    gROI = findobj(APIgroup.Children,'Type','hggroup');
    
    % haxROIS.Children(1).Children(1).Color
    % haxROIS.Children(1).DisplayName
    % sROI(1).Parent.DisplayName
    
    
    for nn = 1:length(sROI)
        
        ROInameID = sROI(nn).Parent.DisplayName;
        
        sROIpos = sROI(nn).Vertices;
        sROIarea = polyarea(sROIpos(:,1),sROIpos(:,2));
        sROImask = poly2mask(sROIpos(:,1),sROIpos(:,2), ...
                             size(intensity,1), size(intensity,2));


        ROI_LIFETIME  = lifetime .* AllGood .* sROImask;
        ROI_INTENSITY = intensity .* AllGood .* sROImask;
        ROI_CHI       = chi .* AllGood .* sROImask;

        ROI_imgG      = imgG .* AllGood .* sROImask;
        ROI_imgR      = imgR .* AllGood .* sROImask;



        ROI_LIFETIME_MEAN  = mean(ROI_LIFETIME(ROI_LIFETIME > 0));
        ROI_INTENSITY_MEAN = mean(ROI_INTENSITY(ROI_INTENSITY > 0));
        ROI_CHI_MEAN       = mean(ROI_CHI(ROI_CHI > 0));
        ROI_imgG_MEAN      = mean(ROI_imgG(ROI_imgG > 0))*1000;
        ROI_imgR_MEAN      = mean(ROI_imgR(ROI_imgR > 0))*1000;



        flimdats{nn} = {ROI_LIFETIME_MEAN, ...
                        ROI_INTENSITY_MEAN, ...
                        ROI_CHI_MEAN, ...
                        sROIarea, ...
                        ROI_imgG_MEAN, ...
                        ROI_imgR_MEAN, ...
                        ROInameID};
    
    
    end
    % ------
    
    
    
    
    cd(datfilefdir);

    saveDatafilename = inputdlg('Enter a filename to save data','file name',1,...
                                {datfilef(1:end-4)});

    Datafilename = char(strcat(saveDatafilename));
    
    for nn = 1:size(flimdats,2)

        % flimdat(nn,:) = [flimdats{1,nn}{1:6} maglevel dendritesize]; 
        flimdat(nn,:) = [flimdats{1,nn}{1:6}];
        
        ROInames{nn} = flimdats{1,nn}{7};
        if strcmp(ROInames{nn},'')
            ROInames{nn} = num2str(nn);
        end
    end
    
    
    flimtable = array2table(flimdat);
    flimtable.Properties.VariableNames = {'LIFETIME' 'INTENSITY' 'CHISQR' 'VOLUME' 'GFP' 'RFP'};
    flimtable.Properties.RowNames = ROInames;
    
    writetable(flimtable,[Datafilename '_ROIANALYSIS.csv'],'WriteRowNames',true)
    memocon('ROI analysis data saved successfully!')
    
    
    
    
    %------------------------------------------------------
    memocon(' ')
    memocon('Preparing LINE profile data for export...')
    
    
    sROI = findobj(TRACEgroup,'Type','Line','-regexp','Tag','end \w+');
    gROI = findobj(TRACEgroup.Children,'Type','hggroup');
    %lineROIS = findobj(haxROIS,'-regexp','DisplayName','L\w+');
    
    linedat = zeros(length(gROI)*3,100);
    linedatnames = cell(length(gROI)*3,1);
    mm=1;
    if length(gROI) < 1
        memocon('No LINE profile data found.')
    else
    for nn = 1:length(gROI)
        
        L_ROInameID = gROI(nn).DisplayName;
        
        Lpos = [gROI(nn).Children(2).XData  gROI(nn).Children(2).YData;...
                gROI(nn).Children(1).XData  gROI(nn).Children(1).YData];
            
        L_spineextent = sqrt((Lpos(1,1)-Lpos(2,1))^2 + (Lpos(1,2)-Lpos(2,2))^2);
        
        [SPIx,SPIy,SPIf] = improfile(imgF, Lpos(:,1), Lpos(:,2), round(L_spineextent));
        [SPIx,SPIy,SPIg] = improfile(imgG, Lpos(:,1), Lpos(:,2), round(L_spineextent));
        [SPIx,SPIy,SPIr] = improfile(imgR, Lpos(:,1), Lpos(:,2), round(L_spineextent));
        
        linedat(mm,1:numel(SPIf)) = SPIf(:); linedatnames{mm} = [L_ROInameID '_FLI']; mm=mm+1;
        linedat(mm,1:numel(SPIg)) = SPIg(:); linedatnames{mm} = [L_ROInameID '_GRN']; mm=mm+1;
        linedat(mm,1:numel(SPIr)) = SPIr(:); linedatnames{mm} = [L_ROInameID '_RED']; mm=mm+1;
        % linedatnames{mm} = ['___' mm '___']; mm=mm+1;
        
    end
    
    linetable = array2table(linedat);
    linetable.Properties.RowNames = linedatnames;
    
    writetable(linetable,[Datafilename '_LINEPROFILES.csv'],'WriteRowNames',true)
    memocon('Line profiles saved successfully!')
    end
    
    memocon('DATA EXPORT FINISHED!')

end






%----------------------------------------------------
%        EXPORT ROI LOCATIONS
%----------------------------------------------------
function saveROIfun(hObject, eventdata)
    
    memocon('Preparing to export ROI location data.')
    memocon('Enter a name for the ROI export file.')
    memocon(' ''_ROIs.mat'' will be appended to this filename.')
    
    cd(datfilefdir);

    ROIfna = inputdlg('Name of ROI export file.','file name',1,{datfilef(1:end-4)});
    ROIfname = char(strcat(ROIfna));
    
    save([ROIfname '_ROIs.mat'],'APIgroup', 'TRACEgroup', 'ROIgroup')
    
    memocon([ROIfname '_ROIs.mat  was created.'])
    memocon('ROI LOCATION DATA EXPORTED SUCCESSFULLY!')

end




%----------------------------------------------------
%        LOAD ROI DATA
%----------------------------------------------------
function loadROIfun(hObject, eventdata)

    
    recountROI
    axes(haxROIS)
    memocon('Select .csv or .mat file with ROI data')
    
    
    [ROIFileName,ROIPathName,ROIFilterIndex] = uigetfile({'*.csv','*.mat'});
    ROIfullpath = [ROIPathName ROIFileName];
    
    
    if strcmp(ROIFileName(end-2:end),'csv')
        
        ROIcsv = csvread(ROIfullpath);
        
        ROIcsv(1,:) = [];
        ROIcsv(:,1) = [];
        ROIcsv(:,1:2:31) = [];
        ROIcsv(ROIcsv(:,1) == 0 & ROIcsv(:,2) == 0,:) = [];
        
        SVx = [ROIcsv(:,1), ROIcsv(:,3), ROIcsv(:,5), ROIcsv(:,7), ROIcsv(:,1), NaN(size(ROIcsv,1),1)];
        SVy = [ROIcsv(:,2) ROIcsv(:,4) ROIcsv(:,6) ROIcsv(:,8) ROIcsv(:,2) NaN(size(ROIcsv,1),1)];
        DVx = [ROIcsv(:,9) ROIcsv(:,11) ROIcsv(:,13) ROIcsv(:,15) ROIcsv(:,9) NaN(size(ROIcsv,1),1)];
        DVy = [ROIcsv(:,10) ROIcsv(:,12) ROIcsv(:,14) ROIcsv(:,16) ROIcsv(:,10) NaN(size(ROIcsv,1),1)];
        SVx = SVx';
        SVy = SVy';
        DVx = DVx';
        DVy = DVy';
        
        
        for nn = 1:size(SVx,2)
            
            
            %--- PROCESS SPINE ROI IMPORTS ---
            ROI_IDh.String = int2str(str2num(ROI_IDh.String) + 1);
            
            axes(haxAPI); APIgroup.Visible='on'; pause(0.02);
            
            SPI_hROI = imellipse(haxAPI, [[mean(SVx(1:4,nn)),mean(SVy(1:4,nn))]-round(stampSize/2) stampSize stampSize]);
            SPIpos = SPI_hROI.getPosition;
            setColor(SPI_hROI,[.4 .4 1]);
            
            haxAPI.Children(1).DisplayName = ['S' ROI_IDh.String];
            set(SPI_hROI,'Parent',APIgroup)
            
            
            APIgroup.Visible = 'off';
            axes(haxROIS); pause(0.02);
        
            hREC = rectangle('Parent',haxROIS,'Position',SPIpos,'Curvature',1,'LineWidth',2,'EdgeColor',[.4 .4 1]);
            set(hREC,'Parent',ROIgroup)
            
            
            
            %--- PROCESS DENDRITE ROI IMPORTS ---
            axes(haxAPI); APIgroup.Visible='on'; pause(0.05);
            
            DEN_hROI = imellipse(haxAPI, [[mean(DVx(1:4,nn)),mean(DVy(1:4,nn))]-round(stampSize/2) stampSize stampSize]);
            DENpos = DEN_hROI.getPosition;
            setColor(DEN_hROI,[1 1 0]);
            
            haxAPI.Children(1).DisplayName = ['D' ROI_IDh.String];
            set(DEN_hROI,'Parent',APIgroup)
            
            
            APIgroup.Visible = 'off';
            axes(haxROIS); pause(0.02);
        
            hREC = rectangle('Parent',haxROIS,'Position',DENpos,'Curvature',1,'LineWidth',2,'EdgeColor',[1 1 0]);
            set(hREC,'Parent',ROIgroup)
                        
            
            % nROI = numel(ROI)+1;
            % ROI(nROI).POS     = SPIpos;
            % ROI(nROI).POSd    = DENpos;
            % ROI(nROI).PH      = SPI_hROI;
            % ROI(nROI).PHd     = DEN_hROI;
            % ROI(nROI).TYPE    = 'CSVIMPORT';

            pause(.02)
        
        end

        memocon('Finished importing ROIs.')     
        
    elseif strcmp(ROIFileName(end-2:end),'mat')
        memocon('The code to import mat files is not finished.')
        memocon('Brad will finish this in the next update.')
    else
        memocon('ROI import file must be .csv or .mat')
        return
    end

    

end















%----------------------------------------------------
%        FORMAT IMAGE
%----------------------------------------------------
function [IMG] = imformat(imgpath)
%% imformat.m

% filedir = uigetdir();
% imgFiles = dir([filedir,'/*Y.jpg']); 
% imgFiles = dir([filedir,'/*Y.*']);
% imgFileNames = {imgFiles(:).name}';
% imgpath = [filedir, '/', imgFileNames{1}];

if ischar(imgpath) || iscellstr(imgpath)

    iminfo = imfinfo(imgpath);
    [im, map] = imread(imgpath);
    im_ctype = iminfo.ColorType;

else
    im = imgpath;
    im_ctype = 'truecolor';
end


im_size = size(im);
im_nmap = numel(map);



if strcmp(im_ctype, 'truecolor') || numel(im_size) > 2

    IMG = rgb2gray(im);
    IMG = im2double(IMG);

elseif strcmp(im_ctype, 'indexed')

    IMG = ind2gray(im,map);
    IMG = im2double(IMG);

elseif strcmp(im_ctype, 'grayscale')

    IMG = im2double(im);

else

    IMG = im;

end


end




%----------------------------------------------------
%        LINTRANS LINE TRACE PROFILE
%----------------------------------------------------
function [y] = lintrans(x,a,b,c,d)

    %lintrans = @(x,a,b,c,d) (c*(1-(x-a)/(b-a)) + d*((x-a)/(b-a)));    
    y = (c.*(1-(x-a)./(b-a)) + d.*((x-a)./(b-a)));
    
end    




%----------------------------------------------------
%        CSUS DROPDOWN MENU CALLBACK
%----------------------------------------------------
function getStampSize(hObject, eventdata)
    
    PopValue = stampSizeH.Value;
    stampSize = str2num(stampSizeH.String{PopValue});

    % set(stampSizeH, 'String', stampvals);
    
    memocon(sprintf('Stamp size set to: % s ', num2str(stampSize)));
    % memocon(sprintf('stampSizeH.Value: % 0.1g ', stampSizeH.Value));
    % memocon(sprintf('stampSizeH.ListboxTop: % 0.1g ', stampSizeH.ListboxTop));
    
end






%----------------------------------------------------
%        MEMO CONSOLE
%----------------------------------------------------
function memocon(spf,varargin)
    
  
    if iscellstr(spf)
        spf = [spf{:}];
    end
    
    if iscell(spf)
        return
        keyboard
        spf = [spf{:}];
    end
    
    if ~ischar(spf)
        return
        keyboard
        spf = [spf{:}];
    end
    
    

    memes(1:end-1) = memes(2:end);
    memes{end} = spf;
    conboxH.String = memes;
    pause(.02)
    
    if nargin == 3
        
        vrs = deal(varargin);
                
        memi = memes;
                 
        memes(1:end) = {' '};
        memes{end-1} = vrs{1};
        memes{end} = spf;
        conboxH.String = memes;
        
        conboxH.FontAngle = 'italic';
        conboxH.ForegroundColor = [.9 .4 .01];
        pause(vrs{2})
        
        conboxH.FontAngle = 'normal';
        conboxH.ForegroundColor = [0 0 0];
        conboxH.String = memi;
        pause(.02)
        
    elseif nargin == 2
        vrs = deal(varargin);
        pause(vrs{1})
    end
    
    
    disp(spf)

end





%----------------------------------------------------
%        ENABLE AND DISABLE GUI BUTTONS
%----------------------------------------------------
function enableButtons()

    changeimgh.Enable = 'on';
    getROIh.Enable = 'on';
    loadROISh.Enable = 'on';
    exploreFLIMh.Enable = 'on';
    saveROISh.Enable = 'on';
    savefileh.Enable = 'on';
    loadIMGh.Enable = 'on';
    setintenh.Enable = 'on';
    dftintenh.Enable = 'on';
    
end
function disableButtons()
    
    changeimgh.Enable = 'off';
    getROIh.Enable = 'off';
    loadROISh.Enable = 'off';
    exploreFLIMh.Enable = 'off';
    saveROISh.Enable = 'off';
    savefileh.Enable = 'off';
    loadIMGh.Enable = 'off';
    setintenh.Enable = 'off';
    dftintenh.Enable = 'off';

end




end
%% EOF
