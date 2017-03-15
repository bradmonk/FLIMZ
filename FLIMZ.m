function varargout = FLIMZ(varargin)
%% FLIMIM.m USAGE NOTES
%{

Syntax
-----------------------------------------------------
    FLIMCCD()
    FLIMCCD(datfilefdir)


Description
-----------------------------------------------------
    FLIMCCD() can be run with no arguments passed in. In this case user
    will be prompted to select a directory which contains the FLIM dat 
    file along with the corresponding CCD images. Optionally this function can 
    be called using FLIMCCD(datfilefdir) where the full path to the data directory
    is explicitly provided.
    

Useage Definitions
-----------------------------------------------------


    FLIMCCD()
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

thisfile = 'FLIMIM.m';
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


global haxROIS haxTABH
global LifeImageFile FLIMcmap
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
global ROIcsv CLIMslider

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
haxes = {};
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
    'Position', [0.05 0.15 0.8 0.8], 'OuterPosition', [-.2 0 1 1],...
    'PlotBoxAspectRatio', [1 1 1],'XColor','none','YColor','none'); 
    % ,'XDir','reverse',...

haxCCDR = axes('Parent', MAINGUIFIG, 'NextPlot', 'Add',...
    'Position', [0.05 0.15 0.8 0.8], 'OuterPosition', [-.2 0 1 1],...
    'PlotBoxAspectRatio', [1 1 1],'XColor','none','YColor','none'); 
    % ,'XDir','reverse',...

haxFLIM = axes('Parent', MAINGUIFIG, 'NextPlot', 'Add',...
    'Position', [0.05 0.15 0.8 0.8], 'OuterPosition', [-.2 0 1 1],...
    'PlotBoxAspectRatio', [1 1 1],'XColor','none','YColor','none');

haxROIS = axes('Parent', MAINGUIFIG, 'NextPlot', 'Add','Color','none',...
    'Position', [0.05 0.15 0.8 0.8], 'OuterPosition', [-.2 0 1 1],...
    'PlotBoxAspectRatio', [1 1 1],'XColor','none','YColor','none');

linkaxes([haxFLIM,haxCCDG,haxCCDR,haxROIS],'xy')
haxes = {haxFLIM, haxCCDG, haxCCDR, haxROIS};


changeimgh = uicontrol('Parent', MAINGUIFIG, 'Units', 'normalized', ...
    'Position', [0.10 0.03 0.1 0.04], 'String', 'Next Image', 'FontSize', 11,...
    'Callback', @changeimg);


CLIMsliderTxt = uicontrol('Parent', MAINGUIFIG, 'Style', 'Text', 'Units', 'normalized', ...
    'Position', [.01 .97 .09 .02], 'FontSize', 10, 'String', 'IMG Intensity','BackgroundColor',[1 1 1]);
CLIMslider = uicontrol('Parent', MAINGUIFIG, 'Units', 'normalized','Style','slider',...
	'Max',150,'Min',1,'Value',50,'SliderStep',[0.01 0.10],...
	'Position', [.10 .97 .30 .02], 'Callback', @setclim);




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
    'Position', [.78 .87 .20 .11], 'FontSize', 10, 'String', 'Stamp Size','BackgroundColor',[1 1 1]);


stampvals = {'9','10','11','12','14','16','18','20'};
stampSizeH = uicontrol('Parent', IPpanelH,'Style', 'popup',...
    'Units', 'normalized', 'String',stampvals,...
    'Position', [.78 .78 .20 .12],'BackgroundColor',[.9 .9 .9],...
    'Callback', @getStampSize, 'Enable','on');
stampSizeH.Value = 3;
stampSizeH.ListboxTop = 3;


%---------------------
boxtypeh = uibuttongroup('Parent', IPpanelH, 'Visible','off',...
                  'Units', 'normalized','BackgroundColor',[1 1 1],...
                  'Position',[.02 .50 .95 .24],...
                  'SelectionChangedFcn',@boxselection);
              
% Create three radio buttons in the button group.
boxtypeh1 = uicontrol(boxtypeh,'Style','radiobutton','FontSize', 11,...
                  'String','freehand','BackgroundColor',[1 1 1],...
                  'Units', 'normalized',...
                  'Position',[0.05 0.05 0.3 0.9],...
                  'HandleVisibility','off');
              
boxtypeh2 = uicontrol(boxtypeh,'Style','radiobutton','FontSize', 11,...
                  'String','rectangle','BackgroundColor',[1 1 1],...
                  'Units', 'normalized',...
                  'Position',[0.3 0.05 0.3 0.9],...
                  'HandleVisibility','off');

boxtypeh3 = uicontrol(boxtypeh,'Style','radiobutton','FontSize', 11,...
                  'String','elipse','BackgroundColor',[1 1 1],...
                  'Units', 'normalized',...
                  'Position',[0.55 0.05 0.3 0.9],...
                  'HandleVisibility','off');

boxtypeh4 = uicontrol(boxtypeh,'Style','radiobutton','FontSize', 11,...
                  'String','stamp','BackgroundColor',[1 1 1],...
                  'Units', 'normalized',...
                  'Position',[0.77 0.05 0.3 0.9],...
                  'HandleVisibility','off');              
boxtypeh.Visible = 'on';
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
    
    boxtypeh.SelectedObject = boxtypeh4; % Set radiobutton to stamp
    
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

    axes(haxROIS)
    ROI_IDh.String = int2str(str2num(ROI_IDh.String) + 1);
    
    %---------------------------
    % GET SPINE ROI
    %---------------------------
    if strcmp(boxtypeh.SelectedObject.String,'rectangle')
        
        hROI = imrect(haxROIS);
        ROIpos = hROI.getPosition;
        ROIarea = ROIpos(3) * ROIpos(4);
        setColor(hROI,[0 0 1]);
        
    elseif strcmp(boxtypeh.SelectedObject.String,'elipse')
        
        hROI = imellipse(haxROIS);
        ROIpos = hROI.getPosition;
        ROIarea = pi * (.5*ROIpos(3)) * (.5*ROIpos(4));
        setColor(hROI,[0 0 1]);
        
    elseif strcmp(boxtypeh.SelectedObject.String,'stamp')
        
        % [x,y] = FLIMginput(2,'custom');
        hROI = impoint;
        ROIpos = hROI.getPosition;
        delete(hROI)
        hROI = imellipse(haxROIS, [ROIpos-round(stampSize/2) stampSize stampSize]);
        ROIpos = hROI.getPosition;
        ROIarea = pi * (stampSize/2)^2;
        setColor(hROI,[0 0 1]);
        
    else % strcmp(boxtypeh.SelectedObject.String,'freehand')
        hROI = imfreehand(haxROIS);
        ROIpos = hROI.getPosition;
        ROIarea = polyarea(ROIpos(:,1),ROIpos(:,2));
        setColor(hROI,[0 0 1]);
    end
    
    haxROIS.Children(1).DisplayName = ['S' ROI_IDh.String];
    
    nROI = numel(ROI)+1;
    ROI(nROI).POS   = ROIpos;
    ROI(nROI).PH    = hROI;
    ROI(nROI).TYPE  = boxtypeh.SelectedObject.String;



    %---------------------------
    % GET DENDRITE ROI
    %---------------------------
    Q_getd = questdlg('Select dendrite ROI?', 'Select dendrite ROI?', 'Yes', 'No', 'Yes');
    if strcmp(Q_getd,'Yes')
        axes(haxROIS)
        
        if strcmp(boxtypeh.SelectedObject.String,'rectangle')

            hROI = imrect(haxROIS);
            ROIpos = hROI.getPosition;
            ROIarea = ROIpos(3) * ROIpos(4);
            setColor(hROI,[1 1 0]);

        elseif strcmp(boxtypeh.SelectedObject.String,'elipse')

            hROI = imellipse(haxROIS);
            ROIpos = hROI.getPosition;
            ROIarea = pi * (.5*ROIpos(3)) * (.5*ROIpos(4));
            setColor(hROI,[1 1 0]);

        elseif strcmp(boxtypeh.SelectedObject.String,'stamp')

            % [x,y] = FLIMginput(2,'custom');
            hROI = impoint;
            ROIpos = hROI.getPosition;
            delete(hROI)
            hROI = imellipse(haxROIS, [ROIpos-round(stampSize/2) stampSize stampSize]);
            ROIpos = hROI.getPosition;
            ROIarea = pi * (stampSize/2)^2;
            setColor(hROI,[1 1 0]);

        else % strcmp(boxtypeh.SelectedObject.String,'freehand')
            hROI = imfreehand(haxROIS);
            ROIpos = hROI.getPosition;
            ROIarea = polyarea(ROIpos(:,1),ROIpos(:,2));
            setColor(hROI,[1 1 0]);
        end
        
        haxROIS.Children(1).DisplayName = ['D' ROI_IDh.String];
        
    end
    
    
    %---------------------------
    % GET LINE-TRACE ROI
    %---------------------------
    Q_getline = questdlg('Draw ROI Line?', 'Draw ROI Line?', 'Yes', 'No', 'Yes');
    if strcmp(Q_getline,'Yes')
        axes(haxROIS)

        hROI = imline(haxROIS);
        dpos = hROI.getPosition;
        setColor(hROI,[0 1 0]);
        haxROIS.Children(1).DisplayName = ['L' ROI_IDh.String];
        
        spineextent = sqrt((dpos(1,1)-dpos(2,1))^2 + (dpos(1,2)-dpos(2,2))^2);
        
        [SPIx,SPIy,SPIf] = improfile(imgF, dpos(:,1), dpos(:,2), round(spineextent));
        [SPIx,SPIy,SPIg] = improfile(imgG, dpos(:,1), dpos(:,2), round(spineextent));
        [SPIx,SPIy,SPIr] = improfile(imgR, dpos(:,1), dpos(:,2), round(spineextent));
        
        % scatter(haxTABH, cx,cy, dotsz,'MarkerFaceColor', [1 0 0])
        axes(haxTABH); hold off;
        plot(haxTABH, SPIf); hold on;
        plot(haxTABH, SPIg); hold on;
        plot(haxTABH, SPIr)
        axes(haxROIS)
        
    end



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
function setclim(hObject, eventdata)
% Hints: hObject.Value returns position of slider
%        hObject.Min and hObject.Max determine range of slider
% sunel = get(handles.CLIMslider,'value'); % Get current light elev.
% sunaz = get(hObject,'value');   % Varies from -180 -> 0 deg


lowerinten = str2num(intThreshMin.String);
upperinten = str2num(intThreshMax.String);
lowerintenPCT = prctile(intensity(:),lowerinten);
upperintenPCT = prctile(intensity(:),upperinten);


slideVal = ceil(CLIMslider.Value);

haxFLIM.CLim = [(lowerintenPCT / ((slideVal/50)*(slideVal/50))),...
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
function boxselection(source,callbackdata)
    
    % callbackdata.OldValue.String
    % boxtypeh.SelectedObject.String
    % boxtype = callbackdata.NewValue.String;
    
    memocon(['ROI Selector: ' callbackdata.NewValue.String])
    
    display(['Previous: ' callbackdata.OldValue.String]);
    display(['Current: ' callbackdata.NewValue.String]);
    display('------------------');

end






%----------------------------------------------------
%           GET DENDRITE SIZE
%----------------------------------------------------
function getdendsize(boxidselecth, eventdata)


   hline = imline;
   dpos = hline.getPosition();
    
   dendritesize = sqrt((dpos(1,1)-dpos(2,1))^2 + (dpos(1,2)-dpos(2,2))^2);

   disp(['dendrite size:' num2str(dendritesize)])

end






%----------------------------------------------------
%           COMPILE DATA FILE
%----------------------------------------------------
function compilefile(hObject, eventdata)
%Compile function triggers the user input for three datafiles. (Will return
%error if files are not of the right size, etc.) Compile then writes a new 
%file that is the stacked version of all three files. 

    set(initmenuh, 'Visible', 'Off');

    home = cd;
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

                                                        disp(strcat(savefilename, ' was successfully compiled.'));
                                                    else
                                                        savefilename = mat2str(strcat(chifilename, '.dat'));
                                                        savefilename = savefilename(2:end-1);
                                                        disp('Error compiling ', savefilename, '. One or more of the component files may be incorrect');

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
    msgbox('All files successfully compiled.');



    cd(home);

    set(initmenuh, 'Visible', 'On');


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
%        LOAD FILE
%----------------------------------------------------
function loadfile(hObject, eventdata)
%Load file triggers uiresume; the initial menu is set to invisible. Prompts
%user for file to load, copies the datastack from the file; sets the image 
%windows to visible, and plots the images.    

    set(initmenuh, 'Visible', 'Off');

    
    if numel(datfilef) < 1
        [datfilef, datfilefdir] = uigetfile('*.dat', 'load a .dat file');
    end
    
    if numel(ccdfileg) < 1
        [ccdfileg, ccdfilegdir] = uigetfile({'*.tif*';'*.bmp';'*.jpg';'*.png'}, 'load green-channel .bmp file');
    end

    if numel(ccdfiler) < 1
        [ccdfiler, ccdfilerdir] = uigetfile({'*.tif*';'*.bmp';'*.jpg';'*.png'}, 'load red-channel .dat file');
    end
    
    
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
    
    
    imgG = ccdget([ccdfilegdir ccdfileg]);
    imgR = ccdget([ccdfilerdir ccdfiler]);
    
    
    
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
        
        
        
    
    
%{
FLIMsets{1} = 'SPCI2/';
basedir = '/Users/bradleymonk/Documents/MATLAB/myToolbox/LAB/FRET_FLIM/FRETdata/ActinProfilin/2016_02_21/';

%---------------
pixdir = FLIMsets{1};
filedir = [basedir pixdir];
regexpStr = '((\S)+(\.tif|\.jpg|\.bmp|\.png+))';
allfileinfo = dir(filedir);
allfilenames = {allfileinfo.name};
filepaths    = fullfile(filedir,allfilenames(~cellfun('isempty',regexp(allfilenames,regexpStr))));
fprintf('%s \r',filepaths{:})
%---------------


% IMGfiles = {IMGfilenames.name};
for nn = 1:numel(filepaths)
    
    IMGs{nn} = imgcon(filepaths{nn});

end

IMG = IMGs{1};
    
%}
    
    

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

    set(magh, 'String', num2str(magnification));

    set(MAINGUIFIG, 'Name', datfilef);
    set(ROI_IDh, 'String', int2str(1));
    set(haxCCDG, 'XLim', [1 size(imgG,2)], 'YLim', [1 size(imgG,1)]);
    % set(haxCCDG, 'YLim', [1 ydim]);
    set(haxCCDR, 'XLim', [1 size(imgR,2)], 'YLim', [1 size(imgR,1)]);
    % set(haxCCDR, 'YLim', [1 ydim]);
    set(haxFLIM, 'XLim', [1 xdim]);
    set(haxFLIM, 'YLim', [1 ydim]);
    
    
    set(stampSizeH, 'String', num2str(stampSize));
    boxtypeh.SelectedObject = boxtypeh4; % Set radiobutton to stamp
    % boxtype = boxtypeh.SelectedObject.String;
    %----------------------------------------------------
    
    
    
    
    %----------------------------------------------------
    %                   DRAW IMAGE
    %----------------------------------------------------
    
    axes(haxCCDG)
    colormap(haxCCDG,hot)
    phCCDG = imagesc(imgG , 'Parent', haxCCDG);
              pause(1)
              
    axes(haxCCDR)
    colormap(haxCCDR,hot)
    phCCDR = imagesc(imgR, 'Parent', haxCCDR);
        pause(1)
        
    axes(haxFLIM)    
    colormap(haxFLIM,hot)
    phFLIM = imagesc(intensity, 'Parent', haxFLIM,...
                  [prctile(intensity(:),intenseThreshMIN) prctile(intensity(:),intenseThreshMAX)]);
        pause(1)
        axes(haxFLIM)
    
    pause(.2)
    imXlim = haxFLIM.XLim;
    imYlim = haxFLIM.YLim;
    


end





%----------------------------------------------------
%        LOAD FILE
%----------------------------------------------------
function saveDataset(hObject, eventdata)
    
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
    
    
    
    
    
    %------
    sROI = findobj(haxROIS,'Type','patch');
    
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
    memocon('ROI ANALYSIS data saved successfully!')
    
    
    
    
    
    memocon(' ')
    memocon('Now saving line profiles...')
    
    lineROIS = findobj(haxROIS,'-regexp','DisplayName','L\w+');
    
    linedat = zeros(length(lineROIS)*3,100);
    linedatnames = cell(length(lineROIS)*3,1);
    mm=1;
    for nn = 1:length(lineROIS)
        
        L_ROInameID = lineROIS(nn).DisplayName;
        
        Lpos = [lineROIS(nn).Children(2).XData  lineROIS(nn).Children(2).YData;...
                lineROIS(nn).Children(1).XData  lineROIS(nn).Children(1).YData];
            
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






%----------------------------------------------------
%        PREP FOR SAVE
%----------------------------------------------------
function prepForSave(savefileh, eventData)
    
    % ------
    
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

    
    
    
    
    % ------
    sROI = findobj(haxFLIM,'Type','patch');
    
    for nn = 1:length(sROI)
        
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
                        ROI_imgR_MEAN};
    
    
    end
    % ------
        
end







%----------------------------------------------------
%        SAVE DATA
%----------------------------------------------------
function saveFile(savefileh, eventData)
    
    
    prepForSave(savefileh, eventData)
    

    cd(datfilefdir);

    saveDatafilename = inputdlg('Enter a filename to save data','file name',1,...
                                {datfilef(1:end-4)});

    Datafilename = char(strcat(saveDatafilename));

    maglevel = str2double(magh.String);
    
    if numel(dpos) < 1; % If dendrite size was manually selected, numel(dpos) > 1
        dendritesize = maglevel*5;
    end
    
    
    
    for nn = 1:size(flimdats,2)
        
        VxD = flimdats{1,nn}{4} ./ (.5 .* dendritesize).^2;
        
        flimdat(nn,:) = [flimdats{1,nn}{1:6} maglevel dendritesize VxD];        
        ROInames{nn} = num2str(nn);        
    end
    
    
    
    flimtab = array2table(flimdat);
    flimtab.Properties.VariableNames = {'LIFETIME' 'INTENSITY' 'CHISQR' 'VOLUME' ...
                                        'GFP' 'RFP' 'MAG' 'DSIZE' 'VxD'};
    flimtab.Properties.RowNames = ROInames;
    
    
    
    writetable(flimtab,[Datafilename '.csv'],'WriteRowNames',true)
    disp('Data saved successfully!')
    % msgbox('Data saved successfully');

    cd(home);


end





%----------------------------------------------------
%        SAVE ROI DATA
%----------------------------------------------------
function saveROIfun(hObject, eventdata)

    
%{
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
        
        for nn = 1:size(ROIcsv,1)
            Sh{nn} = line('XData',SVx(:),'YData',SVy(:), 'Color', 'r', 'Parent', gca);
            Sh{nn}.DisplayName = num2str(nn); Sh{nn}.Tag = num2str(nn);
            Dh{nn} = line('XData',DVx(:),'YData',DVy(:), 'Color', 'y', 'Parent', gca);
            Dh{nn}.DisplayName = num2str(nn); Dh{nn}.Tag = num2str(nn);
        end
        
        % processLoadedROI(SVx, SVy, DVx, DVy)
        
        
    elseif strcmp(ROIFileName(end-2:end),'mat')
        ROIloaded = load([ROIPathName ROIFileName],'MORPHdata','ROISAVES');

        MORPHdata = ROIloaded.MORPHdata;
        ROISAVES = ROIloaded.ROISAVES;

        memocon('ROI data loaded into workspace!')

        lwd = 4;
        axes(haxCCD)
        for nn = 1:length(ROISAVES)
            line(ROISAVES(nn).SpinePos(:,1), ROISAVES(nn).SpinePos(:,2),'Color',[.95 .95 .10],'LineWidth',lwd)
            line(ROISAVES(nn).SpineExtentPos(:,1), ROISAVES(nn).SpineExtentPos(:,2),'Color',[.10 .95 .95],'LineWidth',lwd)
            line(ROISAVES(nn).SpineHeadPos(:,1), ROISAVES(nn).SpineHeadPos(:,2),'Color',[.95 .10 .95],'LineWidth',lwd)
            line(ROISAVES(nn).SpineNeckPos(:,1), ROISAVES(nn).SpineNeckPos(:,2),'Color',[.95 .10 .10],'LineWidth',lwd)
            line(ROISAVES(nn).DendritePos(:,1), ROISAVES(nn).DendritePos(:,2),'Color',[.10 .95 .10],'LineWidth',lwd)
            line(ROISAVES(nn).SpineNearPos(:,1), ROISAVES(nn).SpineNearPos(:,2),'Color',[.10 .10 .95],'LineWidth',lwd)
        end
    else
        memocon('ROI import file must be .csv or .mat')
        return
    end

    
%}
end




%----------------------------------------------------
%        LOAD ROI DATA
%----------------------------------------------------
function loadROIfun(hObject, eventdata)

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
            
            SPI_hROI = imellipse(haxROIS, [[mean(SVx(1:4,nn)),mean(SVy(1:4,nn))]-round(stampSize/2) stampSize stampSize]);
            SPIpos = SPI_hROI.getPosition;
            setColor(SPI_hROI,[0 0 1]);
            
            
            DEN_hROI = imellipse(haxROIS, [[mean(DVx(1:4,nn)),mean(DVy(1:4,nn))]-round(stampSize/2) stampSize stampSize]);
            DENpos = DEN_hROI.getPosition;
            setColor(DEN_hROI,[1 1 0]);
            
            
            %haxROIS.Children(1).Children(1).Color
            haxROIS.Children(1).DisplayName = ['D' num2str(nn)];
            haxROIS.Children(2).DisplayName = ['S' num2str(nn)];

            
            nROI = numel(ROI)+1;
            ROI(nROI).POS     = SPIpos;
            ROI(nROI).POSd    = DENpos;
            ROI(nROI).PH      = SPI_hROI;
            ROI(nROI).PHd     = DEN_hROI;
            ROI(nROI).TYPE    = 'CSVIMPORT';

            pause(.05)
        
        end

        ROI_IDh.String = int2str(str2num(ROI_IDh.String) + size(SVx,2));
        memocon('Finished importing ROIs.')
        
        % keyboard
        % for nn = 1:size(SVx,2)
        %     disp(haxROIS.Children(nn).DisplayName)
        % end

%         for nn = 1:size(ROIcsv,1)
%             Sh{nn} = line('XData',SVx(:),'YData',SVy(:), 'Color', 'r', 'Parent', gca);
%             Sh{nn}.DisplayName = num2str(nn); Sh{nn}.Tag = num2str(nn);
%             Dh{nn} = line('XData',DVx(:),'YData',DVy(:), 'Color', 'y', 'Parent', gca);
%             Dh{nn}.DisplayName = num2str(nn); Dh{nn}.Tag = num2str(nn);
%         end
        
        % processLoadedROI(SVx, SVy, DVx, DVy)
        
        
    elseif strcmp(ROIFileName(end-2:end),'mat')
        ROIloaded = load([ROIPathName ROIFileName],'MORPHdata','ROISAVES');

        MORPHdata = ROIloaded.MORPHdata;
        ROISAVES = ROIloaded.ROISAVES;

        memocon('ROI data loaded into workspace!')

        lwd = 4;
        axes(haxROIS)
        for nn = 1:length(ROISAVES)
            line(ROISAVES(nn).SpinePos(:,1), ROISAVES(nn).SpinePos(:,2),'Color',[.95 .95 .10],'LineWidth',lwd)
            line(ROISAVES(nn).SpineExtentPos(:,1), ROISAVES(nn).SpineExtentPos(:,2),'Color',[.10 .95 .95],'LineWidth',lwd)
            line(ROISAVES(nn).SpineHeadPos(:,1), ROISAVES(nn).SpineHeadPos(:,2),'Color',[.95 .10 .95],'LineWidth',lwd)
            line(ROISAVES(nn).SpineNeckPos(:,1), ROISAVES(nn).SpineNeckPos(:,2),'Color',[.95 .10 .10],'LineWidth',lwd)
            line(ROISAVES(nn).DendritePos(:,1), ROISAVES(nn).DendritePos(:,2),'Color',[.10 .95 .10],'LineWidth',lwd)
            line(ROISAVES(nn).SpineNearPos(:,1), ROISAVES(nn).SpineNearPos(:,2),'Color',[.10 .10 .95],'LineWidth',lwd)
        end
    else
        memocon('ROI import file must be .csv or .mat')
        return
    end

    

end






%----------------------------------------------------
%        CSUS DROPDOWN MENU CALLBACK
%----------------------------------------------------
function recountROIs(hObject, eventdata)
    
    ROI_IDh.String = int2str(size(haxROIS.Children,1));
    
    memocon(sprintf('ROI Count: % s ', ROI_IDh.String));
    
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
