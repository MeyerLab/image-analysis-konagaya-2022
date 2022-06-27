%% Rename IX micro images
IX_Micro_rename('H:\Nalin\201901003-C198-live\10696','H:\8TB4\Data\C198-live\Raw',{'CFP','YFP','RFP'},'COPY',false);
IX_Micro_rename('H:\Nalin\201901003-C198-IF1\2019-10-05\10699','H:\8TB4\Data\C198-live-IF1\Raw',{'DAPI','YFP','FarRed'});
IX_Micro_rename('H:\Nalin\201901003-C198-IF2\2019-10-07\10701','H:\8TB4\Data\C198-live-IF2\Raw',{'DAPI','YFP','FarRed'});

% or if need to combine multiple nikon 
% nikon_rename_time('H:\Nikon\C197_live\Raw2','H:\Nikon\C197_live\Raw1')

%% Calculate live image bias
clear p
%%%Directories
p.biaspath='H:\8TB4\Data\C198-live\Bias';     % Output path for bias
p.imagepath='H:\8TB4\Data\C198-live\Raw';     % Image directory. If nikon, single folder for all sites
p.shadingpath='F:\MATLAB\BGimages\IX_Liu\cmosoffset_bin1.mat'; %cmos offset file for IX micro
p.names={
    'CFP';
    'YFP';
    'RFP';
    };
p.row_mat = [2:6];
p.col_mat = [2:11];
p.site_mat = [1:4];
p.frame_mat=[1];
p.biasAll = 0;    % Save averaged bias across sites
                  % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
p.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'; 
p.microscope = 'IXM';

%%% Settings
p.maskforeground = 1;  % 1 to mask nuclei
p.nucradius=12;%12;
p.blur_radius = 5;
p.adaptive = [0 0 0];  % use add channels together to segment foreground
p.dilate = {.75,1,1};  % multiples of nucradius, default .5, to expand mask
p.foreground_calc = 0; % 1 for foreground calc (average foreground signal)
p.method = 'block';    % or pixel
p.blocknum= 15;
p.prctilethresh = 50;
p.compress = .25;
p.prctile_thresh=[0 100]; % pixel remove outliers
p.sigma = 25;             % pixel smooth window


biasCalc(p,0)

%% Calculate Bias IF (empty well)
clear p
%%%Directories
p.biaspath='H:\8TB4\Data\C198-IF1\Bias';                       % Output path for bias
p.imagepath='H:\8TB4\Data\C198-IF1\Raw';                       % Image directory. If nikon, single folder for all sites
p.shadingpath='F:\MATLAB\BGimages\IX_Liu\cmosoffset_bin1.mat'; % cmos offset file for IX micro
p.names={
    'CFP';
    'YFP';
    'RFP';
    };
p.row_mat = [1:4];
p.col_mat = [1];
p.site_mat = [1:16];
p.frame_mat=[1];
p.biasAll = 0;  % Save averaged bias across sites
                % Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
p.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'; 
p.microscope = 'IXM';

%%% Settings
p.maskforeground = 0;  % 1 to mask nuclei
p.nucradius=12;%12;
p.blur_radius = 5;
p.adaptive = [0 0 0];  % use add channels together to segment foreground
p.dilate = {.75,1,1};  % multiples of nucradius, default .5, to expand mask
p.foreground_calc = 0; % 1 for foreground calc (average foreground signal)
p.method = 'block';    % or pixel
p.blocknum= 15;
p.prctilethresh = 50;
p.compress = .25;
p.prctile_thresh=[0 100]; % pixel remove outliers
p.sigma = 25;             % pixel smooth window

biasCalc(s,0);

%% Prep processing
clear s
%%% Paths
s.experiment_name='C198-live';                % Experiment name
s.image_drive = 'H:\8TB4\Data\C198-live\';    % Location of images base folder
s.savepath=fullfile('F:\Data\C-Cdt1\',s.experiment_name,'Data'); % Output location of mat files
s.imagepath=fullfile(s.image_drive,'Raw');                       % All sites in same folder for nikon
s.biaspath=fullfile(s.image_drive,'Bias');                       % Location of bias files
s.maskpath=fullfile(s.image_drive,'Mask');                       % Future location of mask saving
s.bgcmospath='C:\Users\Meyerlab\Documents\GitHub\cell-cycle-tracking\BGimages\IX_Liu\cmosoffset_bin1.mat';
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
s.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff'; 
s.microscope = 'IXM';

%%% General parameters
s.startFrame = 1; 
%s.endFrame = 100;
s.endFrame = getMaxFrame(s.imagepath,'Time(\d+)_');

s.magnification=10;               % 10=10x or 20=20x
s.binsize=1;                      % 1=bin 1 or 2=bin 2
s.postbin = 0;             % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
s.signals={'CFP_','YFP_','RFP_'};
s.maskwrite = 1;                  % 1 = save an image with the nucmask; 0 = dont save an image with the nucmask
s.maskname = 'nucedge_';
s.register = 0;                   % Register images between timepoints for all
s.register_exception = [];        % If don't register images, manually put in sites to register

%%% Quantification parameters
s.bgcmoscorrection = 1;
s.bias = [1 1 1];
s.signal_foreground = [0 0 0];
s.bgsubmethod = {'global nuclear','global nuclear','global nuclear'}; 
                                  % Options:'global nuclear','global cyto','tophat','semi-local nuclear', 'none'
s.compression = 4;                % For bg subtraction
s.bgprctile = [25 25 25];         % Set bg as perctile for each channel
s.sigblur = [3 3 3]; 
s.localbgmeasure = [0 0 0];       % Measure local background in box around cell
s.ringcalc = [0 0 0];
s.ringthresh = [0 50 50];         % Threshold for foreground ring pixels
s.punctacalc = 0;                 % Analyze puncta in nucleus
s.punctaThresh = [125 150 200];   % Threshold for tophat filter for puncta
s.varThresh = [75 100 125];       % Threshold for variance filter puncta 

%%% Segmentation parameters
s.firstsegmethod = 'concavity';  % Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
s.secondsegmethod = 'concavity'; % Options: 'concavity', 'log', 'single adapt', 'double adapt', 'single', 'double'
s.nucr = 12;                     % 10x bin1 = 12 and 20x bin2 = 12
s.blurradius  =  3;              % 10x: 3
s.soliditythresh = 0.8;          % Minimum solidity of nucleus allowed
s.debrisarea = 100;              % 10x = 100 % 150 if no mitosis
s.boulderarea = 1500;            % 10x = 1500
s.blobthreshold = -0.02;         % For blobdector 'log' segmentation

%%% Tracking parameters
s.maxjump = s.nucr*3;
s.masschangethreshold = 0.30;
s.areachangethreshold = 0.60;
s.daughtervariance = 0.10;

%% Run timelapse on sample site
Timelapse(s,2,2,3,0)

%%% Run timelapse parallelized
rows = [2:7];
cols = [2:11];
sites = [1:4];
Parallelwell({@Timelapse},s,rows,cols,sites)

%% Prep match IF to live processing
%%% Paths
clear s;
s.data_path = 'F:\Data\C-Cdt1\C198-live\Data\';      % Output location of mat files
s.IF_imagesessions = {'H:\8TB4\Data\C198-IF1\','H:\8TB4\Data\C198-IF2\'};    
                                                            % Cell array containing all imaging sessions to be matched
s.live_imagepath = 'H:\8TB4\Data\C198-live\Raw\';    % Live imaging raw image path
s.bgcmospath = 'C:\Users\Meyerlab\Documents\GitHub\cell-cycle-tracking\BGimages\IX_Liu\cmosoffset_bin1.mat';
s.crop_save = 'H:\8TB4\Data\C198-IFcrop\';           % Output cropped images, leave empty if no save
s.mask_save = 'H:\8TB4\Data\C198-IFcrop\';           % Output mask, leave empty if no save

%%% General parameters
s.microscope = 'IXM';
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
s.formatCodeLive = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
s.formatCodeFixed = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';
s.magnification = 20;     % 10 = 10x or 20 = 20x
s.binsize = 1;            % 1 = bin 1 or 2 = bin 2
s.postbin = .5;           % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
s.match10x20x = 1;
s.numsites = 4;         % Number of sites imaged for 20x to 10x (double check order is correct)
s.scaleLive = 1;          % Scale live to have same pixel size as final IF image
s.signals = {{'DAPI_','YFP_','FarRed_'}, {'DAPI_','YFP_','FarRed_'}};  % Each imaging session in a cell inside cell array
s.nucLive = 'CFP_';       % Label for IX micro live nucleus (assume first tiff in stack for nikon)
s.manualLiveFrame = [];   % Manually set last frame for live imaging, otherwise use tracedata if empty

%%% Quantification parameters
s.bias = {[1 1 1], [1 1 1]};
s.biasall = {[0 0 0], [0 0 0]};
s.sigblur = {[0 0 0], [0 0 0]};
s.signal_foreground = {[0 0 0], [0 0 0]};
s.bleedthrough = {[0 0 0], [0 0 0]};   % Calculate bleedthrough for channel
s.bleedthroughslope = {};              % Cell array of cell arrays
s.bleedthroughoff = {};                % Cell array of cell arrays
s.ringcalc = {[1 1 1], [0 1 1]};
s.ringthresh = {[0 0 0 ], [0 0 0]};    % Threshold for foreground ring pixels
s.punctacalc = {[0 0 0], [0 0 0]};     % Analyze puncta in nucleus
s.punctaThresh = {{[],[],[]},{}};      % Threshold for tophat filter for puncta
s.punctatopsize = 2;                   % Top hat filter size
s.localbg = {[0 0 0], [0 0 0]};        % Measure local background in box around cell
s.minringsize = 100;
s.bgcmoscorrection = 1;                % 1 = correct for shading; 0 = dont correct for shading;
s.compression = 4;
s.bgsubmethod = {{'global nuclear','global nuclear', 'global nuclear'}, ...
    {'global nuclear', 'global nuclear', 'global nuclear'}}; % 'global nuclear','global cyto','tophat','semi-local nuclear'
s.bgperctile = {[25 25 25], [25 25 25]};  % Set bg as perctile for each channel
s.frameIF=1;
 
%%% Segmentation parameters
s.segmethod = {'concavity','thresh'};  % 'log' or 'single' or 'double'
s.nucr = 12;                           % 10x bin1 = 12 and 20x bin2 = 12
s.debrisarea = 100;                    % 100
s.boulderarea = 1500;                  % 1500 
s.blobthreshold = -0.03;               % For blobdector 'log' segmentation
s.blurradius = 3;                      % Blur for blobdetector
s.soliditythresh = 0.80;               % Minimum solidity of nucleus allowed
s.badFrameCheck = .25;                 % Min fraction of cells untracked before logging error, 0 for no check

%%% Tracking parameters
s.distthresh = 3*s.nucr;
s.arealowthresh = -.4;
s.areahighthresh = .5;

Timelapse_addIF(s,5,2,1,0)


%% Run timelapse
rows = [2:7];
cols = [2:11];
sites = [1:4];
Parallelwell({@Timelapse_addIF},s,rows,cols,sites)