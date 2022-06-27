%% Calculate Bias IF (empty well)
clear p
%%%Directories
p.biaspath='H:\Nikon\C197_control_IF1\Bias';                       % Output path for bias
p.imagepath='H:\Nikon\C197_control_IF1\Raw';                       % Image directory. If nikon, single folder for all sites
p.shadingpath=''; % cmos offset file for IX micro
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
p.microscope = 'Nikon';

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

%% Prep fixed imaging
%%% Paths
settings.data_path = 'F:\Data\C-Cdt1\C197-control\Data\';   % Output location of mat files
settings.IF_imagesessions = {'H:\Nikon\C197_control_IF1\'}; % Cell array containing all imaging sessions to be matched
settings.bgcmospath = '';
settings.maskpath = 'H:\Nikon\C197_control_IFcrop\';        % Output mask, leave empty if no save
settings.crop_save = 'H:\Nikon\C197_control_IFcrop\';       % Output cropped images, leave empty if no save
settings.microscope = 'nikon';
% Formatting code for nikon files, 1$ is timepoint, 2$ is well, 3$ is point
settings.formatCode = 'Time%1$05d_Well%2$s_Point%2$s_%3$04d_*.tiff';

%%% File Options
settings.maskname = 'nucedge_';

%%% General parameters
settings.magnification = 20;           % 10 = 10x or 20 = 20x
settings.binsize = 1;                  % 1 = bin 1 or 2 = bin 2
settings.postbin = 0;                  % 0 for no software bin, number for scaling factor (e.g. 0.5 for 2x2 bin);
settings.maxjit = 500;
settings.signals = {{'DAPI_','Cy5_'}};

%%% Quantification parameters
settings.bias = {[1 1], [1 1 1]};
settings.biasall = {[0 0], [0 0 0]};
settings.sigblur = {[0 0], [0 0 0]};
settings.signal_foreground = {[0 0], [0 0 0]};
settings.bleedthrough = {[0 0], [0 0 0]};
settings.bleedthroughslope = {};
settings.bleedthroughoff = {};
settings.ringcalc = {[1 1], [0 0 1]};
settings.ringthresh = {[0 0 ], [0 0 0]};
settings.punctacalc = {[0 0], [0 0 0]};
settings.punctaThresh = {{[],[],[]},{}};
settings.punctatopsize = 2;
settings.localbg = {[0 0], [0 0 0]};
settings.minringsize = 100;
settings.bgcmoscorrection = 1; % 1 = correct for shading; 0 = dont correct for shading;
settings.bgsubmethod = {{'global nuclear', 'global nuclear'}, ...
    {'global nuclear', 'global nuclear', 'global nuclear'}}; 
                               % 'global nuclear','global cyto','tophat','semi-local nuclear'
settings.bgperctile = {[25 25], [25 25 25]}; % Set bg as perctile for each channel

%%% Segmentation parameters
settings.segmethod = {'concavity','thresh'}; % 'log' or 'single' or 'double'
settings.nucr = 24;                          % 10x bin1 = 12 and 20x bin2 = 12
settings.debrisarea = 400;                   % 100
settings.boulderarea = 4500;                 % 1500
settings.blobthreshold = -0.03;              % For blobdector 'log' segmentation
settings.blurradius = 3;                     % Blur for blobdetector
settings.soliditythresh = 0.80;              % Minimum solidity of nucleus allowed
settings.compression = 4;                    % Min fraction of cells untracked before logging error, 0 for no check

%% Test on single well 
Fixed_IF(settings,7,2,1,1)

%% Run timelapse parallelized
rows = [7];
cols = [2:11];
sites = [1:16];
Parallelwell({@Fixed_IF},settings,rows,cols,sites)
