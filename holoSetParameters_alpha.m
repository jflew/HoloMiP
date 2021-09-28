% holoSetParameters_alpha.m
%
% James Flewellen 2021 - HoloMiP Suite
%
% This alpha version is a final version for hologram processing.
% 24-Sep-2018 JF
%
% This is STEP ONE out of THREE for 3D localisation in the HoloMiP Suite.
%
% 1. holoSetParameters_alpha
% 2. holoSetThreshold_alpha
% 3. holoAllPositions_alpha
%
% Step 1 sets parameters. Step 2 assesses and sets threshold, subvolume 
% halfwidth and Gaussian blur filter radius. Step 3 processes the holograms
% to output 3D positions.
%
% See holoSetParameters_Aug2018.m for the original script this version is
% based on.
%
% ----------------------------------------------------------------------
%                  ******** 1. SET PARAMETERS ********
%
% STEP ONE in suite to determine 3D positions of microscopic particles via
% holographic reconstruction.
%
% This script is used to select the directory containing data of interest 
% and to set the recording and reconstruction parameters.
%
% There is an option to save the parameters into a .mat file, which can be
% called by subsequent functions.
%
% ----------------------------------------------------------------------
%                  ******** DATA ORGANISATION ********
%
% Raw holograms need to be contained in a sub-directory called 'holograms'.
% A background file is stored in a sub-directory called 'background'.
% Pre-normalised holograms (optional) may be stored in a sub-directory 
% called 'hologramsNorm'.
% Parameter files are stored in sub-directory called 'parameters'.
% The 'output' subdirectory is used to store all output from processing.
%
% ----------------------------------------------------------------------
%                  ********    ASSUMPTIONS    ********
%
% Holograms are of an even dimension. Use evenSizer.m script to process
% this a priori.
%
% ----------------------------------------------------------------------
%                  ******** INPUT & OUTPUT ********
% 
% User selects the base directory containing the data. Then follows the
% text prompts to input experimental parameters.
%
% The script outputs a 'holoParams.mat' file to be called later.
%
% ----------------------------------------------------------------------
%                  ********     NOTES     ********
% 
% This technique uses a 'bottom' Sobel-like kernel to apply a gradient
% filter. To change this, you'll need to hard-code it.
%
% =======================================================================

%% Prelimininaries

clc; close all; clear all;

% User selects directory:
disp('***** Select base directory containing data *****')
baseDirectory = uigetdir;
disp('Thank you.')
disp(['Directory = ',baseDirectory])
disp(' ')

% Check whether using pre-normalised holograms:
preNormTF = input('Is the dataset pre-normalised? [Default = N]: ','s');
if preNormTF == 'Y'
    preNormTF = 1;
else
    preNormTF = 0;
end

if preNormTF == 1
    holoDir = 'hologramsNorm';
else
    holoDir = 'holograms';
end

% Check whether 'parameters' directory exists and make it if not:
paramDirExist = exist(fullfile(baseDirectory,'parameters'),'dir');
if paramDirExist ~= 7
    mkdir(baseDirectory,'parameters');
end

% Check whether 'output' directory exists and make it if not:
paramDirExist = exist(fullfile(baseDirectory,'output'),'dir');
if paramDirExist ~= 7
    mkdir(baseDirectory,'output');
end

%% Determine number and size of holograms in dataset
% Directory with hologram data:
dataDirectory = fullfile(baseDirectory,holoDir);

% List of hologram files:
list_of_holograms = dir(dataDirectory);

% Remove entries that begin with '.' 
% (directories: '.' and '..', and '.DSstore')
while list_of_holograms(1).name(1)== '.'
    list_of_holograms(1)=[];
end

% Number of holograms in sub-directory:
nHolos = length(list_of_holograms); 

% Filename of first hologram in sub-directory:
holo1name = list_of_holograms(1).name;
holo1file = fullfile(dataDirectory,holo1name);
holo1 = imread(holo1file);

% Size
[nR,nC] = size(holo1); %returns # rows & cols

disp(['The full hologram size is ',num2str(nR), ' x ', num2str(nC),'.'])
overwriteTF = input('Would you like to overwrite these values? [Default N]: ','s');
if overwriteTF == 'Y'
    overwriteTF = 1;
else
    overwriteTF = 0;
end

if overwriteTF == 1
    nR = input('New value for number of rows: ');
    nC = input('New value for number of columns: ');
end


%% Obtain experimental recording parameters from user:

% ======== LIGHT SETTINGS ========
% Wavelength of illumination in µm:
% (Red LED used on MT set-up has central lambda of 625nm)
lambdaVacuum = input('Vacuum wavelength in µm [Default = 0.625µm]: ');
if isempty(lambdaVacuum)
    lambdaVacuum = 0.625;
end

% Refractive index of medium:
nRefractive = input('Refractive index of medium [Default = 1.333]: ');
if isempty(nRefractive)
    nRefractive = 1.333;
end

% Calculate wavelength in the medium:
lambda = lambdaVacuum/nRefractive;

% ======== VIDEO SETTINGS ========
% Time interval between frames:
frameInt = input('Time interval between frames in ms [Default = 40ms]: ');
if isempty(frameInt)
    frameInt = 40;
end

% Calculate framerate (in frames per second):
framerate = (1/frameInt)*1000;

% ======== MAGNIFICATION SETTINGS ========
% Pixel spacing of camera in µm:
% (Hamamatsu Orca Flash 4 has 6.5µm pixels)
camPx = input('Camera pixel spacing in µm [Default = 6.5µm]: ');
if isempty(camPx)
    camPx = 6.5;
end

% Magnification used in experimental acquisition:
mag = input('Total magnification in image acquisition [Default = 100x]: ');
if isempty(mag)
    mag = 100;
end

% Calculating effective pixel size in µm:
effPx = camPx/mag; %µm; effective pixel size

% Calculating dimensions of hologram in µm:
sizeR = nR*effPx; sizeC = nC*effPx;

% ======== RECONSTRUCTION SETTINGS ========
% Floor of reconstruction volume:
dFloor = input('Floor position of reconstruction in µm [Default = 0.1µm]: ');
if isempty(dFloor)
    dFloor = 0.1;
end

% Ceiling of reconstruction volume:
dCeiling = input('Ceiling position of reconstruction in µm [Default = 30.0µm]: ');
if isempty(dCeiling)
    dCeiling = 30.0;
end

% Spacing between reconstruction planes:
dSpacing = input('Spacing between reconstruction planes in µm [Default = 0.1µm]: ');
if isempty(dSpacing)
    dSpacing = 0.1;
end

%Evaluating list of reconstruction planes:
dList = (dFloor:dSpacing:dCeiling)'; %in µm

% 3D kernel for gradient filter:
% This is a 'Bottom' Sobel:
S = [-1 -2 -1; -2 -4 -2; -1 -2 -1];
S(:,:,2) = [0 0 0; 0 0 0; 0 0 0];
S(:,:,3) = [1 2 1; 2 4 2; 1 2 1];


%% Print variables to screen:
fprintf('\n')
fprintf('========================================== \n')
fprintf('\n')
fprintf('***** EXPERIMENTAL PARAMETERS ***** \n')
fprintf('Data directory: \n')
disp(baseDirectory)
fprintf('Illumination wavelength (vacuum) = %1.4fµm \n',lambdaVacuum) 
fprintf('Illumination wavelength (medium) = %1.4fµm \n',lambda)
fprintf('Refractive index of medium = %1.4fµm \n',nRefractive)
fprintf('Interval between frames = %1.1fms \n',frameInt)
fprintf('Framerate = %1.2ffps \n',framerate)
fprintf('Effective pixel spacing = %1.4fµm \n',effPx)
fprintf('Row dimension of hologram = %ipx \n',nR)
fprintf('Column dimension of hologram = %ipx \n',nC)
fprintf('Row size of hologram = %1.4fµm \n',sizeR)
fprintf('Column size of hologram = %1.4fµm \n',sizeC)
fprintf('Number of planes to reconstruct: %1i \n',length(dList))
fprintf('The kernel filter is a ''Bottom'' Sobel: \n')
disp(S)
fprintf('If you don''t like it, change in the script \n')


%% Option to save the recording parameters
% holoParameters = {};
saveTF = input('Save these parameters? [Default Y; type N to avoid saving]: ','s');
if saveTF == 'N'
    saveTF = 0;
else
    saveTF = 1;
end

if saveTF == 1
    % Save parameters into structs: ('p' = 'Parameters')
    
    pRecording.baseDirectory = baseDirectory;
    pRecording.nHolos = nHolos;
    pRecording.lambdaVacuum = lambdaVacuum;
    pRecording.nRefractive = nRefractive;
    pRecording.camPx = camPx;
    pRecording.mag = mag;
    pRecording.frameInt = frameInt;
    pRecording.framerate = framerate;
    
    pRecon.nR = nR;
    pRecon.nC = nC;
    pRecon.sizeR = sizeR;
    pRecon.sizeC = sizeC;
    pRecon.lambda = lambda;
    pRecon.effPx = effPx;
    pRecon.dFloor = dFloor;
    pRecon.dCeiling = dCeiling;
    pRecon.dSpacing = dSpacing;
    pRecon.dList = dList;
    pRecon.preNormTF = preNormTF;
    
    pLocalise.S = S; 
    

    %% Save to the 'output' directory:
    savePath = fullfile(baseDirectory,'parameters','holoParams.mat');
    save(savePath,'pRecording','pRecon','pLocalise');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

