% holoMedian_tifStack.m
% 
% 5-Feb-2019 James Flewellen
% Modifying holoMedian.m to output a median background image from a tif
% stack, rather than individual files.
%
% INPUT: User selects tif stack file.
%
% OUTPUT: background image saved into a sub-directory called 'background',
% placed in same directory as input. File is a 16-bit uint tif, named 
% 'background.tif'
%
% See holoMedian.m for more details.
% 
% ========================================================================

clc; clear all; close all;

%% User selects tif stack:
disp('***** Select tif stack containing background images: *****')
[filename,baseDirectory] = uigetfile({'*.tif'},'Chose image data sequence (.tif stack):');

%% Get info on tif stack:
tifInfo = imfinfo(fullfile(baseDirectory,filename));
nFrames = length(tifInfo); % number of frames in sequence.
% Dimensions of images:
nR = tifInfo(1).Height;
nC = tifInfo(1).Width;


%% Initialise 
% Initialise 3D image stack to take median through:
holoStack = double(zeros(nR,nC,nFrames));

%% Iterate to read images in selected tif stack:
for ff = 1:nFrames %ff = frame 
    holo = imread(fullfile(baseDirectory,filename),ff);
    holoStack(:,:,ff) = double(holo);
end

%% Take median along the stack dimension:
background = median(holoStack,3);

% ===== DISPLAY OUTPUT: =====
figure('Name','Median background','Position',[50 50 1200 750]);
imagesc(background); colormap gray; axis image;
title(['Background. Median of ',num2str(nFrames),' holos.'])
    
%% Convert to 16-bit for tif saving
background = uint16(background);   

%% Give option to save file:
disp(' ')
saveTF = input('Save background file? [Y/N; Default Y]: ','s');
if (saveTF ~= 'Y')
    saveTF = 0;
else
    saveTF = 1;
end
    
%% Save background as tif file:
if saveTF == 1
    % Make 'background' subdirectory:
    mkdir(baseDirectory,'background1')
    % Save to the 'background1' directory:
    savePath = fullfile(baseDirectory,'background1','background.tif');
    imwrite(background,savePath,'tif')
end
 
% ========================================================================