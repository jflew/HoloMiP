% holoMedian_gold.m
%
% 20-May-2021 James Flewellen
%
% (based on holoMedian.m 2017-20)
% ------------------------
% Generates a MEDIAN image from a set of holograms to be used as a 
% background image for normalising a hologram dataset. 
% Returns an even-sized array if clipping is necessary.
% Saves the median as an 8- or 16-bit uint tif image file named
% 'background.tif'  in the 'background' directory of the dataset parent 
% directory. 
%
% ----------------------------------------------------------------------
%                  ******** DATA ORGANISATION ********
%
% Raw holograms need to be contained in a sub-directory called 'holograms'.
% Holograms for median generation to be used as a background need to be
% contained in a sub-directory called 'findbkgr'.
% The background file is stored in a sub-directory called 'background'.
% The 'output' subdirectory is used to store all output from processing.
% ----------------------------------------------------------------------
%
% User selects the base directory containing the data. This script loads
% all images in the 'findbkgr' sub-directory.
%
% ======================
%
% INPUT: folder of hologram images to be processed. Folder must ONLY have
% the images to be used for median generation.
%
% Hologram images should be .tif format. Can be UINT as 'double' command
% converts to double precision.
%
% =======
%
% OUTPUT: a 8- or 16-bit uint tif image file named 'background.tif'.
% 
% ========================================================================

clc; clear all; close all;

%% User selects directory:
disp('***** Select base directory containing data *****')
baseDirectory = uigetdir;

% Directory with hologram images to find median:
medianDirectory = fullfile(baseDirectory,'findbkgr');

%% Generate list of files
list_of_median_files = dir(medianDirectory);
% Remove entries that begin with '.' 
% (directories: '.' and '..', and '.DSstore')
while list_of_median_files(1).name(1)== '.'
    list_of_median_files(1)=[];
end
    
nIms2Median = length(list_of_median_files); %number of images to find mean of

%% Give option to exit:
disp(['***** ',num2str(nIms2Median),' images found. *****'])
contTF = input('Continue? [Y/N; Default Y]: ','s');
if (contTF ~= 'Y')
    contTF = 0;
else
    contTF = 1;
end

%% Continue with median generation:
if contTF == 1
    
    
    %% Read one image from the corresponding hologram dataset:
    dataDirectory = fullfile(baseDirectory,'holograms');
    list_of_data_files = dir(dataDirectory);
    % Remove entries that begin with '.'
    % (directories: '.' and '..', and '.DSstore')
    while list_of_data_files(1).name(1)== '.'
        list_of_data_files(1)=[];
    end

    holo1name = list_of_data_files(1).name;
    holo1file = fullfile(dataDirectory,holo1name);
    holo1 = imread(holo1file);
    
    % detect class, 8- or 16-bit
    imclassString = class(holo1);
    if imclassString(5) == '8'
       imclass = 8
    elseif imclassString(5) == '1'
        imclass = 16
    end
    
    %% Initialise background array
    background = double(zeros(size(holo1)));
    
    %% Initialis 3D image stack to take median through:
    [Rdim,Cdim]=size(holo1);
    holoStack = double(zeros(Rdim,Cdim,nIms2Median));
        
    %% Iterate to read images in selected directory
    for k = 1:nIms2Median
        holoName = list_of_median_files(k).name;
        holoFile = fullfile(medianDirectory,holoName);
        holo = imread(holoFile);
        holoStack(:,:,k) = double(holo);
        
%         % Add to background
%         background = background+holo;
%         
%         clear holoName holoFile holo
    end
    
    
    %% Take median along the stack dimension:
    background = median(holoStack,3);
    
%     %% Divide background by number of files
%     background = background./nIms2Mean;
%     
    % ===== DISPLAY OUTPUT: =====
    figure('Name','Median background','Position',[50 50 1200 750]);
    subplot(1,2,1)
    imagesc(holo1); colormap gray; axis image;
    title(['Hologram. Frame 1. '])
    subplot(1,2,2)
    imagesc(background); colormap gray; axis image;
    title(['Background. Median of ',num2str(nIms2Median),' holos.'])
    % ==============================
    
    %% Convert to 8- or 16-bit for tif saving
    if imclass == 16
        background = uint16(background);
    elseif imclass == 8
        background = uint8(background);
    end
    % background = background';
    
    %% Invoke evenSizer function to ensure row & col dimensions are even:
    [backgroundEven,nRnew,nCnew] = evenSizer(background);
    [nRold,nCold]=size(background);
    
    evenSizerFlag = 0;
    if (nRnew ~= nRold) || (nCnew ~= nCold)
        evenSizerFlag = 1;
    end
    if evenSizerFlag == 1
        disp('***** Background image has been clipped to ensure size is even. *****')
    end

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
        % Save to the 'background' directory:
        savePath = fullfile(baseDirectory,'background','background.tif');
        imwrite(backgroundEven,savePath,'tif')
    end
    
    else
    return
end

%clearvars -except 
% ========================================================================