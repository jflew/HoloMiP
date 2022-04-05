% holoAllPositions_alphaX.m
%
% James Flewellen 2021 - HoloMiP Suite
%
% This alpha version is a final version for hologram processing.
% 27-Sep-2021 JF
%
% This is STEP THREE out of THREE.
%
% 1. holoSetParameters_alpha
% 2. holoSetThreshold_alpha
% 3. holoAllPositions_alphaX
%
% Step 1 sets parameters. Step 2 assesses and sets threshold, subvolume 
% halfwidth and Gaussian blur filter radius. Step 3 processes the holograms
% to output 3D positions.
%
% This particular _alpha version is based on the _xenon version of the
% development suite (hence _alphaX).
%
% After this script has run, you can then run your own routines to analyse
% the 3D positions of the candidates through time.
%
% ----------------------------------------------------------------------
%               ******** 3. PROCESS ALL POSITIONS ********
%
% STEP THREE in suite to determine 3D positions of microscopic particles 
% via holographic reconstruction.
%
% This script loads the parameter file stored in STEP TWO (or the AUXILIARY
% script, if the parameters have been changed).
%
% The user inputs some parameters for the wholseale processing of all
% candidate positions throughout a dataset, then the script iterates
% through each frame in the dataset, localising each microparticle
% candidate using the parabolic masking technique.
%
% This script can also have the position of the parfor loop adjusted for
% processing timing tests.
%
% The output of positions, as well as other metadata and timing
% information, is saved if the user selects this option.

% ----------------------------------------------------------------------
%                   ******** OPERATING NOTES ********
%
% To determine the z position of microparticle candidates, this algorithm
% identifies a likely z candidate by generating a line profile along z of
% the subvolumes of interest, and then deploying a peak-finder on this
% line. *In some cases* the peak returned is a spurious object, not the
% genuine object of interest. To avoid missing any genuine peaks, the
% algorithm stores both the MAXIMUM peak returned to find a z-output AND 
% the one furthest away. This is to ensure something obvious isn't missed.
% 
% There are options for truncating the parabolic mask, the truncation
% factor if so, and whether to renormalise the masked sub-volumes. I have
% found that truncating the parabolic mask by a factor of 1, but not
% renormalising the masked subvolumes works best. These are the default
% options in the user-interface questions.
% 
% This script also gives options for testing whether there are spurious
% peaks within a certain radius (rTest) and also for excluding candidates
% that are deemed to be 'too close' to one another (within radius defined
% by rExcl). 
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
% holoSetThreshold_alpha.m (and / or holo1FrameVisualise_alpha.m) has been
% run to store a 'holoParamsFull.mat' file in the 'parameters' 
% sub-directory
%
% ----------------------------------------------------------------------
%                  ******** INPUT & OUTPUT ********
% 
% User selects the base directory containing the data. The script loads the
% holoParamsFull.mat file output from STEP TWO.
%
% The user then needs to input:
%    - a maximal number of particles per field of view
%    - the tolerance value for parabolic masking convergence
%    - the max nunmber of iterations to perform
%    - whether to use the maximum pixel method for xy localisation
%    - which method to use for z location initial guess
%    - whether to save output or not
%
% The script outputs:
%    - a 3D array of positions of candidates in each frame (allPosns)
%       1. row dimension = index of frames
%       2. col dimension = index of particles (aka candidates)
%       3. 3rd dimension stores the salient information as:
%           xPosnGuess  | yPosnGuess  | zPosnGuess  | ...
%           xPosnParab1 | yPosnParab1 | zPosnParab1 | ...
%           xPosnParab2 | yPosnParab2 | zPosnParab2 | ...
%           xPosnFinal  | yPosnFinal  | zPosnFinal  | ...
%           nIterations before convergence | ...
%           clippingFlag
%    - times taken for various stages in the processing
%    - user inputs of processing parameters as metadata
%
%   The relevant outputs for subsequent analysis are allPosns(:,:,10:12),
%   which are the final x,y,z coordinates for all candidates through time.
%
% ----------------------------------------------------------------------
%                  ******* PROCESSING OPTIONS *******
% 
% It is highly recommended to run this script on a CPU with multiple cores
% operating in parallel. The number of workers able to be run
% simultaneously will depend on the available RAM of the machine and the
% size of the reconstructed hologram volume.
%
% You can change the number of workers (search for nWorkers) for your own
% specifications. You might need to experiment with your own set-up to
% determine the optimal value for each dataset. (line 402)
%
% If you are not running on a system with multiple processors, change the
% parfor loop to a for loop. (line 409)
%
% =======================================================================


%% Preliminaries
clc;
close all;
clear all;

%% User inputs required
% Select directory:
disp('***** Select base directory containing data *****')
baseDirectory = uigetdir;
disp(baseDirectory)
disp(' ')

% ------------------------------------------------------------------------
% Load parameter file:
disp(' ')
disp('**********           Thank you              **********')
disp('**********           Loading ...            **********')
load(fullfile(baseDirectory,'parameters','holoParamsFull.mat'));
disp(' ')

% ------------------------------------------------------------------------
%% Processing parameters:
% Request processing parameters file:
paramFileTF = input('Would you like to import a processing parameters file? [Default = N]: ','s');

if paramFileTF == 'Y'
    paramFileTF = 1;
else
    paramFileTF = 0;
end

if paramFileTF == 1 %Load processing parameter file:
    [paramFilename, paramPathname] = uigetfile('*.mat', 'Select .mat processing parameters file');
    load(fullfile(paramPathname,paramFilename))
    disp('*** Processing parameters loaded ***')
    disp(' ')
else %Or else input parameters here:
    % Confirm max number of particles for storage array:
    nParticlesMax = input('What is the maximum number of particles in a given frame? [Default = 200]: ');
    if isempty(nParticlesMax)
        nParticlesMax = 200;
    end
    
    % ==== Modification 25-Mar-2020 JF ====
    % Radius for testing multiple peaks of a single object:
    rTest = input('What is the testing radius for multiple peaks (µm)? [Default = 1]: ');
    if isempty(rTest)
        rTest = 1;
    end 
    
    % Threshold above which multiple candidates are considered viable
    % multiple objects vs phantom objects:
    maxPxTh = input('What is the max px value threshold (%)? [Default = 85%]: ');
    if isempty(maxPxTh)
        maxPxTh = 85;
    end
    maxPxTh = maxPxTh/100; %Convert to decimal
    
    % =====================================
    
    % ==== Modification 17-Mar-2020 JF ====
    % Radius for excluding identified objects that are too close:
    rExcl = input('What is the exclusion radius for object identification (µm)? [Default = 2]: ');
    if isempty(rExcl)
        rExcl = 2;
    end 
    % =====================================
    
    % Confirm value of tolerance for parabolic masking convergence:
    tol = input('What is the error tolerance for convergence? [Default = 0.005]: ');
    if isempty(tol)
        tol = 0.005;
    end
    
    % Confirm max number of iterations to attempt parabolic masking:
    nItsMaxForConvergence = input('What is the maximum number of parabolic masking iterations allowed? [Default = 300]: ');
    if isempty(nItsMaxForConvergence)
        nItsMaxForConvergence = 300;
    end
    
    % Confirm user would like to use (the faster) maxPixel method for xy localisation:
    maxPxTF = input('Use the Max-Pixel method for xy localisation? [Default Y]: ', 's');
    if (maxPxTF == 'N')
        maxPxTF = 0;
    else
        maxPxTF = 1;
    end
    
    % Ask which z-guess method the user would like:
    disp(' ')
    disp('Options for the initial z guess:')
    disp('    1: z-line through kerneled data at maxPx location.')
    disp('    2: Max-Pixel through kerneled sub-volume.')
    disp('    3: Sum-Projection through kerneled sub-volume.')
    
    zGuessOptionHappy = 0;
    while zGuessOptionHappy == 0
        zGuessOption = input('Select a method (1, 2 or 3) [Default = 1]: ');
        if isempty(zGuessOption)
            zGuessOption = 1;
        end
        if (zGuessOption==1) || (zGuessOption==2) || (zGuessOption==3)
            zGuessOptionHappy = 1;
        end
    end
    
    disp(' ')
    % Ask if user would like to use parabolic fitting and recentering of subvolume:
    recentreTF = input('Would you like to recentre the subvolumes? [Default Y]: ', 's');
    if (recentreTF == 'N')
        recentreTF = 0;
    else
        recentreTF = 1;
    end
    
    % Ask if user would like to truncate the parabolic mask function:
    truncateTF = input('Would you like to truncate the parabolic masking function? [Default Y]: ', 's');
    if (truncateTF == 'N')
        truncateTF = 0;
        truncFactor = NaN;
    else
        truncateTF = 1;
    end
    
    % If truncation selected, ask what the factor should be:
    if truncateTF == 1
        truncFactorHappy = 0;
        while truncFactorHappy == 0
            truncFactor = input('Please input truncation factor: [Default = 1]: ');
            if isempty(truncFactor)
                truncFactor = 1;
            end
            if (truncFactor>0)
                truncFactorHappy = 1;
            end
        end
    end
    
    disp(' ')
    % Ask if user would like to use normalised subvolume & parabolic mask:
    normSubVolTF = input('Would you like to normalise the subvolumes and parabolic masks? [Default N]: ','s');
    if (normSubVolTF == 'Y')
        normSubVolTF = 1;
    else
        normSubVolTF = 0;
    end

end

disp(' ')
% Confirm the user would like to save the parameters:
saveTF = input('Save position output as holoAllRawPositions.mat? [Default Y]: ', 's');
if (saveTF == 'N')
    saveTF = 0;
else
    saveTF = 1;
end

if saveTF == 0
    disp('Output will *not* be saved automatically.')
end


%% Load hologram information:
% Directory with hologram data:
% Set directory for hologram data (normalised or not)
if pRecon.preNormTF == 1
    holoDir = 'hologramsNorm';
else
    holoDir = 'holograms';
end
dataDirectory = fullfile(baseDirectory,holoDir);

% List of hologram files:
list_of_holograms = dir(dataDirectory);

% Remove entries that begin with '.'
% (directories: '.' and '..', and '.DSstore')
while list_of_holograms(1).name(1)== '.'
    list_of_holograms(1)=[];
end

% Output just the filenames from the dir operation:
list_of_holograms_cell = struct2cell(list_of_holograms);
names_of_holograms = list_of_holograms_cell(1,:);

% ------------------------------------------------------------------------
% Confirm user wants to use all frames:
disp(' ')
disp(['There are ',num2str(pRecording.nHolos), ' holograms in this dataset.'])
allHolosTF = input('Would you like to use all of them? [Default Y]: ', 's');
if (allHolosTF == 'N')
    allHolosTF = 0;
else
    allHolosTF = 1;
end

if allHolosTF == 0
    startEndFlag = 0; %Flag to check start and end frames are appropriate
    while startEndFlag == 0
        startFrame = input('Input start frame: ');
        endFrame = input('Input end frame: ');
        
        if (startFrame <= endFrame) && (startFrame > 0) && (endFrame <= pRecording.nHolos)
            startEndFlag = 1;
        else
            disp('** Warning: Re-evaluate start & end frames. **')
        end
    end
else
    startFrame = 1;
    endFrame = pRecording.nHolos;
end

nHolos_to_work_on = endFrame-startFrame+1;

disp(['Script will run from frame ',num2str(startFrame),' to frame ',num2str(endFrame),'.'])

% ------------------------------------------------------------------------
%% Begin running:
%Begin total timing:
totalTime = tic;
%Begin prelims timing:
prelimTime = tic;

% Test for whether data are pre-normalised:
if pRecon.preNormTF == 0
    % Load background image for normalisation:
    bkgrnd = imread(fullfile(baseDirectory,'background','background.tif'));
    bkgrnd = double(bkgrnd); % Convert to class double
end

% Initialise storage arrays:
% This is a 3-dimensional array:
% 1. row dimension = index of frames
% 2. col dimension = index of particles (aka candidates)
% 3. 3rd dimension stores the salient information as:
% xPosnGuess | yPosnGuess | zPosnGuess | ...
% xPosnParabFit | yPosnParabFit | zPosnParabFit | ...
% xPosnParabFit2 | yPosnParabFit2 | zPosnParabFit2 | ...
% xPosnFinal | yPosnFinal | zPosnFinal | ...
% nIterations before convergence | ...
% clippingFlag

allPosns = zeros(nHolos_to_work_on,nParticlesMax,14);
frames = (startFrame:endFrame)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ********************       ALL FRAMES         *********************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Processing times
nCandInit = zeros(length(frames),1); %no. candidates initially identified
nCandProc = zeros(length(frames),1); %no. candidates actually processed
tInit = zeros(length(frames),1);
tRecon = zeros(length(frames),1);
tXYGuess = zeros(length(frames),1);
tLocAllCand = zeros(length(frames),1);
tLocPerCand = zeros(length(frames),1);
frameList = zeros(length(frames),1);

nWorkers = NaN; %Leave this in even if just using for loop as variable is saved.

%% Initialise parallel pool
% Delete any existing parallel pools extant
delete(gcp('nocreate'));
% Number of available workers:
nWorkers = 14; % <<=============================== adjust this to work on max number of clusters (max 14) ***************
% Initialise:
parpool(nWorkers);

tPrelim = toc(prelimTime);

%% Iteration through hologram frames in the video sequence:
parfor ff = 1:nHolos_to_work_on %ff = frame
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ********************        INITIALISE        *********************
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    initTime = tic; % Time initialisation
    %disp(['ff = ',num2str(ff)]);
    frameID = frames(ff);
    %disp(['frameID = ',num2str(frameID)]);
    %% Display frame progress:
    %disp(' ')
    disp(['Working on frame ',num2str(frameID),'. ',num2str(ff),' / ',num2str(nHolos_to_work_on),'...']);
    
    %% Create a temporary storage array for each frame:
    % This is necessary to copy into the master allPosns array at the end
    % of the parfor loop: 
    allPosns1frame = zeros(nParticlesMax,14);
    
    tInit(ff) = toc(initTime); % Time initialisation
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ********************      RECONSTRUCTION      *********************
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    reconTime = tic; % Time reconstruction (incl normalisation etc)
    
    %% Read hologram frame:
    holoName = names_of_holograms{frameID};
    holoFile = fullfile(dataDirectory,holoName);
    holo = imread(holoFile);
    % Convert to class double
    holo = double(holo);
    
    %% Normalise if necessary:
    if pRecon.preNormTF == 0
        % Subtract background
        sub = holo-bkgrnd;
        % Divide
        holoNorm = sub./(2*sqrt(bkgrnd));
    elseif pRecon.preNormTF == 1
        holoNorm = holo;
        %holoNorm = holo-median(median(holo)); % Pre-normalised case
    end
    
    %% Initialise storage array for 3D reconstruction:
    reconArray = zeros(pRecon.nR,pRecon.nC,length(pRecon.dList));
    
    %% Reconstruct using Angular Spectrum method
    
    % Fourier transform of the hologram
    holoFT=fftshift(fft2(holoNorm));
    % Multiply hologram FT with Gaussian low-pass filter
    holoFTFilt = holoFT.*pRecon.H;
    % Re-assign variable name incorporating special case for no low-pass filtering:
    if pRecon.D0 == 0
        holoFT_to_use = holoFT;
    else
        holoFT_to_use = holoFTFilt;
    end
    
    % Iterate through each slice
    for ii = 1:length(pRecon.dList)
        d=pRecon.dList(ii);
        % The Angular Spectrum transfer function:
        angSpec = exp((2*pi*1i*d/pRecon.lambda)*pRecon.sqrtPartOfTF);
        % Convolve transfer fn with hologram FT and take IFT
        reconstruction = ifft2(holoFT_to_use.*angSpec);
        % Save to array
        reconArray(:,:,ii) = reconstruction;
    end
    
    tRecon(ff) = toc(reconTime); % Time reconstruction (incl normalisation etc)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % **************      INITIAL XY CANDIDATE GUESS      ***************
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xyGuessTime = tic; % Time xy guess
    
    %% Find intensity of reconstruction:
    reconArrayInt = abs(reconArray).^2;
    
    %% Use either Max-Px or Sum-Projection method:
    if maxPxTF == 1
        % ======= Max-Pixel method: =======
        %% Assign to a 2D plane the value of the max pixel along z
        reconMaxPx = max(reconArrayInt,[],3);
             
        reconMaxPx_orig = reconMaxPx; % Original for mesh display
        % Set pixels below the threshold to zero
        maxVal = max(max(reconMaxPx));
        reconMaxPx(reconMaxPx<(pLocalise.thPc/100)*maxVal) = 0;
        
        %% Identify local maxima for x & y estimate candidate positions:
        localMax = imregionalmax(reconMaxPx); %Binary array
        linCandidates = find(localMax); %Linear indices
        % Convert to subscript indices:
        [yCandidates,xCandidates] = ind2sub(size(reconMaxPx),linCandidates);
        nCandidates = length(linCandidates);
        
        % ==== Modification 25-Mar-2020 JF ====
        % Store the maxPx value for each candidate:
        maxPxVals = zeros(nCandidates,1);
        for mm = 1:nCandidates
            maxPxVals(mm) = reconMaxPx(yCandidates(mm),xCandidates(mm));
        end
        % =====================================
        
    else
        % ======= Sum-Projection method: =======
        %% Generate binary mask of thresholded voxels:
        maskArray = zeros(size(reconArrayInt));
        % Mask: everything ABOVE threshold shows as white (=1).
        % Threshold is thPc x the max identified voxel in this frame's
        % reconstructed volume.
        maskArray(reconArrayInt>(pLocalise.thPc/100)*(max(max(max(reconArrayInt)))))=1;
        % Convert mask to logical array
        maskArray=maskArray~=0;
        
        %% Masked array
        reconArrayIntMasked = zeros(size(reconArrayInt),'like',reconArrayInt);
        reconArrayIntMasked(maskArray) = reconArrayInt(maskArray);
        
        %% Sum along optical (z) axis:
        xySumProj = sum(reconArrayIntMasked,3);
        %xyMeanProj = mean(reconArrayIntMasked,3);
        
        %% Identify local maxima for x & y estimate candidate positions:
        localMax = imregionalmax(xySumProj); %Binary array
        linCandidates = find(localMax); %Linear indices
        % Convert to subscript indices:
        [yCandidates,xCandidates] = ind2sub(size(xySumProj),linCandidates);
        nCandidates = length(linCandidates);
        
    end
    
    %% Test for candidate particles close to edge of hologram and remove:
    % 'Close to edge' means 'within defined halfwidth value + 1'
    % --> change to 2*hw just in case of recentering.
    for cc_a = nCandidates:-1:1 %cc_a = 'all' candidates
        if ((pLocalise.hw_xy+1) >= xCandidates(cc_a)) || ...
                ((pLocalise.hw_xy+1) >= pRecon.nC-xCandidates(cc_a)) || ...
                ((pLocalise.hw_xy+1) >= yCandidates(cc_a)) || ...
                ((pLocalise.hw_xy+1) >= pRecon.nR-yCandidates(cc_a))
            % Remove candidate data from storage array:
            xCandidates(cc_a) = []; yCandidates(cc_a) = []; linCandidates(cc_a) = [];
        end
    end
    nCandidates = length(linCandidates);
    
    nCandInit(ff) = nCandidates; % Store number of candidates
    
    % Store original (x,y) candidates:
    xCandidatesOrig = xCandidates; yCandidatesOrig = yCandidates;
    nCandidatesOrig = nCandidates;
    
    % ==== Modification 25-Mar-2020 JF ====
    
    % Testing radius in pixels:
    rTestPx = rTest/pRecon.effPx;
    
    %storage array for candidates that pass the test:
    candidateTest = zeros(nCandidates,1);
    
    for cc_A = 1:nCandidates
        rSep = zeros(size(xCandidatesOrig)); % Temp array listing all separation radii
        for cc_B = 1:nCandidates
            rSep(cc_B) = sqrt(  ( xCandidatesOrig(cc_A) - xCandidatesOrig(cc_B) )^2 + ...
                ( yCandidatesOrig(cc_A) - yCandidatesOrig(cc_B))^2  );   
        end
        rSep = abs(rSep); % absolute values of all separation radii values
        
        % First test to see if there are any candidates within the testing
        % radius; if none then the candidate of interest gets a pass:
        
        % Copy rSep and remove zero value:
        rSepCopy = rSep;
        rSepCopy(rSepCopy==0) = [];
        
        if min(rSepCopy) > rTestPx
            candidateTest(cc_A) = 1;
        else
            
            % Now, identify all candiates that are within the testing radius:
            candWithin_rTest = find(rSep < rTestPx);
            % Remove the value corresponding to the candidate being investigated:
            candWithin_rTest(candWithin_rTest==cc_A) = [];
            % Create array to test maxPx values:
            maxPxValsTest = zeros(length(candWithin_rTest),1);
            % Populate this with maxPx values:
            for vv = 1:length(maxPxValsTest)
                maxPxValsTest(vv) = maxPxVals(candWithin_rTest(vv));
            end
            
            % ========== JF modification 3-June-2020 ==========
            % If maxPxValsTest is an empty array, set it to zero. This
            % avoids a null result in the event of only a single particle
            % detected
            
            if isempty(maxPxValsTest)
                maxPxValsTest = 0;
            end
            % =================================================
            
            
            % Test to see if maxPx value of the investigated candidate is greater than the others:
            if maxPxVals(cc_A) > max(maxPxValsTest)
                candidateTest(cc_A) = 1;
                % If not greater value, thn test to see if still within the
                % threshold value of the largest value in the set:
            elseif maxPxVals(cc_A) > maxPxTh*max(maxPxValsTest)
                candidateTest(cc_A) = 1;
            end
        end
    end
    
    % Remove all candidates that fail the test:
    xCandidatesPass = xCandidates(candidateTest == 1);
    yCandidatesPass = yCandidates(candidateTest == 1);
    nCandidatesPass = length(xCandidatesPass);

    % ==== Modification 17-Mar-2020 JF ====
    
    % Make new storage array for all candidates that are NOT within exclusion radius:    
    xCandidates = zeros(size(xCandidatesPass));
    yCandidates = xCandidates;
    
    % Exclusion radius in pixels:
    rExclPx = rExcl/pRecon.effPx;
    
    % Iterate through inter-candidate distances identifying candidates that
    % are NOT within excl radius and store in new array:
    
    for cc_A = 1:nCandidatesPass
        rSep = zeros(size(xCandidatesPass)); % Temp array listing all separation radii
        for cc_B = 1:nCandidatesPass
            rSep(cc_B) = sqrt(  ( xCandidatesPass(cc_A) - xCandidatesPass(cc_B) )^2 + ...
                ( yCandidatesPass(cc_A) - yCandidatesPass(cc_B))^2  );    
        end
        
        rSep = abs(rSep); % absolute values
        rSep(rSep==0) = []; %remove zero value
        %[cc_A xCandidatesOrig(cc_A) yCandidatesOrig(cc_A) min(rSep)] %for testing
        
        % ========== JF modification 3-June-2020 ==========
        % If rSep is an empty array (occurs when only a single bead in the
        % FoV) then set x & y candidates to the 'passed' candidates from
        % before.
        
        if isempty(rSep)
            xCandidates(cc_A) = xCandidatesPass(cc_A);
            yCandidates(cc_A) = yCandidatesPass(cc_A);
        elseif min(rSep) > rExclPx
            xCandidates(cc_A) = xCandidatesPass(cc_A);
            yCandidates(cc_A) = yCandidatesPass(cc_A);
        end
        % =================================================
        

    end
    
    xCandidates(xCandidates==0) = []; %remove zero values
    yCandidates(yCandidates==0) = []; %remove zero values  
    nCandidates = length(xCandidates);
    
    % Double-up the x and y candidates s.t. the max z position can also be saved: <<<<<<<<<<======================================= MODIFIED
    xCandidates = [xCandidates' xCandidates']';
    yCandidates = [yCandidates' yCandidates']';
    
    %% Store xy candidate guess in temporary output array:
    allPosns1frame(1:2*nCandidates,1) = xCandidates;
    allPosns1frame(1:2*nCandidates,2) = yCandidates;
    
    tXYGuess(ff) = toc(xyGuessTime); % Time xy guess
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % *********      APPLY 3D KERNEL FILTER TO SUBVOLUMES      **********
    % *********                Z CANDIDATE GUESS                *********
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    locTime = tic; % Time zGuess + parabolic masking for all candidates in this frame
    
    % Counter to keep track of number of candidates NOT processed
    %... i.e. deliver a NaN in initial z-guess
    nNaNs = 0; 
    
    for cc = 1:2*nCandidates
        
        %disp(['Candidate ',num2str(cc)])
        % Extract a subvolume 'chimney pipe' along z
        chimneyInt = reconArrayInt(yCandidates(cc)-pLocalise.hw_xy:yCandidates(cc)+pLocalise.hw_xy,...
            xCandidates(cc)-pLocalise.hw_xy:xCandidates(cc)+pLocalise.hw_xy,:);
        
        %% Apply kernel filter
        % i.e. the volume has been 'Kerneled':
        chimneyKrnld = zeros(size(chimneyInt));
        
        % Dimensions over which to apply kernel:
        nRvol = size(chimneyInt,1);
        nCvol = size(chimneyInt,2);
        nZvol = size(chimneyInt,3);
        
        % Iterate to pass the kernel across the sub-volume:
        for mm=1:nRvol-2 %results in a border of zeroes around the array
            for nn=1:nCvol-2
                for oo=1:nZvol-2
                    S1=sum(sum(sum(pLocalise.S.*chimneyInt(mm:mm+2,nn:nn+2,oo:oo+2))));
                    chimneyKrnld(mm+1,nn+1,oo+1)=S1;
                end
            end
        end
        
        %% Find z-Guess:
        % Options 1, 2 or 3 as described at top of script:
        if zGuessOption == 1
            % Taking the z-line through the identified xy centre:
            zLine = squeeze(chimneyKrnld(pLocalise.hw_xy+1,pLocalise.hw_xy+1,:));
            zLine(zLine<0)=0;
        elseif zGuessOption == 2
            % Use the Max-Px method through the kerneled subvolume:
            chimneyKrnldMaxPx2 = squeeze(max(chimneyKrnld,[],2));
            chimneyKrnldMaxPx = max(chimneyKrnldMaxPx2,[],1);
            zLine = chimneyKrnldMaxPx;
        elseif zGuessOption == 3
            % Use the Sum-Projection method through the kerneled subvolume:
            chimneyKrnldSumOnce = sum(chimneyKrnld,1);
            chimneyKrnldSumOnce = squeeze(chimneyKrnldSumOnce);
            chimneyKrnldSumTwice = sum(chimneyKrnldSumOnce,1);
            zLine = chimneyKrnldSumTwice;
        end
        
        % Simple peak finder, sorted from tallest to shortest:
        [pks, locs, w, p] = findpeaks(zLine,'SortStr','descend',...
            'MinPeakHeight',0.10*max(zLine),... % Ensure peaks are all positive and at least 10% of the max value in the krnl
            'NPeaks',3,... % Ensure only top 3 peaks are returned
            'MinPeakProminence',pLocalise.minPeakProm,... % Minimum peak prominence
            'MinPeakDistance',20); % Peaks must be at least 20z-px apart
        
        % Remove all peaks that are within a subvolume halfwidth of the floor or ceiling
        pks(locs<(pLocalise.hw_z+1)) = []; %floor
        locs(locs<(pLocalise.hw_z+1)) = []; %floor
        
        pks(locs>(nZvol-pLocalise.hw_z-1)) = []; %ceiling
        locs(locs>(nZvol-pLocalise.hw_z-1)) = []; %celing
        
        % Conditions to save BOTH the highest peak AND the furthest peak <<<<<<<<<<================ MODIFIED
        
        if cc <= nCandidates
            % Check for NaNs:
            if isempty(pks)
                zGuess = NaN;
                %disp(' ')
                disp(['No peaks found for candidate ',num2str(cc), ' in frame ',num2str(frameID),'.'])
                nNaNs = nNaNs+1; % Counting number of NaNs
            else
                % Take the furthest peak:
                [maxLocs,maxLocsPos] = max(locs);
                zGuess = locs(maxLocsPos);
            end
            
        elseif cc > nCandidates
            % Check for NaNs:
            if isempty(pks)
                zGuess = NaN;
                %disp(' ')
                disp(['No peaks found for candidate ',num2str(cc), ' in frame ',num2str(frameID),'.'])
                nNaNs = nNaNs+1; % Counting number of NaNs
            else
                % Take highest peak (this is a backup):
                zGuess = locs(1);
            end
        end
        
        % ============================================================== <<<<<<<<<<================ MODIFIED
        
        %% Store returned peak value:
        allPosns1frame(cc,3) = zGuess;
             
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % **************          PARABOLIC MASKING           ***************
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% Parabolic masking to determine accurate 3D position...
        % ... as long as the zGuess is not NaN
        
        if ~isnan(zGuess)
            
            % Flag to indicate clipping of recentred subvolume:
            svClipFlag = 0;
            
            % Cuboid subvolume of the kerneled data specified by the initial
            % x,y,z position guesses and the prescribed halfwidth.
            subvolKrnld = chimneyKrnld(:,:,zGuess-pLocalise.hw_z:zGuess+pLocalise.hw_z);
            
            % Clip off the zero-boundaries from the sub-volume:
            subvolKrnld(2*pLocalise.hw_xy+1,:,:) = []; %rows
            subvolKrnld(1,:,:) = [];
            
            subvolKrnld(:,2*pLocalise.hw_xy+1,:) = []; %cols
            subvolKrnld(:,1,:) = [];
            
            % Translate so all values are positive and begin at zero
            subvolKrnld = subvolKrnld-max(max(max(subvolKrnld)));
            subvolKrnld = subvolKrnld + abs(min(min(min(subvolKrnld))));
            % Normalise if desired:
            if normSubVolTF == 1
                subvolKrnld = subvolKrnld/sum(sum(sum(subvolKrnld)));
            end
            
            % Create arrays of x, y & z positions relating to the subarray:
            [Xs,Ys,Zs] = meshgrid(xCandidates(cc)-(pLocalise.hw_xy-1):xCandidates(cc)+(pLocalise.hw_xy-1),...
                yCandidates(cc)-(pLocalise.hw_xy-1):yCandidates(cc)+(pLocalise.hw_xy-1),...
                zGuess-pLocalise.hw_z:zGuess+pLocalise.hw_z);

            % x,y,z pixel location values:
            xx = (xCandidates(cc)-(pLocalise.hw_xy-1):xCandidates(cc)+(pLocalise.hw_xy-1))';
            yy = (yCandidates(cc)-(pLocalise.hw_xy-1):yCandidates(cc)+(pLocalise.hw_xy-1))';
            zz = (zGuess-pLocalise.hw_z:zGuess+pLocalise.hw_z)';
            
            % Parabolic fit to determine new centre of subvolume:
            if recentreTF == 1
                % Extract lines along x,y,z through the centre voxel:
                xLineForFitting = squeeze(subvolKrnld(pLocalise.hw_xy,:,pLocalise.hw_z+1))';
                yLineForFitting = squeeze(subvolKrnld(:,pLocalise.hw_xy,pLocalise.hw_z+1));
                zLineForFitting = squeeze(subvolKrnld(pLocalise.hw_xy,pLocalise.hw_xy,:));
                 
                % Fit a parabola along each axis:
                [fx,gx] = fit(xx,xLineForFitting,'poly2'); %fx = x-fit; gx = goodness of x-fit
                [fy,gy] = fit(yy,yLineForFitting,'poly2'); %fy = y-fit; gy = goodness of y-fit
                [fz,gz] = fit(zz,zLineForFitting,'poly2'); %fz = z-fit; gz = goodness of z-fit
                
                % New centre:
                xCentreParabFit1 = -fx.p2/(2*fx.p1);
                yCentreParabFit1 = -fy.p2/(2*fy.p1);
                zCentreParabFit1 = -fz.p2/(2*fz.p1);
                
                %% Store xyz parabolic fitting outputs in temporary array:
                allPosns1frame(cc,4) = xCentreParabFit1;
                allPosns1frame(cc,5) = yCentreParabFit1;
                allPosns1frame(cc,6) = zCentreParabFit1;
                
                %% Recentre subvolume on the new centres from parabolic fit:
                % Round new centres to find integer values:
                xCentreParabFitRound = round(xCentreParabFit1);
                yCentreParabFitRound = round(yCentreParabFit1);
                zCentreParabFitRound = round(zCentreParabFit1);
                
                %% Test for clipping:
                % Test if the recentred z-value is now too close to the
                % floor or ceiling of the sample:
                if zCentreParabFitRound <= (pLocalise.hw_z+1)
                    zCentreParabFitRound = pLocalise.hw_z+2;
                    svClipFlag = 1;
                end
                if zCentreParabFitRound > (length(pRecon.dList)-(pLocalise.hw_z+1))
                    zCentreParabFitRound = (length(pRecon.dList)-(pLocalise.hw_z+1));
                    svClipFlag = 1;
                end
                % Test if the recentred x- or y-value is now too close to 
                % the edge of the hologram:
                if xCentreParabFitRound <= (pLocalise.hw_xy)
                    xCentreParabFitRound = pLocalise.hw_xy+1;
                    svClipFlag = 1;
                end
                if xCentreParabFitRound >= (pRecon.nR-(pLocalise.hw_xy-1))
                    xCentreParabFitRound = (pRecon.nR-(pLocalise.hw_xy));
                    svClipFlag = 1;
                end
                if yCentreParabFitRound <= (pLocalise.hw_xy)
                    yCentreParabFitRound = pLocalise.hw_xy+1;
                    svClipFlag = 1;
                end
                if yCentreParabFitRound >= (pRecon.nC-(pLocalise.hw_xy-1))
                    yCentreParabFitRound = (pRecon.nC-(pLocalise.hw_xy));
                    svClipFlag = 1;
                end
                
                % Store clipping flag:
                allPosns1frame(cc,14) = svClipFlag;
                
                
                %% Newly-centred subvolume of reconstructed intensities:
                newSVInt = reconArrayInt(yCentreParabFitRound-pLocalise.hw_xy:yCentreParabFitRound+pLocalise.hw_xy,...
                    xCentreParabFitRound-pLocalise.hw_xy:xCentreParabFitRound+pLocalise.hw_xy,...
                    zCentreParabFitRound-(pLocalise.hw_z+1):zCentreParabFitRound+(pLocalise.hw_z+1));
                % Apply kernel filter to this new subvolume:
                newSVKrnld = zeros(size(newSVInt));
                
                % Dimensions over which to apply kernel:
                nRvol = size(newSVInt,1);
                nCvol = size(newSVInt,2);
                nZvol = size(newSVInt,3);
                
                % Iterate to pass the kernel across the sub-volume:
                for mm=1:nRvol-2 %results in a border of zeroes around the array
                    for nn=1:nCvol-2
                        for oo=1:nZvol-2
                            S1=sum(sum(sum(pLocalise.S.*newSVInt(mm:mm+2,nn:nn+2,oo:oo+2))));
                            newSVKrnld(mm+1,nn+1,oo+1)=S1;
                        end
                    end
                end
                
                % Clip off the zero-boundaries from the sub-volume:
                newSVKrnld(2*pLocalise.hw_xy+1,:,:) = []; %rows
                newSVKrnld(1,:,:) = [];
                
                newSVKrnld(:,2*pLocalise.hw_xy+1,:) = []; %cols
                newSVKrnld(:,1,:) = [];
                
                newSVKrnld(:,:,2*pLocalise.hw_z+3,:) = []; %z
                newSVKrnld(:,:,1) = [];
                % (z is +3 as we added 2 to the newSVInt definition
                % above to keep the dimensions the same.)
                
                % Translate so all values are positive and begin at zero
                newSVKrnld = newSVKrnld-max(max(max(newSVKrnld)));
                newSVKrnld = newSVKrnld + abs(min(min(min(newSVKrnld))));
                % Normalise if desired:
                if normSubVolTF == 1
                    newSVKrnld = newSVKrnld/sum(sum(sum(newSVKrnld)));
                end
                
                % Create arrays of x, y & z positions relating to the new subvolume:
                [Xs,Ys,Zs] = meshgrid(xCentreParabFitRound-(pLocalise.hw_xy-1):xCentreParabFitRound+(pLocalise.hw_xy-1),...
                    yCentreParabFitRound-(pLocalise.hw_xy-1):yCentreParabFitRound+(pLocalise.hw_xy-1),...
                    zCentreParabFitRound-pLocalise.hw_z:zCentreParabFitRound+pLocalise.hw_z);
                
                %% Repeat fitting on new subvolume for updated centre positions:
                % x,y,z pixel location values, for plotting (and fitting)
                xx = (xCentreParabFitRound-(pLocalise.hw_xy-1):xCentreParabFitRound+(pLocalise.hw_xy-1))';
                yy = (yCentreParabFitRound-(pLocalise.hw_xy-1):yCentreParabFitRound+(pLocalise.hw_xy-1))';
                zz = (zCentreParabFitRound-pLocalise.hw_z:zCentreParabFitRound+pLocalise.hw_z)';
                
                %% Repeat parabola fit on the re-centred subvolume: 
                % Extract lines along x,y,z through the centre voxel:
                xLineForFitting = squeeze(newSVKrnld(pLocalise.hw_xy,:,pLocalise.hw_z+1))';
                yLineForFitting = squeeze(newSVKrnld(:,pLocalise.hw_xy,pLocalise.hw_z+1));
                zLineForFitting = squeeze(newSVKrnld(pLocalise.hw_xy,pLocalise.hw_xy,:));
                
                % Fit a parabola along each axis:
                [fx,gx] = fit(xx,xLineForFitting,'poly2'); %fx = x-fit; gx = goodness of x-fit
                [fy,gy] = fit(yy,yLineForFitting,'poly2'); %fy = y-fit; gy = goodness of y-fit
                [fz,gz] = fit(zz,zLineForFitting,'poly2'); %fz = z-fit; gz = goodness of z-fit
                
                % New centre:
                xCentreParabFit2 = -fx.p2/(2*fx.p1);
                yCentreParabFit2 = -fy.p2/(2*fy.p1);
                zCentreParabFit2 = -fz.p2/(2*fz.p1);
                
                %% Store xyz parabolic fitting outputs in temporary array:
                allPosns1frame(cc,7) = xCentreParabFit2;
                allPosns1frame(cc,8) = yCentreParabFit2;
                allPosns1frame(cc,9) = zCentreParabFit2;
                            
            end  

            %% Initialise values for parabolic masking
            nIts = 1;                       %Iteration loop counter
            error1 = 1;                     %Convergence criterion
            %disp(' ')
            %disp(['Candidate ',num2str(cc)])
            
            % Kerneled subvolume and initial position seed:
            if recentreTF == 1
                SVKrnld2Use = newSVKrnld;
                xCentre = allPosns1frame(cc,7);
                yCentre = allPosns1frame(cc,8);
                zCentre = allPosns1frame(cc,9);
            else
                SVKrnld2Use = subvolKrnld;
                xCentre = allPosns1frame(cc,1);
                yCentre = allPosns1frame(cc,2);
                zCentre = allPosns1frame(cc,3);
            end
            
            % Paraboloid constants:
            A = -1;
            B = -1;
            C = -1;
            
            while (error1>tol || nIts<11)
                %disp(['Iteration ',num2str(nIts),': xCentre = ',num2str(xCentre),'. yCentre = ',num2str(yCentre),'. zCentre = ',num2str(zCentre)]);
                
                % ----- PARABOLA MASK -----
                % Create parabola mask
                parabMask = A*(Xs-xCentre).^2 + ...
                    B*(Ys-yCentre).^2 + ...
                    C*(Zs-zCentre).^2;

                if truncateTF == 1
                    % If mask truncation selected, translate so all values 
                    % are positive, and set to zero all voxels below 
                    % truncFactor * median value 
                    parabMask = parabMask - max(max(max(parabMask)));
                    parabMask = parabMask + truncFactor*abs(median(median(median(parabMask))));
                    parabMask(parabMask<0) = 0;
                else
                    % Otherwise translate so all values are positive and 
                    % begin at zero
                    parabMask = parabMask - max(max(max(parabMask)));
                    parabMask = parabMask + abs(min(min(min(parabMask))));
                end

                % Multiply the kerneled subvolume by the parabola mask:
                SV = SVKrnld2Use.*parabMask;

                % Calculate revised estimates of x, y & z centres:
                xCentreNew = sum(sum(sum(SV.*Xs)))/sum(sum(sum(SV)));
                yCentreNew = sum(sum(sum(SV.*Ys)))/sum(sum(sum(SV)));
                zCentreNew = sum(sum(sum(SV.*Zs)))/sum(sum(sum(SV)));
                % Difference new-old
                deltaX = xCentreNew-xCentre;
                deltaY = yCentreNew-yCentre;
                deltaZ = zCentreNew-zCentre;
                % Determine error1 value to decide upon convergence
                error1 = sqrt(deltaX^2 + deltaY^2 + deltaZ^2);
                
                % Update values for next iteration:
                xCentre = xCentreNew;
                yCentre = yCentreNew;
                zCentre = zCentreNew;
                nIts = nIts+1; %iteration count
                % ----------------------------------
                % Prevent infinite loop:
                if nIts>nItsMaxForConvergence
                    noConvergenceFlag = 1;
                    break
                else
                    noConvergenceFlag = 0;
                end
            end
            
            %% Returning output from parabola masking:
            % Checking convergence is reached
            if noConvergenceFlag == 1
                allPosns1frame(cc,10) = xCentre;
                allPosns1frame(cc,11) = yCentre;
                allPosns1frame(cc,12) = zCentre;
                allPosns1frame(cc,13) = NaN; % <<==================================== Flag to show no convergence
            else
                % Final output
                allPosns1frame(cc,10) = xCentre;
                allPosns1frame(cc,11) = yCentre;
                allPosns1frame(cc,12) = zCentre;
                allPosns1frame(cc,13) = nIts;
            end
        else
            % Final output filled with NaNs indicating not a good starting
            % z candidate
            allPosns1frame(cc,10) = NaN;
            allPosns1frame(cc,11) = NaN;
            allPosns1frame(cc,12) = NaN;
            allPosns1frame(cc,13) = NaN;
            allPosns1frame(cc,14) = NaN;
        end
    end
    % Store number of candidates actually processed:
    nCandProc(ff) = nCandidates-nNaNs; 
    % Timing:
    tLocAllCand(ff) = toc(locTime);
    tLocPerCand(ff) = toc(locTime)/(nCandidates-nNaNs);
    
    %% Now store these vectors in the allPosns storage array:
    allPosns(ff,:,:)=allPosns1frame;
    % And store the order in which frames were worked upon:
    frameList(ff)=frameID;
end

%% Confirm processing completion:
disp(' ')
disp('**********      Processing complete.        **********')
disp(' ')

% Final timing
tTotal = toc(totalTime);
tProces = table(frames,frameList, nCandInit, nCandProc, tInit, tRecon, tXYGuess, tLocAllCand, tLocPerCand);
disp(['Total time = ', num2str(tTotal),'s.'])

% Print to screen:
disp(['Preliminaries time = ', num2str(tPrelim),'s.'])
disp(['Mean initialisation time = ', num2str(mean(tInit)),'s.'])
disp(['Mean reconstruction time (',num2str(pRecon.nR),'px x ',...
    num2str(pRecon.nC),'px x ',num2str(length(pRecon.dList)),'z-px) = ',...
    num2str(mean(tRecon)),'s.'])
disp(['Mean initial xy guess finding time = ',num2str(mean(tXYGuess)),'s.'])
disp(['Over all ',num2str(nHolos_to_work_on),' frames, ',...
    num2str(sum(nCandInit)),' candidates were initally identfied,'])
disp(['out of which ',num2str(sum(nCandProc)),' were processed,']);
disp(['taking a mean time of ',num2str(mean(tLocPerCand)),'s,'])
disp(['and a median time of ', num2str(median(tLocPerCand)),'s.'])
disp(' ')
disp('=============')
disp(' ')

% Print user parameters:
disp('Processing metadata parameters:')
fprintf('nParticlesMax = %i \n',nParticlesMax)
fprintf('rTest = %i \n',rTest)
fprintf('maxPxTh = %i \n',maxPxTh)
fprintf('rExcl = %i \n',rExcl)
fprintf('tol = %i \n',tol)
fprintf('nItsMaxForConvergence = %i \n',nItsMaxForConvergence)
fprintf('maxPxTF = %i \n',maxPxTF)
fprintf('zGuessOption = %i \n',zGuessOption)
fprintf('recentreTF = %i \n',recentreTF)
fprintf('normSubVolTF = %i \n',normSubVolTF)
fprintf('truncateTF = %i \n',truncateTF)
fprintf('truncFactor = %f \n',truncFactor)

pProces.nParticlesMax = nParticlesMax;
pProces.rTest = rTest;
pProces.maxPxTh = maxPxTh;
pProces.rExcl = rExcl;
pProces.tol = tol;
pProces.nItsMaxForConvergence = nItsMaxForConvergence;
pProces.maxPxTF = maxPxTF; 
pProces.zGuessOption = zGuessOption;
pProces.recentreTF = recentreTF;
pProces.normSubVolTF = normSubVolTF;
pProces.truncateTF = truncateTF;
pProces.truncFactor = truncFactor;

pProces.time_of_completion = datetime('now');
pProces.dataDir = dataDirectory;


%% Save output

if saveTF == 1
    % Save to the 'output' directory:
    savePath = fullfile(baseDirectory,'output','holoAllRawPositions.mat');
    save(savePath,'pRecording','pRecon','pLocalise','pProces',...
        'tTotal','tProces','tPrelim','nWorkers',...
        'frameList',...
        'allPosns');
    disp(' ')
    disp('**********      Data have been saved.       **********')
    disp('================================================================')
else
    disp('================================================================')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

