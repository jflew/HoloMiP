% holoSetThreshold_alpha.m
%
% James Flewellen 2021 - HoloMiP Suite
%
% This alpha version is a final version for hologram processing.
% 24-Sep-2018 JF
%
% This is STEP TWO out of THREE.
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
% based on. Also holoAllPositions_Sep2018_C.m for aspects of the
% reconstruction process.
%
% ----------------------------------------------------------------------
%                  ******** 2. SET THRESHOLD ETC ********
%
% STEP TWO in suite to determine 3D positions of microscopic particles via
% holographic reconstruction.
%
% This script loads the parameter file stored in STEP ONE, then loads a
% hologram of choice to to assess the effects of setting values for:
%   - intensity threshold for MaxPx masking
%   - Gaussian blur radius (D0)
%   - halfwidth value of subvolume box for parabolic masking
%
% The user can also investigate the kernel lines along the optical axis to
% check the correct peak is being detected.
%
% The script saves the parameters into a .mat file to be called by STEP
% THREE in the suite.
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
% holoSetParameters_alpha.m has been run to store a 'holoParams.mat' file
% in the 'parameters' sub-directory.
%
% ----------------------------------------------------------------------
%                  ******** INPUT & OUTPUT ********
% 
% User selects the base directory containing the data. The script loads the
% holoParams.mat file output from STEP ONE.
%
% User then follows the text prompts to input reconstruction parameters.
%
% The script outputs a 'holoParamsFull.mat' file to be called later.
%
% ----------------------------------------------------------------------
%                  ********     NOTES     ********
% 
% There is no hard and fast rule for setting intensity threshold, Gaussian
% blur radius, or the xy or z halfwidths.
% 
% A 3D mesh map of the un-thresholded normalised hologram is produced,
% along with the maximum pixel intensity value. This can help in setting an
% appropriate threshold value.
%
% To see the effect of the Gaussian filter, set D0 to zero (i.e. no
% filter). For a severe Gaussian smoothing, set a low value of D0; 
% for a subtle Gaussian smoothing, set a high value. 
% Ideally, choose the subtlest Gaussian filter that allows effective 
% object localisation (i.e. smoothes out noisy pixels, but still presents 
% precise peaks). What a 'good' value is also depends on the pixel 
% dimension of the hologram.
%
% The output provided when selecting individual object candidates to
% investigate will help with this. I generally use option 1 (taking a line
% along z through the kerneled data at the maxPx location), and in fact I
% have commented out the other options in this version. Feel free to
% uncomment these lines and to experiment with different options.
%
% Halfwidth values can be inferred from this output. You are looking,
% ideally, for a parabolic profile.
%
% In general, several passes through the script are required, adjusting 
% for threshold, D0, and halfwidth - as all variables are dependent on each 
% other. 
%
% You may also want to compare the optimal values found in one frame, with
% a frame at a later point in the experiment, which may have different
% diffraction patterns (e.g. if the objects are further away from, or closer
% to, the focus at a later point in time.
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

% Script loads extant parameters file
load(fullfile(baseDirectory,'parameters','holoParams.mat'))

% Set directory for hologram data (normalised or not)
if pRecon.preNormTF == 1
    holoDir = 'hologramsNorm';
else
    holoDir = 'holograms';
end

% Ask user if they would like graphical output:
graphicalTF = input('Would you like graphical output? [Default Y]: ','s');
if graphicalTF == 'N'
    graphicalTF = 0;
else
    graphicalTF = 1;
end

disp(' ')

%% Get input for reconstruction and localisation parameters:
% Initialise value to test for satisfaction with parameters
happyTF = 0;

while happyTF == 0
    disp('The current reconstruction parameters are:')
    disp(['    Floor = ',num2str(pRecon.dFloor),'µm'])
    disp(['    Ceiling = ',num2str(pRecon.dCeiling),'µm'])
    disp(['    Spacing = ',num2str(pRecon.dSpacing),'µm'])
    
    changeReconParamTF = input('Would you like to change these? [Default = N]: ','s');
    if changeReconParamTF == 'Y'
        changeReconParamTF = 1;
    else
        changeReconParamTF = 0;
    end
    if changeReconParamTF == 1
        pRecon.dFloor = input('New value (in µm) for reconstruction floor: ');
        pRecon.dCeiling = input('New value (in µm) for reconstruction ceiling: ');
        pRecon.dSpacing = input('New value (in µm) for reconstruction spacing: ');
        
        pRecon.dList = (pRecon.dFloor:pRecon.dSpacing:pRecon.dCeiling)'; %in µm 
    end
    
    disp(' ')
    disp('***** Input localisation parameters: *****')
    
    % Threshold value for mask identification of initial position guess
    thPc = input('Value of threshold as a % of maximum pixel magnitude [Default = 50%]: ');
    if isempty(thPc)
        thPc = 50;
    end
    pLocalise.thPc = thPc;
    
    % Radius value of Gaussian low-pass filter in Fourier plane during recon
    D0 = input('Gaussian low-pass filter radius (in inverse px) [Default = 100]: ');
    if isempty(D0)
        D0 = 100;
    end
    pRecon.D0 = D0;
    
%     % Ask whether user would like to use (the faster) maxPixel method for xy localisation:
%     maxPxTF = input('Use the Max-Pixel method for xy localisation? [Default Y]: ', 's');
%     if (maxPxTF == 'N')
%         maxPxTF = 0;
%     else
%         maxPxTF = 1;
%     end
    
    % Halfwidth of cube around particle for parabolic localisation
    hw_xy = input('xy-halfwidth of cube around particle for parabolic masking [Default = 7px]: ');
    if isempty(hw_xy)
        hw_xy = 7;
    end
    pLocalise.hw_xy = hw_xy;
    
    hw_z = input('z-halfwidth of cube around particle for parabolic masking [Default = 14px]: ');
    if isempty(hw_z)
        hw_z = 14;
    end
    pLocalise.hw_z = hw_z;
    
    % Value of minimum peak prominence
    minPeakProm = input('Value of the minimum peak prominence for z-guess? [Default = 100]: ');
    if isempty(minPeakProm)
        minPeakProm = 100;
    end
    pLocalise.minPeakProm = minPeakProm;
    
    %% Generating auxiliary transfer function:
    % Generate the pixel arrays for the transfer function
    % Generate meshgrid
    [r,c]=meshgrid(1:pRecon.nR,1:pRecon.nC);
    % Recentre on zero
    rCentred = r-pRecon.nR/2;%-1; %Some implementations have a -1 here...
    cCentred = c-pRecon.nC/2;%-1; %...Perhaps depends on whether indexing begins at 0 or 1
    % Now take the transpose to ensure dimensions are same in subsequent
    % muliplication
    rCentred=rCentred'; cCentred=cCentred';
    % Now multiply these arrays by wavelength / holo size and then square:
    varRsq=((pRecon.lambda/pRecon.sizeR)*rCentred).^2; %var for 'variable'
    varCsq=((pRecon.lambda/pRecon.sizeC)*cCentred).^2;
    % Combine everything that is under the sqrt in the transfer function:
    underSqRt=1-varRsq-varCsq;
    % Take square root:
    sqrtPartOfTF=sqrt(underSqRt);
    
    pRecon.sqrtPartOfTF = sqrtPartOfTF;
    
    %% Generating low-pass Gaussian function to convolve with hologram FT
    % Generate Gaussian function (H)
    % Use same meshgrid as above
    D = sqrt(cCentred.^2 + rCentred.^2);
    H = exp((-D.^2)./(2*(pRecon.D0^2)));
    
    pRecon.H = H;
    
    %% Display graphical output if selected
    if graphicalTF == 1
        disp(' ')
        disp(['*** There are ',num2str(pRecording.nHolos), ' holograms in this dataset. ***'])
        
        holoFrameHappy = 0;
        while holoFrameHappy == 0
            holoFrame = input('Which frame would you like to look at? [Default = 1]: ');
            if isempty(holoFrame)
                holoFrame = 1;
            end
            
            if (holoFrame <= pRecording.nHolos) && (holoFrame > 0)
                holoFrameHappy = 1;
            else
                disp('** Warning: Frame out of range. Reselect. **')
            end
        end
        disp(' ')
        disp(['***** Working on frame ',num2str(holoFrame),'... *****'])
        
        
        %% Loading hologram file 'holograms' or 'hologramsNorm' sub-directory:
        if pRecon.preNormTF == 0
            dataDirectory = fullfile(baseDirectory,'holograms');
        else
            dataDirectory = fullfile(baseDirectory,'hologramsNorm');
        end
        
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
        
        % Read hologram:
        holoName = names_of_holograms{holoFrame};
        holoFile = fullfile(dataDirectory,holoName);
        holo = imread(holoFile);
        % Convert to class double
        holo = double(holo);
        
        % Load background image for normalisation, if needed:
        if pRecon.preNormTF == 0
            bkgrnd = imread(fullfile(baseDirectory,'background','background.tif'));
            bkgrnd = double(bkgrnd); % Convert to class double
            % Subtract
            sub = holo-bkgrnd;
            % Divide
            holoNorm = sub./(2*sqrt(bkgrnd));
        else
            holoNorm = holo; % Pre-normalised case
        end
        
        %% Reconstruction via Angular Spectrum method:
        % Initialise storage array for 3D reconstruction:
        reconArray = zeros(pRecon.nR,pRecon.nC,length(pRecon.dList));
        
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
        
        reconArrayInt = abs(reconArray).^2;
        
%         %% Use either Max-Px or Sum-Projection method:
%         if maxPxTF == 1
%             % ======= Max-Pixel method: =======
%             %% Assign to a 2D plane the value of the max pixel along z
%             reconMaxPx = max(reconArrayInt,[],3);
%             % Set pixels below the threshold to zero
%             maxVal = max(max(reconMaxPx));
%             reconMaxPx(reconMaxPx<(pLocalise.thPc/100)*maxVal) = 0;
%             
%             %% Identify local maxima for x & y estimate candidate positions:
%             localMax = imregionalmax(reconMaxPx); %Binary array
%             linCandidates = find(localMax); %Linear indices
%             % Convert to subscript indices:
%             [yCandidates,xCandidates] = ind2sub(size(reconMaxPx),linCandidates);
%             nCandidates = length(linCandidates);
%             
%         else
%             % ======= Sum-Projection method: =======
%             %% Generate binary mask of thresholded voxels:
%             maskArray = zeros(size(reconArrayInt));
%             % Mask: everything ABOVE threshold shows as white (=1).
%             % Threshold is thPc x the max identified voxel in this frame's
%             % reconstructed volume.
%             maskArray(reconArrayInt>(pLocalise.thPc/100)*(max(max(max(reconArrayInt)))))=1;
%             % Convert mask to logical array
%             maskArray=maskArray~=0;
%             
%             %% Masked array
%             reconArrayIntMasked = zeros(size(reconArrayInt),'like',reconArrayInt);
%             reconArrayIntMasked(maskArray) = reconArrayInt(maskArray);
%             
%             %% Sum along optical (z) axis:
%             xySumProj = sum(reconArrayIntMasked,3);
%             %xyMeanProj = mean(reconArrayIntMasked,3);
%             
%             %% Identify local maxima for x & y estimate candidate positions:
%             localMax = imregionalmax(xySumProj); %Binary array
%             linCandidates = find(localMax); %Linear indices
%             % Convert to subscript indices:
%             [yCandidates,xCandidates] = ind2sub(size(xySumProj),linCandidates);
%             nCandidates = length(linCandidates);
%         end
         
        
        %% Using Max-Px method to find xy initial guesses:
        % Assign to a 2D plane the value of the max pixel along z
        reconMaxPx = max(reconArrayInt,[],3);
              
        reconMaxPx_orig = reconMaxPx; % Original for mesh display
        % Set pixels below the threshold to zero
        maxVal = max(max(reconMaxPx));
        reconMaxPx(reconMaxPx<(pLocalise.thPc/100)*maxVal) = 0;
        
        % Identify local maxima for x & y estimate candidate positions:
        localMax = imregionalmax(reconMaxPx); %Binary array
        linCandidates = find(localMax); %Linear indices
        % Convert to subscript indices:
        [yCandidates,xCandidates] = ind2sub(size(reconMaxPx),linCandidates);
        nCandidates = length(linCandidates);
        
        %% Display:
        % Display normalised hologram & results of xy-localisation
        f = figure;
        set(f, 'Name','HoloNorm',...
            'Position', [100 100 1200 750])
        subplot(1,2,1);
        imagesc(holoNorm); colormap gray; axis image
        title({['Normalised hologram, frame ',num2str(holoFrame),'.'];...
            baseDirectory}, 'interpreter', 'none')
        
        subplot(1,2,2);
        imagesc(real(log(holoFT_to_use))); colormap gray; axis image
        title({['Fourier transform used, frame ',num2str(holoFrame),'.'];...
            ['D0 = ',num2str(pRecon.D0),'.']})
        
        f = figure;
        set(f, 'Name','Localisation',...
            'Position', [100 100 1200 750]); 
        subplot(1,2,1);
        imagesc(holoNorm); colormap gray; axis image; hold on
        plot(xCandidates,yCandidates,'o',...
            'MarkerSize',8,...
            'LineWidth',2,...
            'Color',[0.35 0.65 0.45])
        for cc = 1:nCandidates
            text(xCandidates(cc)+8,yCandidates(cc)+8,num2str(cc,'%02d'),...
                'FontSize',10,...
                'FontName','Calibri',...
                'Color',[0.35 0.65 0.45]);
        end
        title(['Maximum Pixel projection, frame ',num2str(holoFrame),'.'])
        
        subplot(1,2,2);
        imagesc(reconMaxPx); colormap gray; axis image; hold on
        plot(xCandidates,yCandidates,'o',...
            'MarkerSize',8,...
            'LineWidth',2,...
            'Color',[0.35 0.65 0.45])
        for cc = 1:nCandidates
            text(xCandidates(cc)+8,yCandidates(cc)+8,num2str(cc,'%02d'),...
                'FontSize',10,...
                'FontName','Calibri',...
                'Color',[0.35 0.65 0.45]);
        end
        title(['Maximum Pixel projection, frame ',num2str(holoFrame),'.'])
        
        figure('Name', 'Orig MaxPx Mesh');
        mesh(reconMaxPx_orig);
        xlabel('x'); ylabel('y'); zlabel('intensity (A.U.)');
        title({['Original maxPx output, frame ',num2str(holoFrame)];...
            baseDirectory;...
            ['Max pixel value = ',num2str(maxVal),'.']}, 'interpreter', 'none');
        
        figure('Name', 'Thresholded MaxPx Mesh');
        mesh(reconMaxPx);
        xlabel('x'); ylabel('y'); zlabel('intensity (A.U.)');
        title({['Thresholded maxPx output, frame ',num2str(holoFrame)];...
            baseDirectory;...
            ['Max pixel value = ',num2str(maxVal),'.']}, 'interpreter', 'none');
        
        %% Ask if user wants to investigate zProfiles of any candidates:
        zInvestigateTF = input('Would you like to look at z-profiles of any candidates? [Default = Y]: ','s');
        if zInvestigateTF == 'N'
            zInvestigateTF = 0;
        else
            zInvestigateTF = 1;
        end

        if zInvestigateTF == 1
            % Initialise a test to ensure candidates are appropriate
            candidatesHappy = 0;
            while candidatesHappy == 0
                disp('Please input candidate numbers in form: [a b c ...]')
                candToLookAt = input('Input here: ');
                
                if (isempty(candToLookAt)) || (max(candToLookAt) > nCandidates)
                    disp('Incorrect format')
                else
                    candidatesHappy = 1;
                end
            end
            
            % Ask which z-guess method the user would like:
            disp('Which method would you like to use for the initial z-guess localisation?')
            disp('Options are:')
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

            %% Investigating candidates' z profile etc
            for ll = 1:length(candToLookAt)
                % Selects desired candidate:
                cc = candToLookAt(ll);
                % Extract a subvolume 'chimney pipe' along z
                chimneyInt = reconArrayInt(yCandidates(cc)-pLocalise.hw_xy:yCandidates(cc)+pLocalise.hw_xy,...
                    xCandidates(cc)-pLocalise.hw_xy:xCandidates(cc)+pLocalise.hw_xy,:);
                
                %% max pixel trick on intensity, for display purposes
                chimneyIntMaxPx2 = squeeze(max(chimneyInt,[],2));
                chimneyIntMaxPx = max(chimneyIntMaxPx2,[],1);
                
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
                
                % Just take the furthest peak:
                if isempty(pks)
                    zGuess = NaN;
                    disp(['No peaks found for candidate ',num2str(cc), ' in frame ',num2str(holoFrame),'.'])
                else
                    [maxLocs,maxLocsPos] = max(locs);
                    zGuess = locs(maxLocsPos);
                end
                
                if isnan(zGuess)
                    f = figure;
                    set(f, 'Name',['Subvolume output. Frame ',num2str(holoFrame),'. Candidate ',num2str(cc)],...
                        'Position', [100 100 1200 750])
                    % Original (normalised) hologram:
                    subplot(3,3,1)
                    imagesc(holoNorm(yCandidates(cc)-pLocalise.hw_xy:yCandidates(cc)+pLocalise.hw_xy,...
                        xCandidates(cc)-pLocalise.hw_xy:xCandidates(cc)+pLocalise.hw_xy)); colormap gray; axis image
                    title(['Holo, candidate ',num2str(cc),'.'])
                    
                    % Intensity trace along z
                    subplot(3,3,2)
                    plot(chimneyIntMaxPx,...
                        'Color', [0.35 0.35 0.75],...
                        'LineWidth',2); grid on
                    xlabel('z index (z-px)'); ylabel('Intensity (A.U.)')
                    title({'Intensity (Max-Px)';['D0 = ',num2str(pRecon.D0),'px']})
                    % Kerneled trace along z
                    subplot(3,3,5)
                    plot(zLine,...
                        'Color', [0.75 0.35 0.75],...
                        'LineWidth',2); grid on;
                    xlabel('z index (z-px)'); ylabel('Intensity (A.U.)')
                    title({['Kerneled zLine (option=',num2str(zGuessOption),')'];...
                        ['zGuess = ',num2str(zGuess)]})
                    
                elseif ~isnan(zGuess)
                        f = figure;
                        set(f, 'Name',['Subvolume output. Frame ',num2str(holoFrame),'. Candidate ',num2str(cc)],...
                            'Position', [100 100 1200 750])
                    % Original (normalised) hologram:
                    subplot(3,3,1)
                    imagesc(holoNorm(yCandidates(cc)-pLocalise.hw_xy:yCandidates(cc)+pLocalise.hw_xy,...
                        xCandidates(cc)-pLocalise.hw_xy:xCandidates(cc)+pLocalise.hw_xy)); colormap gray; axis image
                    title(['Holo, candidate ',num2str(cc),'.'])
                    % Reconstructed intensity at z-guess:
                    subplot(3,3,4)
                    imagesc(chimneyInt(:,:,zGuess)); colormap gray; axis image
                    title(['Recon Int, c=',num2str(cc),'; z=',num2str(zGuess),'.'])
                    % Kerneled values at z-guess:
                    subplot(3,3,7)
                    imagesc(chimneyKrnld(:,:,zGuess)); colormap gray; axis image
                    title(['Recon Krnld, c= ',num2str(cc),'; z=',num2str(zGuess),'.'])
                    
                    % Intensity trace along z
                    subplot(3,3,2)
                    plot(chimneyIntMaxPx,...
                        'Color', [0.35 0.35 0.75],...
                        'LineWidth',2); grid on
                    xlabel('z index (z-px)'); ylabel('Intensity (A.U.)')
                    title('Intensity (Max-Px)')
                    % Kerneled trace along z
                    subplot(3,3,5)
                    plot(zLine,...
                        'Color', [0.75 0.35 0.75],...
                        'LineWidth',2); grid on; hold on
                    plot(zGuess,zLine(zGuess),'o',...
                        'Color', [0.25 0.75 0.45],...
                        'MarkerSize',12,...
                        'LineWidth',2);
                    xlabel('z index (z-px)'); ylabel('Intensity (A.U.)')
                    title({['Kerneled zLine (option=',num2str(zGuessOption),')'];...
                        ['zGuess = ',num2str(zGuess)]})
                    
                    % x trace along sub volume kernel
                    subplot(3,3,3)
                    plot(squeeze(chimneyKrnld(pLocalise.hw_xy,2:2*pLocalise.hw_xy-1,...
                        zGuess)),...
                        'Color', [0.75 0.35 0.75],...
                        'LineWidth',2);
                    grid on
                    xlabel('x'); ylabel('A.U.')
                    title(['x Krnld profile. hw = ',num2str(pLocalise.hw_xy),'.'])
                    % y trace along sub volume kernel
                    subplot(3,3,6)
                    plot(squeeze(chimneyKrnld(2:2*pLocalise.hw_xy-1,pLocalise.hw_xy,...
                        zGuess)),...
                        'Color', [0.75 0.35 0.75],...
                        'LineWidth',2);
                    grid on
                    xlabel('y'); ylabel('A.U.')
                    title(['y Krnld profile. hw = ',num2str(pLocalise.hw_xy),'.'])
                    % z trace along sub volume kernel
                    subplot(3,3,9)
                    plot(squeeze(chimneyKrnld(pLocalise.hw_xy,pLocalise.hw_xy,...
                        zGuess-pLocalise.hw_z+1:zGuess+pLocalise.hw_z+1)),...
                        'Color', [0.75 0.35 0.75],...
                        'LineWidth',2);
                    grid on; hold on
                    
                    xlabel('z'); ylabel('A.U.')
                    title(['z Krnld profile. hw = ',num2str(pLocalise.hw_z),'.'])
                    
                end   
            end
        end
    end
    % Ask if happy with input values, if not the programme repeats
    happyTestTF = input('Are you happy with these values? [Default = N]: ','s');
    if happyTestTF == 'Y'
        happyTestTF = 1;
    else
        happyTestTF = 0;
        
        closeTF = 0;
        closeTF = input('Close all open figures? [Default = Y]: ','s');
        if closeTF == 'N'
            closeTF = 0;
        else
            closeTF = 1;
        end
        
        if closeTF == 1
            close all
        end
    end
    disp(' ')

    if happyTestTF == 1
        happyTF = 1;
    end
    
end

%% Save parameters to 'full' parameter file:
saveTF = input('Save these parameters? [Default Y; type N to avoid saving]: ','s');
if saveTF == 'N'
    saveTF = 0;
else
    saveTF = 1;
end

if saveTF == 1    
    % Save to the 'output' directory:
    savePath = fullfile(baseDirectory,'parameters','holoParamsFull.mat');
    save(savePath,'pRecording','pRecon','pLocalise');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



