function [appData, handles, ERRORFLG] = Confocal2ColorALExSolnMultDataPhotonRecoloring8(...
    appData, handles)
%
% Photon recoloring simulation. Allow simulation of conformational and
% photophysical dynamics using an empirical burst size distribution and
% photon arrival times.
%
%   
%
%   Edited: 
%       KGryte - (2012-08-08) - Created. Used for ES Burst Correlation
%       Analysis
%                             - non-overlapping windows for E/S and correlation
%                             - window of total photons and not Gex
%                             - mixture sub-populations
%       KGryte - (2012-09-02) - Assume static species and replacing any ES
%                               correlations
%
%   References:
%       Gopich, Szabo. Decoding photon colors.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% StartUp/CleanUp:

StartUp();

C = onCleanup(@() CleanUp());


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization:

ALTERNATIONTIME = 50*10^-6; % 50us Alternation Time
ALTERNATIONPERIOD = 100*10^-6; % 100us Alternation period


ERRORFLG = false;

WINDOW = 10; %appData.BVA.Config.Algorithm.Window; % e.g., 5 photons
CORRWINDOW = 5;%10;
sigAlpha = 0.25;
filterFLG = true;
    

PHOTONCUTOFF = floor(WINDOW/2);

CLUSTERBINS = appData.BVA.Config.Algorithm.Clusters.NumberOfBins; % e.g., 20
CLUSTERWIDTH = (1-0)/CLUSTERBINS;

CLUSTERMINWINDOWS = 1; %appData.BVA.Config.Algorithm.Clusters.Threshold.MinWindows; 

% Define how many photons should be within a burst:
% NUMPHOTONS = 10; % for now...

stoichiometryStateMeans = [...
    0.70;... % S = 0 --> Donor blinking ; S = 1 --> Acceptor blinking
    0.70];


    
% S = 0.7; % Base S state
% stoichiometryStateMeans = [...
%     0.0;...     % donor blinking
%     S/(.67 - (.67-1)*S);... % state with 2/3 AA intensity
%     S;... 0.5 % actual S prob
%     S/(1.33 - (1.33-1)*S);... % state with 4/3 AA intensity
%     1.0]; % 0.7 % acceptor blinking

% efficiencyStateMeans = [...
%     0.5, 0.5, 0.5, 0.5, 0.0;...
%     0.5, 0.5, 0.5, 0.5, 0.0]; % Needs to be a numEStates X numSStates matrix

efficiencyStateMeans = [...
    0.5, 0.5;...
    0.5, 0.5]; % Needs to be a numEStates X numSStates matrix



gamma = 1.0;

% photophysicalInitialProbabilities = [...
%     0.001;...
%     0.05;...
%     0.898;...
%     0.05;...
%     0.001];

photophysicalInitialProbabilities = [...
    0.5;...
    0.5];



conformationalInitialProbabilities = [...
    0.5;...
    0.5];



photophysicalRateMatrix = [...
    -1000,    1000;...
     1000,   -1000];

% photophysicalRateMatrix = [...
%     -0,    0;...
%      0,   -0];

% photophysicalRateMatrix = [...
%      -9000,     3000,       3000,       3000,       0;...   % donor blinking rates [sec^-1]
%       1000,-1000.015,       0.01,          0,   0.005;...       
%       1000,      0.1,  -1000.201,        0.1,   0.001;...   % 10000 [sec^-1]
%       1000,        0,       0.01,  -1000.015,   0.005;...
%          0,      200,       1000,        200,   -1400];     % 10, 50, % acceptor blinking rates [sec^-1]

     
% photophysicalRateMatrix = [...
%     -10000,    10000,       0;...   % donor blinking rates [sec^-1]
%         90,     -100,      10;...   % 10000 [sec^-1]
%          0,     5000,   -5000];     % 10, 50, % acceptor blinking rates [sec^-1]

conformationalRateMatrix = [...
    -1000, 1000;...
    1000, -1000];




mixProbs = [...
    0.0;...
    1.0];




sStateMeans{1} = [...
    0.70;... % S = 0 --> Donor blinking ; S = 1 --> Acceptor blinking
    0.70];

eStateMeans{1} = [...
    0.5, 0.5;...
    0.5, 0.5]; % Needs to be a numEStates X numSStates matrix

sInitialProbabilities{1} = [...
    0.5;...
    0.5];

eInitialProbabilities{1} = [...
    0.5;...
    0.5];

sRateMatrix{1} = [...
    -1000,    1000;...
     1000,   -1000];
 
eRateMatrix{1} = [...
    -1000, 1000;...
    1000, -1000];


 
 


sStateMeans{2} = [...
    0.30;... 0.70
    0.70]; % 1.00

eStateMeans{2} = [...
    0.5, 0.0;...
    0.5, 0.0];


sInitialProbabilities{2} = [...
    0.5;... 0.9
    0.5]; % 0.1


eInitialProbabilities{2} = [...
    0.5;...
    0.5];


sRateMatrix{2} = [...
    -1000,    1000;...
     1000,   -1000];
 
 
eRateMatrix{2} = [...
    -1000, 1000;...
    1000, -1000];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Waitbar Config:

S = what('Seneca');
        
PTH = [S.path, '/matlab/Images'];  

IMS = {...
    '/ALEx1.png',...
    '/ALEx2.png',...
    '/ALEx3a.png',...
    '/ALEx3b.png',...
    '/ALEx4.png',...
    '/ALEx5.png',...
    '/ALEx6.png',...
    '/ALEx7.png'};

STACK = [];

% Cycle through each image and preserve transparency:
for k = 1 : 8


    FULLPTH = [PTH, IMS{k}];

%             IMAGE = imread(FULLPTH, 'backgroundcolor', [1,1,1]);

    [IMAGE, MAP, ALPHA] = imread(FULLPTH);

    IMSIZE = size(IMAGE);

    STACK(1:IMSIZE(1), 1:IMSIZE(2), 1:IMSIZE(3), k) = double(IMAGE) / 255;
    ALPHA1(1:IMSIZE(1), 1:IMSIZE(2), k) = ALPHA;

end % end FOR

STACK(STACK == 0) = NaN;
        
% Assemble animation parameters:
MyAnimation.Animation = STACK;
MyAnimation.Alpha = ALPHA1;
MyAnimation.Background = 'none';


ANIMATION.SEQUENCES{1} = [1; 2; 3];
ANIMATION.SEQUENCES{2} = [1; 2; 4; 5];
ANIMATION.SEQUENCES{3} = [6; 7; 8];

ANIMATION.handle = AnimationProgressBar('initialize',...
    'tag', {'Bar1'},...
    'animation', MyAnimation,...
    'process', {...
        'Getting Arrival Matrix',...
        'Winnowing Data',...
        'Filtering Data',...
        'BCA',...
        'Cluster Statistics'},...
    'title', {'Burst Correlation Analysis (recoloring)'},...
    'color', {[0 0 1; 0 1 0; 1 0 0]},...
    'Pointer', 'hourglass',...
    'info', sprintf('Beginning Burst Correlation Analysis (recoloring)...'));

ANIMATION.SEQFLG = 1;
ANIMATION.SEQ = ANIMATION.SEQUENCES{ANIMATION.SEQFLG};
ANIMATION.SEQLEN = numel(ANIMATION.SEQ);

ANIMATION.EXCITONS = 4;

ANIMATION.Counter = 1;
ANIMATION.Exciton = 1;
ANIMATION.ExMode = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Burst Correlation Analysis:

tic

ArrivalMatrix = [];

%%% TEMP!
appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All = [];
appData.CData.Concatenated.Bursts.EndTime.AllDataSetsValuePerBurst.Unfiltered.All = [];
%%%

TOTALFILES = numel(appData.Filepaths);
for NumFile = 1 : TOTALFILES
    
    FLG{1} = AnimationProgressBar_Update(...
        'Update', ANIMATION.handle,...
        'tag', 'Bar1',...
        'animation', ANIMATION.SEQ(ANIMATION.Counter),...SEQUENCE(FRAME),...
        'value', NumFile/TOTALFILES,...
        'message', '',...
        'process', 'Getting Arrival Matrix',...
        'info', sprintf('Obtaining arrival matrix for file %d', NumFile));

    % ------------------------------------------------------------------- %
    %% Arrival Matrix:
    
    % Get the complete filename:
    FileName = appData.Filepaths{NumFile};
    
    % Convert the *.sm file to a *.txt file, generating an output
    % arrival matrix:
    try
        Temp = sm2txtfast(char(FileName));
        
        if NumFile > 1
            LatestArrival = ALTERNATIONPERIOD * ceil(ArrivalMatrix(end,1) / ALTERNATIONPERIOD); % start afresh, as if the previous file ended perfectly at the end of an alternation period (Gex,Rex)

            Temp(:,1) = Temp(:,1) + LatestArrival;
        end % end IF
        
        ArrivalMatrix = [ArrivalMatrix; Temp];
        
    catch ME
        
        try 
            Temp = sm2txt(char(FileName)); 
            
            if NumFile > 1
                LatestArrival = ALTERNATIONPERIOD * ceil(ArrivalMatrix(end,1) / ALTERNATIONPERIOD); % start afresh, as if the previous file ended perfectly at the end of an alternation period (Gex,Rex)

                Temp(:,1) = Temp(:,1) + LatestArrival;
            end % end IF
            
            ArrivalMatrix = [ArrivalMatrix; Temp];
            
        catch ME2
            keyboard
            ERRORFLG = true;
            return % You will error out...
        end 
        
    end % end TRY/CATCH

%     ArrivalMatrix = load(char(FileName));
    
    %%% HACK!!!! This is necessary due to each photon stream beginning
    %%% again from t = 0 sec; hence, for concatenation, we need to create a
    %%% single time vector, as if only a single stream comprised off all
    %%% data series...
    StartTimes = appData.Data(NumFile).Bursts.StartTime.Unfiltered.All(:,1);
        
    EndTimes = appData.Data(NumFile).Bursts.EndTime.Unfiltered.All(:,1);
        
    
    if NumFile > 1
        StartTimes = StartTimes + LatestArrival;        
        EndTimes = EndTimes + LatestArrival;
        
    end % end IF
    
    appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All =...
        [appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All;...
        StartTimes];
    
    appData.CData.Concatenated.Bursts.EndTime.AllDataSetsValuePerBurst.Unfiltered.All =...
        [appData.CData.Concatenated.Bursts.EndTime.AllDataSetsValuePerBurst.Unfiltered.All;...
        EndTimes];


    %%% END HACK!!!!
    
end % end FOR





% ------------------------------------------------------------------- %
%% Winnow the Data:

FLG{1} = AnimationProgressBar_Update(...
    'Update', ANIMATION.handle,...
    'tag', 'Bar1',...
    'animation', ANIMATION.SEQ(ANIMATION.Counter),...SEQUENCE(FRAME),...
    'value', 0,...
    'message', '',...
    'process', 'Winnowing Data',...
    'info', sprintf('Winnowing data from compiled arrival matrix'));


% Grab the total number of bursts:
TOTALBURSTS = numel(appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All(:,1));

% Create an index array:
BURSTINDICES = 1 : TOTALBURSTS;

% Check!!!! Before binning photons, need to ensure that no StartTimes
% and EndTimes are equal (i.e., no burst ends when another begins), as
% this creates problems when creating the edges (eliminates a non-burst
% bin).
AdjacentBursts = intersect_sorted(...
    appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All(:,1),...
    appData.CData.Concatenated.Bursts.EndTime.AllDataSetsValuePerBurst.Unfiltered.All(:,1));

if isempty(AdjacentBursts) == false
    % Adjacent bursts have been found!!!
    fprintf('%d bursts were found to be immediately adjacent.\n',...
        numel(AdjacentBursts));

    % Cycle through each Adjacent Burst:
    for ADJBURST = 1 : numel(AdjacentBursts)

        % Find the burst:
        AdjBurstIndex = find(...
            appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All(:,1) == AdjacentBursts(ADJBURST),...
            1,...
            'first');

        % Increase the arrival time of any photons which
        % arrived at the previous start time by a fractional amount:
        ArrivalMatrix(ArrivalMatrix(:,1) == appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All(AdjBurstIndex, 1), 1) =...
            ArrivalMatrix(ArrivalMatrix(:,1) == appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All(AdjBurstIndex, 1), 1)...
            +...
            1e-10; % should be less than NI card time resolution...

        % Increase the start time of the next burst by a fractional amount:
        appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All(AdjBurstIndex, 1) =...
            appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All(AdjBurstIndex, 1)...
            +...
            1e-10;

    end % end FOR


end % end IF


% Create Edges:
Edges = union_sorted_rows(...
    appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All(:,1),...
    appData.CData.Concatenated.Bursts.EndTime.AllDataSetsValuePerBurst.Unfiltered.All(:,1));

% Histogram the arrival times:
[Counts, BinIndices] = histc(ArrivalMatrix(:,1), Edges);

% Determine which bins to throw away: (e.g., all those bins containing
% photons which do not belong to a burst)

% Discard the even BinIndices, as these correspond to bins defined by
% EdgeA = EndTime(1); EdgeB = StartTime(2). Hence, this bin (EdgeA,
% EdgeB) corresponds to photons which do not belong to a burst (i.e., photons
% occurring after the end of the previous burst and before the start of
% the next burst):
VALS1 = 0:2:max(BinIndices); % Note: we include 0 to account for photons not falling into any defined bin.
TF = ismember_sorted(BinIndices, VALS1);




% ----- %
%% Background Calculations:


% Before discarding the even BinIndices, let us calculate the background
% count rates in the donor and acceptor channels:
BkgdArray(:,1) = [appData.CData.Concatenated.Bursts.StartTime.AllDataSetsValuePerBurst.Unfiltered.All(:,1); ArrivalMatrix(end,1)];
BkgdArray(:,2) = -1*[0; appData.CData.Concatenated.Bursts.EndTime.AllDataSetsValuePerBurst.Unfiltered.All(:,1)];
% BkgdArray = [ StartOfNextBurst | -1*EndOfPreviousBurst ]

BkgdArray = sum(BkgdArray, 2); % sum across columns; % Time between successive bursts

Bkgd.Time.Total = sum(BkgdArray); % Total time in which we record background counts
Bkgd.Time.Gex = Bkgd.Time.Total / 2; % [sec]
Bkgd.Time.Rex = Bkgd.Time.Gex; % Assume number of alternation periods for Gex and Rex are equal; In the limit of long data sets, this is true, given AltPeriod <<<<<< T

clear BkgdArray

BkgdPhotons = ArrivalMatrix(TF == true, :);

% Separate the backgrond photons into separate photon streams:
Bkgd.DexDem.ArrivalMatrix = BkgdPhotons(BkgdPhotons(:,2) == 0, :);
Bkgd.DexAem.ArrivalMatrix = BkgdPhotons(BkgdPhotons(:,2) == 2, :);
Bkgd.AexDem.ArrivalMatrix = BkgdPhotons(BkgdPhotons(:,2) == 3, :);
Bkgd.AexAem.ArrivalMatrix = BkgdPhotons(BkgdPhotons(:,2) == 1, :);


% Calculate the background count rates:
Bkgd.Rates(1,1)        = size(Bkgd.DexDem.ArrivalMatrix,1) / Bkgd.Time.Gex; % [sec^-1]
Bkgd.Rates(2,1)        = size(Bkgd.DexAem.ArrivalMatrix,1) / Bkgd.Time.Gex;
Bkgd.Rates(3,1)        = size(Bkgd.AexDem.ArrivalMatrix,1) / Bkgd.Time.Rex;
Bkgd.Rates(4,1)        = size(Bkgd.AexAem.ArrivalMatrix,1) / Bkgd.Time.Rex;
Bkgd.totalRate = size(BkgdPhotons,1) / Bkgd.Time.Total;

% Bkgd.Rates(1,1)        = 0;%size(Bkgd.DexDem.ArrivalMatrix,1) / Bkgd.Time.Gex; % [sec^-1]
% Bkgd.Rates(2,1)        = 0;%size(Bkgd.DexAem.ArrivalMatrix,1) / Bkgd.Time.Gex;
% Bkgd.Rates(3,1)        = 0;%size(Bkgd.AexDem.ArrivalMatrix,1) / Bkgd.Time.Rex;
% Bkgd.Rates(4,1)        = 0;%size(Bkgd.AexAem.ArrivalMatrix,1) / Bkgd.Time.Rex;
% Bkgd.totalRate = 0;%size(BkgdPhotons,1) / Bkgd.Time.Total;


% --------- %


% Remove all photons not belonging to a burst, thus creating a 'burst
% arrival matrix':
ArrivalMatrix(TF == true,:) = [];
BinIndices(TF == true) = [];    


% ------------------------------------------------------------------- %
%% Filtering:

FLG{1} = AnimationProgressBar_Update(...
    'Update', ANIMATION.handle,...
    'tag', 'Bar1',...
    'animation', ANIMATION.SEQ(ANIMATION.Counter),...SEQUENCE(FRAME),...
    'value', 0,...
    'message', '',...
    'process', 'Filtering Data',...
    'info', sprintf('Filtering data from compiled arrival matrix'));



% Filter the data:
[appData, REMOVEBURSTS] = Confocal2ColorALExSolnMultDataBVA_FilterData(...
    appData);

% Remove the filtered bursts:
TF = ismember_sorted((BinIndices + 1) ./ 2, REMOVEBURSTS);

ArrivalMatrix(TF == true, :) = [];
BinIndices(TF == true, :) = [];


% For Burst Correlation Analysis, we need to retain all photons, both those
% originating from donor excitation as well as acceptor excitation.

% % Now, we need to rid ourselves of all photons originating from red
% % excitation: DexDem = 0; DexAem = 2; AexAem = 1; AexDem = 3; => all
% % photons labeled 1 or 3... (NOTE: this also removes acceptor-only
% % bursts!)
% LogicalIndexArray = (ArrivalMatrix(:,2) == 1 | ArrivalMatrix(:,2) == 3);
% VALS2 = unique(BinIndices(LogicalIndexArray, :)); % This gives all burst indices which contain acceptor excitation photons
% 
% ArrivalMatrix(LogicalIndexArray,:) = []; % All acceptor excitation photons removed.
% BinIndices(LogicalIndexArray,:) = []; % Note: this array will be important now, as all photons can be grouped by their unique identification with a Bin...
% % NOTE: now VALS and unique(BinIndices) do not necessarily
% % correspond!!!
% 
% % Determine which, if any, bursts have been removed from further
% % analysis:
% REMOVEBURSTS = [REMOVEBURSTS; (setdiff_sorted(VALS2, unique(BinIndices)) + 1) ./ 2];

% REMOVEBURSTS now contains all bursts which either did not meet the
% filtering criteria above or were acceptor-only...

% Frequently, bursts have values which are not numbers; find them:
NaNBURSTS = find(isnan(appData.CData.Concatenated.Bursts.Sraw.AllDataSetsValuePerBurst.Unfiltered.All(:,1)));

REMOVEBURSTS = [REMOVEBURSTS; NaNBURSTS];

% ------------------------------------------------------------------- %   
% Clustering: cluster the bursts according to mean stoichiometry values:

% Define the cluster bounds: % NOTE:this creates an extra bin!!!
CLUSTEREDGES = -CLUSTERWIDTH/2:CLUSTERWIDTH:1+CLUSTERWIDTH/2;%0:CLUSTERWIDTH:1;
CLUSTEREDGES = CLUSTEREDGES';


% ------------------------------------------------------------------- %

% Split the Arrival Matrix into Bursts:
IndexCell = SplitVec(...
    BinIndices,...
    'equal',... % Group all elements which are equal
    'loc'); % Return the locations of these elements in a cellular array

% Each cellular element of IndexCell corresponds to a burst; within each
% cell element, the individual photons are assigned a location
% corresponding to their index in the (filtered) ArrivalMatrix.


% ------------------------------------------------------------------- %
%% BCA:

FLG{1} = AnimationProgressBar_Update(...
    'Update', ANIMATION.handle,...
    'tag', 'Bar1',...
    'animation', ANIMATION.SEQ(ANIMATION.Counter),...SEQUENCE(FRAME),...
    'value', 0,...
    'message', '',...
    'process', 'BCA',...
    'info', sprintf('Beginning sliding window computations...'));



% Initialize a BurstCorrelationAnalysis array:
appData.BCA.Data.Results = nan(TOTALBURSTS, 15);

% Initialize a BurstVarianceAnalysis array:
appData.BVA.Data.Results = nan(TOTALBURSTS, 2);

% Loop through the bursts and perform BCA:
COUNTER = 0;
CLUSTERSTATS = cell(CLUSTERBINS+1, 1);

TRACKER = 0;
for BURST = 1 : TOTALBURSTS

    % Check!!!
    if ismember(BURST, REMOVEBURSTS) % Determine if BURST is blacklisted...
        % Blacklisted burst found! Move along to next burst...
        continue;
    else
        % Update our counter for our non-blacklisted bursts...
        COUNTER = COUNTER + 1;
    end % end IF

    % Grab the current (non-blacklisted) BURST indices:
    INDICES = IndexCell{COUNTER};

    % Remove a 'cutoff' number of photons:
    if numel(INDICES) < PHOTONCUTOFF
        % Go to next iteration...
        REMOVEBURSTS(end+1,1) = BURST;
        TRACKER = TRACKER + 1;
        continue;
    else
%         INDICES(1:PHOTONCUTOFF) = []; % Not strictly necessary due to
%         rounding used below, as likelihood that a photon arrival time
%         falls exactly on laser excitation time is exceedingly low,
%         particularly given the, say, 3us deadtime to prevent cross talk
%         between laser excitations.
    end % end IF/ELSE      
    
    
    % ------------------------------------------------------------------- %
    % New Code! (2012-08-09) - Multiple sub-populations within sample
    
    
    % Determine to which sub-population this burst belongs:
    randVal = rand;
    
    randValVec = randVal * ones(numel(mixProbs),1); 

    % Generate a cumulative histogram from the bkgd rates:
    cumProbHist = cumsum(mixProbs ./ sum(mixProbs));

    % Find where our random value resides among the classes:
    classy = cumProbHist - randValVec;

    % Determine to what photon ID to assign the random value:
    mixID = find(classy > 0, 1, 'first');

    photophysicalInitialProbabilities = sInitialProbabilities{mixID};  
    conformationalInitialProbabilities = eInitialProbabilities{mixID};
    photophysicalRateMatrix = sRateMatrix{mixID};
    conformationalRateMatrix = eRateMatrix{mixID};
    efficiencyStateMeans = eStateMeans{mixID};
    stoichiometryStateMeans = sStateMeans{mixID};
    
    
    
    
    % ------------------------------------------------------------------- %
    % New Code! (2012-08-09)
    
    % Calculate total number of photons:
    TOTALPHOTONS = numel(INDICES);        

    % Check!!!
    if TOTALPHOTONS < WINDOW %NUMPHOTONS
        % Go to next iteration, as not enough photons in 'burst'...
        REMOVEBURSTS(end+1,1) = BURST;
        TRACKER = TRACKER + 1;
        continue;
    end % end IF

    % Abstract a Burst Photon Arrival Matrix from the Arrival Matrix:
    BurstPhotonArrivalMatrix = ArrivalMatrix(INDICES, :); % 2 columns: Arrival Times | Photons
    
    BurstPhotonArrivalMatrix = [BurstPhotonArrivalMatrix, nan(size(BurstPhotonArrivalMatrix,1), 2), zeros(size(BurstPhotonArrivalMatrix,1), 1)]; % tack on an extra three columns for our photophysical and conformational states and background ID
    
    
    
    
    % [1] Assign each photon to photophysical and conformational state from appropriate Markov chains:
        
    
    
    % [1a] Run a Markov Chain simulation for the burst duration, output an
    % array with state ID and absolute transition time: (photophysical
    % states)
    photophysicalMarkovChain = continuousTimeKineticMonteCarlo(...
        [],...
        BurstPhotonArrivalMatrix(end,1) - BurstPhotonArrivalMatrix(1,1),... % total simulation time
        photophysicalInitialProbabilities,...
        photophysicalRateMatrix);
    
    % Set the state transition times according to the burst arrival times:
    photophysicalMarkovChain(:,2) = photophysicalMarkovChain(:,2) + BurstPhotonArrivalMatrix(1,1); 
    
    
    
    
    % [1b] Histogram the photons:
    [temp, IDX] = histc(BurstPhotonArrivalMatrix(:,1), [0; photophysicalMarkovChain(:,2)]); % the last transition in the Markov chain should exceed the last photon arrival
    
    BurstPhotonArrivalMatrix(:,3) = photophysicalMarkovChain(IDX,1); % assigns the relevant states based on the bin index
    
    
    
    
    % [1c] Run a Markov Chain simulation for the burst duration, output an
    % array with state ID and absolute transition time: (conformational
    % states)
    conformationalMarkovChain = continuousTimeKineticMonteCarlo(...
        [],...
        BurstPhotonArrivalMatrix(end,1) - BurstPhotonArrivalMatrix(1,1),... % total simulation time
        conformationalInitialProbabilities,...
        conformationalRateMatrix);
    
    % Set the state transition times according to the burst arrival times:
    conformationalMarkovChain(:,2) = conformationalMarkovChain(:,2) + BurstPhotonArrivalMatrix(1,1); 
    
    
    
    
    % [1d] Histogram the photons:
    [temp, IDX] = histc(BurstPhotonArrivalMatrix(:,1), [0; conformationalMarkovChain(:,2)]); % the last transition in the Markov chain should exceed the last photon arrival
    
    BurstPhotonArrivalMatrix(:,4) = conformationalMarkovChain(IDX,1); % assigns the relevant states based on the bin index
    
    
    
    
    
   
    
    
    % [2] Separate, condense, and merge the Gex and Rex photon streams,
    % retaining the state information.
    
    % Create an edge vector to time window the Burst Photon Arrival Matrix:
    StartEdge = ALTERNATIONTIME * ceil(BurstPhotonArrivalMatrix(1,1) / ALTERNATIONTIME);
    EndEdge = ALTERNATIONTIME * floor(BurstPhotonArrivalMatrix(end,1) / ALTERNATIONTIME);
    
    edgesAlt = StartEdge : ALTERNATIONTIME : EndEdge; 
    
    if rem(edgesAlt, 2)
        edgesAlt(end) = []; % remove the last edge to ensure even number of bins
    end % end IF
    
    % ---------- %
    % Grab the Dex and Aex photons:
    Dex = BurstPhotonArrivalMatrix(BurstPhotonArrivalMatrix(:,2) == 0 | BurstPhotonArrivalMatrix(:,2) == 2,:); % DexDem | DexAem
    Rex = BurstPhotonArrivalMatrix(BurstPhotonArrivalMatrix(:,2) == 1 | BurstPhotonArrivalMatrix(:,2) == 3,:); % AexDem | AexAem
    
    % Bin the Dex and Aex photons by alternation time:
    [Counts, DexIDX] = histc(Dex(:,1), edgesAlt);
    [Counts, RexIDX] = histc(Rex(:,1), edgesAlt);
    
    % Remove photons not in a bin:
    Dex(DexIDX == 0,:) = [];
    Rex(RexIDX == 0,:) = [];
    DexIDX(DexIDX == 0) = [];
    RexIDX(RexIDX == 0) = [];
    
    
    
    % [!] Determine the number of background photons, before condensing and
    % merging...
    bkgdPhotons = poissrnd(Bkgd.totalRate * (max([Dex(end,1), Rex(end,1)])-min([Dex(1,1), Rex(1,1)]))); % \lambda * (\delta T)
        
    
    
    
    % Determine which is first Gex or Rex:
    if (sum(mod(DexIDX,2) == 1) >= 2) && (sum(mod(RexIDX,2) == 0) >= 2) % at least 2 alternation cycles
        altFLG = 0; % Gex in first alternation
        if sum(mod(DexIDX,2) == 0)
            %keyboard; % all bins should be odd!
%             fprintf('%d Dex photon(s) in wrong bin. Burst No: %d\n', sum(mod(DexIDX,2) == 0), BURST);
            Dex(mod(DexIDX,2) == 0,:) = [];
            DexIDX(mod(DexIDX,2) == 0,:) = [];            
        end % end IF
        if sum(mod(RexIDX,2) == 1)
            %keyboard; % all bins should be even!
%             fprintf('%d Rex photon(s) in wrong bin. Burst No: %d\n', sum(mod(RexIDX,2) == 1), BURST);
            Rex(mod(RexIDX,2) == 1,:) = [];
            RexIDX(mod(RexIDX,2) == 1,:) = [];
        end % end IF
    elseif (sum(mod(DexIDX,2) == 0) >= 2) && (sum(mod(RexIDX,2) == 1) >= 2)
        altFLG = 1; % Rex
        if sum(mod(DexIDX,2) == 1)
            %keyboard; % all bins should be even!
%             fprintf('%d Dex photon(s) in wrong bin. Burst No: %d\n', sum(mod(DexIDX,2) == 1), BURST);
            Dex(mod(DexIDX,2) == 1,:) = [];
            DexIDX(mod(DexIDX,2) == 1,:) = [];
        end % end IF
        if sum(mod(RexIDX,2) == 0)
            %keyboard; % all bins should be odd!
%             fprintf('%d Rex photon(s) in wrong bin. Burst No: %d\n', sum(mod(RexIDX,2) == 0), BURST);
            Rex(mod(RexIDX,2) == 0,:) = [];
            RexIDX(mod(RexIDX,2) == 0,:) = [];
        end % end IF
    else
        keyboard; % too few alternation cycles
    end % end IF/ELSE
    
    
    
    % [2a] Loop through each bin and re-assign (artifically move) the photon
    % arrival times: (i.e., we effectively create two continuous excitation
    % streams)
    altBin = 1;
    for k = 1 : 2 : (numel(Counts)-1)
        
        % Two cases:
        if ~altFLG % Dex is first
            Dex(DexIDX == k, 1) = Dex(DexIDX == k, 1) - ((altBin-1)*ALTERNATIONTIME);
            Rex(RexIDX == k+1, 1) = Rex(RexIDX == k+1, 1) - (altBin*ALTERNATIONTIME);
            
        else % Rex is first
            Dex(DexIDX == k+1, 1) = Dex(DexIDX == k+1, 1) - (altBin*ALTERNATIONTIME);
            Rex(RexIDX == k, 1) = Rex(RexIDX == k, 1) - ((altBin-1)*ALTERNATIONTIME);
                      
        end % end IF/ELSE
        
        % Update our counter:
        altBin = altBin + 1;
        
    end % end FOR altBin
    
    
    
    % [2b] With re-assigned photon arrival times, concatenate Gex and Rex
    % photons and sort, creating a merged Poisson process:
    mergedPhotons = [Dex; Rex];
    
    [sortedVec, IDX] = sort(mergedPhotons(:,1));
    
    mergedPhotons = mergedPhotons(IDX, :);
    
    
    
    
    
    % [3] Re-color the photons based on probability S, as a function of a
    % Markov chain. Re-color recolored donor excitation photons based on
    % probability E, as a function of a (conformational) Markov chain.
    recoloredPhotons = mergedPhotons;
    numPhotons = size(recoloredPhotons,1);
    probBkgd = bkgdPhotons / numPhotons; 
    onesVec = ones(4, 1); % 4 photon IDs
    photon = 1;
    while (photon <= numPhotons)
        
        % Determine if photon is background or from process:
        randVal = rand;
        
        if randVal < probBkgd
            % Background photon:
            
            randVal = rand;
            
            randValVec = randVal * onesVec; 
            
            % Generate a cumulative histogram from the bkgd rates:
            cumProbHist = cumsum(Bkgd.Rates ./ sum(Bkgd.Rates));

            % Find where our random value resides among the classes:
            classy = cumProbHist - randValVec;

            % Determine to what photon ID to assign the random value:
            photonID = find(classy > 0, 1, 'first');

            % Update our recolored photon array:
            switch photonID                
                case 1
                    recoloredPhotons(photon,2) = 0; % DexDem
                case 2
                    recoloredPhotons(photon,2) = 2; % DexAem
                case 3
                    recoloredPhotons(photon,2) = 3; % AexDem
                case 4
                    recoloredPhotons(photon,2) = 1; % AexAem
            end % end SWITCH photonID            
            
            % Flag that this photon is background:
            recoloredPhotons(photon,5) = 1; % otherwise 0
            
        else
            % Process photon:
        
            % Draw a random number:
            randVal = rand;

            % Determine if "green" or "black" (i.e., Gex or Rex photon):
            if randVal < stoichiometryStateMeans(mergedPhotons(photon,3))
                recoloredPhotons(photon,2) = 0; % DexDem

                % Determine if "green" photon will be colored "red":
                randVal = rand;
                if randVal < efficiencyStateMeans(mergedPhotons(photon,4), mergedPhotons(photon,3)) % column vector corresponds to the E states associated with a particular S state
                    
                    % Determine if we should 'discard' the photon due to
                    % gamma: (basically, throw the photon back in the pool
                    % to be reassigned; if gamma = 1, we never lose any
                    % energy transferred photons)
                    randVal = rand;
                    if randVal < gamma 
                        recoloredPhotons(photon,2) = 2; % DexAem
                    else
                        continue; % allow the photon to be reassigned by redoing this iteration
                    end % end IF/ELSE
                end % end IF

            else
                recoloredPhotons(photon,2) = 1; % AexAem
            end % end IF/ELSE
            
        end % end IF/ELSE
        
        % Move to next photon:
        photon = photon + 1;
        
    end % end WHILE photon<=numPhotons
    
    
    
    
    
    % [4] Slide a window W and calculate E and S:
    numWindows = floor(size(recoloredPhotons(:,1),1) / WINDOW);% floor(sum((recoloredPhotons(:,2)==0) + (recoloredPhotons(:,2)==2)) / WINDOW); % number of recolored Gex photons
    EfficiencyVector = nan(numWindows, 1);
    StoichiometryVector = nan(numWindows, 1);
    EfficiencyPhotons = nan(numWindows, 1);
    
    temp = recoloredPhotons;
    for win = 1 : numWindows % non-overlapping windows 
        
        % Get the photon streams in the window:
        DexDem = sum(temp(1:WINDOW,2) == 0);
        DexAem = sum(temp(1:WINDOW,2) == 2);
        AexDem = sum(temp(1:WINDOW,2) == 3);
        AexAem = sum(temp(1:WINDOW,2) == 1);
        
        % Calculate E and S:
        EfficiencyVector(win) = DexAem ./ (DexDem + DexAem);
        StoichiometryVector(win) = (DexDem + DexAem) ./ (DexDem + DexAem + AexAem + AexDem);
        
        EfficiencyPhotons(win) = DexDem + DexAem;
        
        % Remove the photons from the stream:
        temp(1:WINDOW,:) = []; % recoloredPhotons(1:IDS(end),:) = []; % makes the array smaller at each iteration, thus allowing us to index from 1 as done above.
                
    end % end FOR win
    
    
    
    
    % ------------------------------------------------------------------- %
    % End
    
    % [4] Calculate the E and S means and standard deviations:
    EfficiencyStDev = nanstd(EfficiencyVector);
    EfficiencyMean = nanmean(EfficiencyVector);
    StoichiometryStDev = nanstd(StoichiometryVector);
    StoichiometryMean = nanmean(StoichiometryVector);
    
    % Calculate the Linear Correlation Coefficient:
    corrWindows = floor(numWindows / CORRWINDOW); % overlapping windows/non-overlapping windows --> numWindows - CORRWINDOW + 1; %
    
    if corrWindows < 1
        % not enough E-S pairs for correlation:
        REMOVEBURSTS(end+1,1) = BURST;
        TRACKER = TRACKER + 1;
        continue;
    end % end IF
    
    rho = nan(corrWindows,1);
    IDX = 1;
    for win = 1 : corrWindows
        
        E = EfficiencyVector(IDX:IDX+CORRWINDOW-1);
        S = StoichiometryVector(IDX:IDX+CORRWINDOW-1);
        rho(win) = correlationCoefficient(E, S);
        
        IDX = IDX + CORRWINDOW; % + 1; % 
        
    end % end FOR win
    
    % Assign over the windowed correlation coefficients:
    winCorrs = rho;
    
    % Calculate the average correlation coefficient:
    rho = nanmean(rho); % sometimes have missing values due to no FRET photons
                    
    % [5] Permutation test for significance:
    numPerms = 120; % number of repeats 
    permTest = nan(numPerms,1);
    for i = 1 : numPerms
        
        % Mix up the indexes:
        indexes = randperm(numWindows); % numel(EfficiencyVector)

        % Get the data permutation:
        xPerm = EfficiencyVector(indexes);
        
        permRho = nan(corrWindows,1);
        IDX = 1;
        for win = 1 : corrWindows
            
            E = xPerm(IDX:IDX+CORRWINDOW-1);
            S = StoichiometryVector(IDX:IDX+CORRWINDOW-1);
            
            % Calculate the correlation coefficient:
            permRho(win) = correlationCoefficient(E, S); 
            
            IDX = IDX + CORRWINDOW; % + 1;
            
        end % end FOR win
        
        % Calculate the average correlation coefficient after permutation:
        permTest(i) = mean(permRho);
        
    end % end FOR i
    
    % Store the pValues:
    pVal = sum(abs(permTest) > abs(rho)) / numPerms; % percentage
    
    
    % [5] Permutation test for significance: 
    numPerms = 120; % number of repeats
    IDX = 1;
    pValues = nan(corrWindows,1);
    for win = 1 : corrWindows
        
        % Grab the E/S data:
        E = EfficiencyVector(IDX:IDX+CORRWINDOW-1);
        S = StoichiometryVector(IDX:IDX+CORRWINDOW-1);
        
        permRho = nan(numPerms,1);
        for i = 1 : numPerms
            
            % Mix up the indexes:
            indexes = randperm(CORRWINDOW); % num data pts (E/S pairs) within a correlation window (e.g., 10)

            % Get the data permutation:
            xPerm = E(indexes);
        
            % Calculate the correlation coefficient:
            permRho(i) = correlationCoefficient(xPerm, S);
            
        end % end FOR i
        
        % Determine the p value:
        pValues(win) = sum(abs(permRho) >= abs(winCorrs(win))) / numPerms; % percentage
        
        IDX = IDX + CORRWINDOW; % 1; %
        
    end % end FOR win

    % Put the pValues and correlations in an array:
    corrMat = [winCorrs, pValues];
    
    % Find the positive and negative correlations:
    sigMat = nan(size(corrMat,1),1);    
    sigMat(corrMat(:,1) < 0, 1) = -1; % all negative correlations
    sigMat(corrMat(:,1) > 0, 1) = +1; % positive correlations
    
    % Find the pValues which are greater than the alpha: (i.e. not
    % significant)
    alpha = 0.05;
    sigMat(corrMat(:,2) > alpha, 1) = 0; % 
    
    % Left with a vector in which significant correlations are
    % non-zero and given by their sign:
    sigMat(isnan(sigMat),1) = 0; % due to searching for S photons, we may have 0 Gex photons, in which case winCorr is 0 as is our pValue. Here, the correlation is not significant.
    
    % Split the vector into `dwells' and find the length of each dwell:
    runs = SplitVec(...
        sigMat,...
        'equal',...
        'length');
    
    vals = SplitVec(...
        sigMat,...
        'equal',...
        'firstval');
    
    % Find number of correlation hot spots beyond our run threshold:
    minNeighbors = 3;
    hotSpots = sum(runs(vals~=0) >= minNeighbors);
    
    % Determine if the number of hot spots (regions of high and similar
    % correlation) meet our minimum threshold:
    minHotSpots = 3; % meaning three regions of high correlation
     
    if hotSpots >= minHotSpots
        % Found multiple hot spots:
        burstFLG = 1;
    else
        burstFLG = 0;
    end % end IF/ELSE
    
    
    % Alternative method: method of fixed effects for testing significance
    % of average correlation coefficients (Hedges & Olkin, 1985):
    zScore = rho ./ sqrt(1./ (corrWindows*(CORRWINDOW-3))); % rho / standardError
    pVal2 = normcdf(-abs(zScore),0,1);
    
    
    
    % Determine the number of states:
    numStates = [numel(unique(mergedPhotons(:,3))), numel(unique(mergedPhotons(:,4)))]; % S | E
     
    % Determine the number of transitions:
    numTransitions = [sum(diff(mergedPhotons(:,3))~=0), sum(diff(mergedPhotons(:,4))~=0)]; % S | E

    % Place the mean and stdev of the stoichiometry in our BCA array:
    appData.BCA.Data.Results(BURST, :) = [...
        EfficiencyMean,...
        EfficiencyStDev,...
        StoichiometryMean,...
        StoichiometryStDev,...
        rho,...
        pVal,...
        numWindows,...
        corrWindows,...
        pVal2,...
        numStates,.... % 1x2
        numTransitions,...]; % 1x2
        burstFLG,...
        hotSpots];
        
    
    
    % ------------------- %
    %% S-BVA:
    
    
    if filterFLG
    
        % [*] Permutation test for significance: (each window)
        corrWindows = floor(numWindows / CORRWINDOW);

        if corrWindows < 1
            % not enough correlation windows:
            EfficiencyVector = nan;
            StoichiometryVector = nan;
            break;
        end % end IF
        
        % Get the mean values:
        meanE = nanmean(EfficiencyVector); % for entire burst; assume static species
        meanS = nanmean(StoichiometryVector);

        numPerms = 120; % number of repeats
        IDX = 1;
        pValues = nan(corrWindows,1);
        for win = 1 : corrWindows

            % Grab the E/S data:
            E = EfficiencyVector(IDX:IDX+CORRWINDOW-1);
            S = StoichiometryVector(IDX:IDX+CORRWINDOW-1);

            permRho = nan(numPerms,1);
            for i = 1 : numPerms

                % Mix up the indexes:
                indexes = randperm(CORRWINDOW); % num data pts (E/S pairs) within a correlation window (e.g., 10)

                % Get the data permutation:
                xPerm = E(indexes);

                % Calculate the correlation coefficient:
                permRho(i) = correlationCoefficient(xPerm, S);

            end % end FOR i

            % Determine the p value:
            pValues(win) = sum(abs(permRho) >= abs(winCorrs(win))) / numPerms; % percentage

            if pValues(win) < sigAlpha
                % We found significant correlation!!!!
                
                % Replace the Efficiency & Stoichiometry Vector values...

                % Draw random numbers from a binomial distribution:
                randVals = binornd(WINDOW, meanS, CORRWINDOW, 1); % col vector                    
                StoichiometryVector(IDX:IDX+CORRWINDOW-1) = randVals ./ WINDOW;

                for i = 1 : CORRWINDOW
                    % Note: the number of Gex photons varies.
                    
                    randVal = binornd(EfficiencyPhotons(IDX+i-1), meanE, 1, 1); % 1x1                   
                    EfficiencyVector(IDX+i-1) = randVal ./ EfficiencyPhotons(IDX+i-1); 
                end 

            end % end IF

            IDX = IDX + CORRWINDOW; % 1; %

        end % end FOR win



        
    end % end IF

    if isnan(StoichiometryVector)
        % Go to next burst, as were not enough photons in 'burst'...
        REMOVEBURSTS(end+1,1) = BURST;
        TRACKER = TRACKER + 1;
        continue;
    end % end IF

   
    
    % [5] Calculate S mean and standard deviation:
    StoichiometryStDev = nanstd(StoichiometryVector);
    StoichiometryMean = nanmean(StoichiometryVector);
    
    % Place the mean and stdev of the stoichiometry in our BVA array:
    appData.BVA.Data.Results(BURST, :) = [...
        StoichiometryMean,...
        StoichiometryStDev];
        
       
    % Histogram the burst by its calculated mean:
    [Counts, WhichCluster] = histc(StoichiometryMean, CLUSTEREDGES);

    % Check!!!
    if WhichCluster == (CLUSTERBINS+1)% +2
        % This accounts for values which are precisely equal to the
        % last edge.  We place this data in the last bin.
        WhichCluster = WhichCluster - 1;
    end % end IF/ELSE

    CLUSTERSTATS{WhichCluster} = [CLUSTERSTATS{WhichCluster}; StoichiometryVector]; % Concatenation
    
    
        
    % ------------------- %
    % Update the waitbar:
    FLG{1} = AnimationProgressBar_Update(...
        'Update', ANIMATION.handle,...
        'tag', 'Bar1',...
        'animation', ANIMATION.SEQ(ANIMATION.Counter),...SEQUENCE(FRAME),...
        'value', BURST/TOTALBURSTS,...
        'message', ['Burst No. ', int2str(BURST)],...
        'process', 'BCA',...
        'info', sprintf('Finishing sliding window calculations for Burst No. %d',...
        	BURST));

    if mod(ANIMATION.Exciton, ANIMATION.EXCITONS+1) == 0
        if ANIMATION.ExMode == 1
            ANIMATION.ExMode = 0;
        else
            ANIMATION.ExMode = 1;
        end % end IF/ELSE

        if ANIMATION.ExMode == 1
            % Donor Excitation:
            CoinFlip = rand;
            if CoinFlip > 0.5
                ANIMATION.SEQFLG = 2;
                ANIMATION.SEQ = ANIMATION.SEQUENCES{ANIMATION.SEQFLG};
                ANIMATION.SEQLEN = numel(ANIMATION.SEQ);
            else
                ANIMATION.SEQFLG = 1;
                ANIMATION.SEQ = ANIMATION.SEQUENCES{ANIMATION.SEQFLG};
                ANIMATION.SEQLEN = numel(ANIMATION.SEQ);
            end % end IF/ELSE

            ANIMATION.Counter = 1;
        else
            % Acceptor Excitation:
            ANIMATION.SEQ = ANIMATION.SEQUENCES{3};
            ANIMATION.SEQLEN = numel(ANIMATION.SEQ);
            ANIMATION.Counter = 1;

        end % end IF/ELSE

    else 
        if ANIMATION.ExMode == 1
            % Donor Excitation
            if ANIMATION.Counter == ANIMATION.SEQLEN
                CoinFlip = rand;
                if CoinFlip > 0.5
                    ANIMATION.SEQFLG = 2;
                    ANIMATION.SEQ = ANIMATION.SEQUENCES{ANIMATION.SEQFLG};
                    ANIMATION.SEQLEN = numel(ANIMATION.SEQ);
                else
                    ANIMATION.SEQFLG = 1;
                    ANIMATION.SEQ = ANIMATION.SEQUENCES{ANIMATION.SEQFLG};
                    ANIMATION.SEQLEN = numel(ANIMATION.SEQ);
                end % end IF/ELSE

                ANIMATION.Counter = 1;
            else
                ANIMATION.Counter = ANIMATION.Counter + 1;
            end % end IF/ELSE

        else
            % Acceptor Excitation:
            if ANIMATION.Counter == ANIMATION.SEQLEN
                ANIMATION.Counter = 1;
            else
                ANIMATION.Counter = ANIMATION.Counter + 1;
            end % end IF/ELSE
        end

    end % end IF

    if ANIMATION.Counter == 1
        ANIMATION.Exciton = ANIMATION.Exciton + 1;
    end 



    EFLAG = sum(strcmpi(FLG, 'Abort'));
    if EFLAG > 0
        return;
    end % end IF
    EFLAG = sum(strcmpi(FLG, 'end'));
    if EFLAG > 0
        break;
    end % end IF


end % end FOR    

fprintf(['%d bursts did not have sufficient photons ',...
    'to perform a sliding window calculation.\n\n'],...
    TRACKER);


IDs = setdiff(BURSTINDICES, REMOVEBURSTS); % Retain only 'good' bursts

appData.BCA.Data.Results = appData.BCA.Data.Results(IDs,:);

% Remove any data in which the mean stoichiometry is either 0 or 1: (why?)
appData.BCA.Data.Results(appData.BCA.Data.Results(:,3) == 0 | appData.BCA.Data.Results(:,3) == 1, :) = [];

% Remove any data in which the standard deviation is zero: (why?)
appData.BCA.Data.Results(appData.BCA.Data.Results(:,4) == 0, :) = [];

% Remove all NaN rows:
appData.BCA.Data.Results(isnan(appData.BCA.Data.Results(:,1)), :) = [];


clear ArrivalMatrix
clear BinIndices





% Tidy-up our Cluster Statistics:
RemoveClusters = cellfun(@isempty, CLUSTERSTATS);
CLUSTERSTATS(RemoveClusters) = [];
CLUSTEREDGES(RemoveClusters) = [];

RemoveClusters = cell2mat(cellfun(@(x) (numel(x)<CLUSTERMINWINDOWS),...
    CLUSTERSTATS,...
    'UniformOutput', false));
CLUSTERSTATS(RemoveClusters) = [];
CLUSTEREDGES(RemoveClusters) = [];


IDs = setdiff(BURSTINDICES, REMOVEBURSTS); % Retain only 'good' bursts

appData.BVA.Data.Results = appData.BVA.Data.Results(IDs,:);

% Remove any data in which the mean stoichiometry is either 0 or 1: (why?)
appData.BVA.Data.Results(appData.BVA.Data.Results(:,1) == 0 | appData.BVA.Data.Results(:,1) == 1, :) = [];

% Remove any data in which the standard deviation is zero: (why?)
appData.BVA.Data.Results(appData.BVA.Data.Results(:,2) == 0, :) = [];

% Remove all NaN rows:
appData.BVA.Data.Results(isnan(appData.BVA.Data.Results(:,1)), :) = [];





% ------------------------------------------------------------------- %
%% Cluster Statistics:

FLG{1} = AnimationProgressBar_Update(...
    'Update', ANIMATION.handle,...
    'tag', 'Bar1',...
    'animation', ANIMATION.SEQ(ANIMATION.Counter),...SEQUENCE(FRAME),...
    'value', 0,...
    'message', '',...
    'process', 'Cluster Statistics',...
    'info', sprintf('Computing the statistics for burst clusters...'));


% ------------ %
% [ ] Group bursts with similar mean values and calculate population
% standard deviation:
ClusterCenters = CLUSTEREDGES(1:end-1) + CLUSTERWIDTH/2;

NumClusters = numel(CLUSTERSTATS);
HighCIvec = zeros(NumClusters, 1);
LowCIvec = zeros(NumClusters, 1);
StDevDistrMean = zeros(NumClusters, 1);

ClusterStats = nan(NumClusters, 2);


% Loop through each 'slice' (or 'cluster' of data) from the 2D Histogram
% and calculate the statistics:
for CLUSTER = 1 : NumClusters
    
    FLG{1} = AnimationProgressBar_Update(...
        'Update', ANIMATION.handle,...
        'tag', 'Bar1',...
        'animation', ANIMATION.SEQ(ANIMATION.Counter),...SEQUENCE(FRAME),...
        'value', CLUSTER/NumClusters,...
        'message', '',...
        'process', 'Cluster Statistics',...
        'info', sprintf('Computing statistics for cluster %d', CLUSTER));


    % Calculate StDev Distribution:
    Edges = 0 : (1-0)/1000 : 1.0; % 0.6; %why 0.6? max St Dev!!

    ClusterStats(CLUSTER,2) = std(CLUSTERSTATS{CLUSTER});
    ClusterStats(CLUSTER,1) = mean(CLUSTERSTATS{CLUSTER});

    NumWindows = numel(CLUSTERSTATS{CLUSTER});
    
    
    % Send off to a Monte Carlo Simulator, which outputs a normalized
    % binomial random number/vector:
    StDevDistr = MonteCarloPrediction(...
        Edges,...
        NumWindows,... % number of windows with WINDOW photons
        WINDOW,... % number of photons within a window
        ClusterCenters(CLUSTER),... % probability of stoichiometry
        appData.BVA.Config.MonteCarlo.Prediction.NumSamples);
    CumDistr = cumsum(StDevDistr);
    INV = 1 - CumDistr;
    StDevDistrMean(CLUSTER) = sum(StDevDistr .* Edges); % What is this?

    %-- CI --
    %-- ADJUST CI BY THE  MULT HYPOTH TEST --
    absvec = abs(CumDistr - appData.BVA.Config.Dynamics.Criterion.ConfidenceInterval);
    [M, closest] = min(absvec);
%     closest = find(absvec == min(absvec),...
%         1, 'first');
%     closest = find(absvec==min(absvec));

    absvecLOW = abs(INV - appData.BVA.Config.Dynamics.Criterion.ConfidenceInterval);    
    closestLOW = find(absvecLOW == min(absvecLOW),...
        1, 'last');
%     closestLOW = find(absvecLOW==min(absvecLOW));

    % Why are we calculating the CI off an index? Index must tell us
    % something of about the number of false positives.
%     try
    HighCIvec(CLUSTER,:) = (closest-1)*.001; % 0.001 is 99.9% Confidence interval
    LowCIvec(CLUSTER,:) = (closestLOW-1)*.001;
%     catch
%         keyboard
%     end

        

end % end FOR


    

% ------------- %
%

appData.BVA.Config.ShotNoise.Prediction.StDev = sqrt(...
        (appData.BVA.Config.ShotNoise.Prediction.OrdinateVals .* (1-appData.BVA.Config.ShotNoise.Prediction.OrdinateVals))...
        /...
       WINDOW);

% Assign over to application data:
appData.BVA.Data.Clusters.Centers = ClusterCenters;
appData.BVA.Data.ConfidenceInterval.Upper = HighCIvec;
appData.BVA.Data.ConfidenceInterval.Lower = LowCIvec;
appData.BVA.Data.Clusters.StDevDistrMean = StDevDistrMean;
appData.BVA.Data.Clusters.Statistics = ClusterStats;
appData.BVA.Data.Clusters.WindowMeans = CLUSTERSTATS;
appData.BVA.Data.Clusters.DistanceFromUpperCI =...
    ClusterStats(:,2) - HighCIvec;
    

% ------------- %
%








% ------------------------------------------------------------------- %
%% Tidy-up:

%     clear ArrivalMatrix
%     clear BinIndices

    
% Closer animated waitbar:
delete(ANIMATION.handle);


TOTALTIME = toc



%% FIGURES:
tagLine = ' Acceptor Blinking Unfiltered ';
statePair = '(2E, 2S)';
states = '(S0.7,1.0, E0.5,0.0)';
rates = '(1000,1000)';
svgLine = 'acceptorblinking_unfiltered_k1000_k1000';


% Correlation versus Significance:
hFig = figure; 
plot(appData.BCA.Data.Results(:,5), appData.BCA.Data.Results(:,6), 'k.')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');
saveas(hFig,...
    ['ES Corr - Correlation versus PVals -', tagLine, statePair, states, rates, '(W', int2str(WINDOW), 'C',int2str(CORRWINDOW),')(2012-08-08)(1).fig'],...
    'fig');
saveas(hFig,...
    ['ES Corr - Correlation versus PVals -', tagLine, statePair, states, rates, '(W', int2str(WINDOW), 'C',int2str(CORRWINDOW),')(2012-08-08)(1).png'],...
    'png');
plot2svg(['corrPVals_', svgLine, '_W',int2str(WINDOW), 'C',int2str(CORRWINDOW),'.svg'], hFig);

delete(hFig);



% Hist of corr values
hFig = figure;
edges = -1.01 : 0.10 : 1.01; 
counts = histc(appData.BCA.Data.Results(:,5), edges);
bar(edges, counts, 'histc');
xlim([-1,1])
%title('');
xlabel('Correlation');
ylabel('Occurrence');

saveas(hFig,...
    ['ES Corr - Correlation Histogram -', tagLine, statePair, states, rates, '(W', int2str(WINDOW), 'C',int2str(CORRWINDOW),')(2012-08-08)(1).fig'],...
    'fig');
saveas(hFig,...
    ['ES Corr - Correlation Histogram -', tagLine, statePair, states, rates, '(W', int2str(WINDOW), 'C',int2str(CORRWINDOW),')(2012-08-08)(1).png'],...
    'png');
plot2svg(['corrHist_', svgLine, '_W',int2str(WINDOW), 'C',int2str(CORRWINDOW),'.svg'], hFig);

delete(hFig);



% Plot of corr values:
hFig = figure; 
plot(appData.BCA.Data.Results(:,5),'k.')
title('Burst Correlation');
xlabel('Burst');
ylabel('Correlation');
ylim([-1,1]);
xlim([0, size(appData.BCA.Data.Results(:,5),1)]);


saveas(hFig,...
    ['ES Corr - Correlation versus Bursts -', tagLine, statePair, states, rates, '(W', int2str(WINDOW), 'C',int2str(CORRWINDOW),')(2012-08-08)(1).fig'],...
    'fig');
saveas(hFig,...
    ['ES Corr - Correlation versus Bursts -', tagLine, statePair, states, rates, '(W', int2str(WINDOW), 'C',int2str(CORRWINDOW),')(2012-08-08)(1).png'],...
    'png');
plot2svg(['corrPlot_', svgLine, '_W',int2str(WINDOW), 'C',int2str(CORRWINDOW),'.svg'], hFig);

delete(hFig);



%% SOUND (indicate done)

% N=10000;
% s=zeros(N,1);
% for a=1:N
% s(a)=tan(a); %*sin(-a/10);
% end
% Fs=2000; %increase value to speed up the sound, decrease to slow it down
% soundsc(s,Fs)


bca = appData.BCA.Data.Results;
sbva = appData.BVA.Data.Results;
save(['Simulation - Photon Recoloring - ES Corr - SBVA - ', tagLine, statePair, states, rates, '(W', int2str(WINDOW), 'C',int2str(CORRWINDOW),')(2012-08-08)(1).mat'],...
    'bca',...
    'eStateMeans',...
    'sStateMeans',...
    'eInitialProbabilities',...
    'sInitialProbabilities',...
    'eRateMatrix',...
    'sRateMatrix',...
    'mixProbs',...
    'gamma',...
    'WINDOW',...
    'CORRWINDOW',...
    'sbva');


%keyboard



return;

myIDS = find(~isnan(appData.BCA.Data.Results(:,6)) & appData.BCA.Data.Results(:,6) < 0.05);
bca = appData.BCA.Data.Results;
save('CorrelatedESbursts.mat', 'myIDS', 'bca');
REMOVEBURSTS = unique([REMOVEBURSTS; myIDS]);


% Correlation versus pVals
figure; 
plot(appData.BCA.Data.Results(:,5), appData.BCA.Data.Results(:,6), 'k.')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');

% Num windows:
figure; 
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,7) <= 25 & appData.BCA.Data.Results(:,7) >= 10, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,7) <= 25 & appData.BCA.Data.Results(:,7) >= 10, 6), 'b*',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,7) > 25, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,7) > 25, 6), 'ro',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,7) < 10, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,7) < 10, 6), 'go')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');

% Num corr windows
figure; 
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) <= 10 & appData.BCA.Data.Results(:,8) >= 5, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) <= 10 & appData.BCA.Data.Results(:,8) >= 5, 6), 'b*',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) > 10, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) > 10, 6), 'ro',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) < 5, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) < 5, 6), 'go')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');


% Num corr windows
figure; 
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) <= 10 & appData.BCA.Data.Results(:,8) >= 5, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) <= 10 & appData.BCA.Data.Results(:,8) >= 5, 9), 'b*',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) > 10, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) > 10, 9), 'ro',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) < 5, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) < 5, 9), 'go')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');

% Num corr windows
figure; 
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) >= 5, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) >= 5, 9), 'ro')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');

% Num states plot
figure; 
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,10) >= 2, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,10) >= 2, 9), 'ro',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,10) < 2, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,10) < 2, 9), 'b*')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');

% Num states plot
figure; 
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,10) >= 2, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,10) >= 2, 6), 'ro',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,10) < 2, 5),...
        appData.BCA.Data.Results(appData.BCA.Data.Results(:,10) < 2, 6), 'b*')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');


% Num transitions plot
figure; 
IDs1 = (appData.BCA.Data.Results(:,12) >= 10);
IDs2 = (appData.BCA.Data.Results(:,12) < 10);
plot(...
    appData.BCA.Data.Results(IDs1, 5),...
        appData.BCA.Data.Results(IDs1, 9), 'ro',...
    appData.BCA.Data.Results(IDs2, 5),...
        appData.BCA.Data.Results(IDs2, 9), 'b*')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');

% Num transitions plot
figure; 
IDs1 = (appData.BCA.Data.Results(:,12) >= 5);
IDs2 = (appData.BCA.Data.Results(:,12) < 5);
plot(...
    appData.BCA.Data.Results(IDs1, 5),...
        appData.BCA.Data.Results(IDs1, 6), 'ro',...
    appData.BCA.Data.Results(IDs2, 5),...
        appData.BCA.Data.Results(IDs2, 6), 'b*')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');


% Num hot spots plot
figure; 
IDs1 = (appData.BCA.Data.Results(:,15) >= 5);
IDs2 = (appData.BCA.Data.Results(:,15) < 5);
plot(...
    appData.BCA.Data.Results(IDs2, 5),...
        appData.BCA.Data.Results(IDs2, 9), 'b*',...
        appData.BCA.Data.Results(IDs1, 5),...
        appData.BCA.Data.Results(IDs1, 9), 'ro')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');


% Num hot spots and num transitions plot
figure; 
IDs1 = (appData.BCA.Data.Results(:,12) < 5);
IDs2 = (appData.BCA.Data.Results(:,15) < 3);
plot(...
    appData.BCA.Data.Results(IDs2, 5),...
        appData.BCA.Data.Results(IDs2, 9), 'ks',...
        appData.BCA.Data.Results(IDs1, 5),...
        appData.BCA.Data.Results(IDs1, 9), 'bo')
hold on
IDs1 = (appData.BCA.Data.Results(:,12) >= 5);
IDs2 = (appData.BCA.Data.Results(:,15) >= 3);
plot(...
    appData.BCA.Data.Results(IDs2, 5),...
        appData.BCA.Data.Results(IDs2, 9), 'g*',...
        appData.BCA.Data.Results(IDs1, 5),...
        appData.BCA.Data.Results(IDs1, 9), 'ro')
hold off
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');

% Hist of hot spot numbers
figure;
edges = 0 : 1 : 30; 
counts = histc(appData.BCA.Data.Results(:,15), edges);
bar(edges, counts, 'histc');
xlim([0,30])


% Hist of corr values
figure;
edges = -1.01 : 0.02 : 1.01; 
counts = histc(appData.BCA.Data.Results(appData.BCA.Data.Results(:,8) >= 5,5), edges);
bar(edges, counts, 'histc');
xlim([-1,1])


% Hist of corr values
figure;
edges = -1.01 : 0.10 : 1.01; 
counts = histc(appData.BCA.Data.Results(:,5), edges);
bar(edges, counts, 'histc');
xlim([-1,1])
%title('');
xlabel('Correlation');
ylabel('Occurrence');

% Plot of corr values:
figure; 
plot(appData.BCA.Data.Results(:,5),'k.')
title('Burst Correlation');
xlabel('Burst');
ylabel('Correlation');
ylim([-1,1]);
xlim([0, size(appData.BCA.Data.Results(:,5),1)]);


% Plot of transitions versus hot spots:
figure; 
plot(appData.BCA.Data.Results(:,12),...
    appData.BCA.Data.Results(:,15),'b*')

% Plot of correlation versus hot spots:
figure; 
plot(appData.BCA.Data.Results(:,5),...
    appData.BCA.Data.Results(:,15),'b*')





figure; 
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,7) > 14 & appData.BCA.Data.Results(:,7) < 31, 5),...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,7) > 14 & appData.BCA.Data.Results(:,7) < 31, 6), 'b*')
title('Correlation and Significance');
xlabel('Correlation');
ylabel('P Value');


figure; 
plot(appData.BCA.Data.Results(:,6), appData.BCA.Data.Results(:,7), 'b*')
title('Significance and NumWindows');
xlabel('P Value');
ylabel('NumWindows');


figure; 
plot(appData.BCA.Data.Results(:,5), appData.BCA.Data.Results(:,7), 'b*',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.15, 5), appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.15,7), 'g.',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.05, 5), appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.05,7), 'ro')
title('Correlation and NumWindows');
xlabel('Correlation');
ylabel('NumWindows');


% Isolating bursts based on correlation: are highly correlated E/S bursts
% more dynamic than uncorrelated?
figure; 
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) > 0.1,1),...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) > 0.1,2), 'k.',...
    'markersize', 15)
E = 0:0.01:1; SN = sqrt(E.*(1-E)/7); hold on; plot(E, SN, 'r'); hold off;
title('E-BCA: uncorrelated');


figure; 
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.05,1),...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.05,2), 'k.',...
    'markersize', 15)
E = 0:0.01:1; SN = sqrt(E.*(1-E)/7); hold on; plot(E, SN, 'r'); hold off;
title('E-BCA: correlated');

figure; 
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.005,1),...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.005,2), 'k.',...
    'markersize', 15)
E = 0:0.01:1; SN = sqrt(E.*(1-E)/7); hold on; plot(E, SN, 'r'); hold off;
title('E-BCA: highly correlated');


figure; 
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) > 0.1,3),...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) > 0.1,4), 'k.',...
    'markersize', 15)
S = 0:0.01:1; SN = sqrt(S.*(1-S)/10); hold on; plot(S, SN, 'r'); hold off;
title('S-BCA: uncorrelated');

figure;
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.05,3),...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.05,4), 'k.',...
    'markersize', 15)
S = 0:0.01:1; SN = sqrt(S.*(1-S)/10); hold on; plot(S, SN, 'r'); hold off;
title('S-BCA: correlated');

figure;
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.005,3),...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.005,4), 'k.',...
    'markersize', 15)
S = 0:0.01:1; SN = sqrt(S.*(1-S)/10); hold on; plot(S, SN, 'r'); hold off;
title('S-BCA: highly correlated');




figure;
plot(...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) > 0.1,3),...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) > 0.1,4), 'b.',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.05,3),...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.05,4), 'k.',...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.005,3),...
    appData.BCA.Data.Results(appData.BCA.Data.Results(:,6) < 0.005,4), 'r.',...
    'markersize', 15)
S = 0:0.01:1; SN = sqrt(S.*(1-S)/10); hold on; plot(S, SN, 'r'); hold off;
title('S-BCA: Species');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% CORRELATION COEFFICIENT:
function corrCoefficient = correlationCoefficient(X, Y)
    %
    %

    % NOTE: the calculated mean is the sample mean, and hence
    % local. Results will differ from the case where the population
    % mean is used. 

    % Calculate the means:
    muX = nanmean(X);
    muY = nanmean(Y);

    % Determine the fluctuation amplitude:
    fluctX = X - muX;
    fluctY = Y - muY;

    % Calculate the variances:
    varX = sum(fluctX.^2);
    varY = sum(fluctY.^2);

    % Calculate the product of X,Y fluctuations:
    fluctXY = fluctX .* fluctY;

    % Calculate the cross-correlation:
    corrCoefficient = sum(fluctXY) / (sqrt(varX*varY) + eps); % denominator includes a very small number to prevent division by zero when Var = 0;
    
    
    
    
%% continuousTimeKineticMonteCarlo:
function markovChain = continuousTimeKineticMonteCarlo(fid, totalTime,...
        initialProbabilities, rateMatrix)
    %
    %

    % Determine the number of states:
    numStates = size(initialProbabilities(:),1);

    % Initialize a vector of ones:
    onesVec = ones(numStates, 1);

    % Let the games begin...


    % ********* %
    % [1] Find out the initial state:

    % Draw a random number:
    randValVec = rand * onesVec; % Faster not to provide arguments to 'rand'

    % Generate a cumulative histogram from the initial probabilities:
    cumProbHist = cumsum(initialProbabilities);

    % Find where our random value resides among the classes:
    classy = cumProbHist - randValVec;

    % Determine to what class to assign the random value, updating our current
    % state:
    currentState = find(classy > 0, 1, 'first');

    % Place the current state in the Markov Chain array:
    markovChain(1,1) = currentState;

    

    % ********* %
    % [2] Find the initial transition time:

    % Get the kinetic rates for the initial state:
    rates = rateMatrix(currentState, :);

    % Remove the current state entry (which has a negative value):
    rates(currentState) = 0;

    % Draw a uniform random number:
    randVal = rand;

    % Generate an exponentially distributed random number using
    % inverse transform sampling (see Wikipedia)
    waitingTime = -1 * log(randVal) / sum(rates); % Dwell drawn from an exponential distribution

    % Place the absolute time in the Markov Chain array:
    markovChain(1,2) = 0+waitingTime; % t = 0
    
    % Place the waiting time in the Markov Chain array:
    markovChain(1,3) = waitingTime;


%     % Write the state to file:
%     fprintf(fid, '%d \t %10.12e \t %10.12e \n', markovChain(1,1), markovChain(1,2), markovChain(1,3));



    % ********* %
    % [3] Simulation for the remaining time:            
    time = waitingTime;
    while time <= (totalTime + 1e-10) % Random error: due to machine precision, sometimes the TimeClock does not equal TotalTime but is very slightly greater than (1e-15). So, we need to add a very small number to ensure that our last INDEX is filled in our vectors.

        % Draw a random number:
        randVal = rand; 

        % Generate a cumulative histogram:
        rates = rateMatrix(currentState, :);

        % Remove the current state, as we are now transitioning away
        % now that the dwell time period has been eclipsed...
        rates(currentState) = 0;

        % Generate a cumulative sum of the rates:
        cumRateHist = cumsum(rates);

        % Multiply the end of our cumulative sum by the drawn random
        % number:
        randVal = cumRateHist(end) * randVal;

        % Find where our random value resides among the classes:
        classy = cumRateHist - randVal;

        % Determine to what class to assign the random value, updating
        % our current state:
        currentState = find(classy > 0, 1, 'first');

        % Get our new rates:
        rates = rateMatrix(currentState, :);
        rates(currentState) = 0;

        % Draw a 'dwell' from an exponential distribution:
        randVal = rand;
        waitingTime = -1 * log(randVal) / sum(rates);

        % Update our Markov Chain:
        markovChain(end+1,:) = [currentState, time+waitingTime, waitingTime];


%         % Write the state to file:
%         fprintf(fid, '%d \t %10.12e \t %10.12e \n', markovChain(1,1), markovChain(1,2), markovChain(1,3));


        % Update our simulation time:
        time = time + waitingTime;

    end % end WHILE

            
        
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% StartUp:
function StartUp()
%
%
%

S = what('Seneca');

% Add path(s):
addpath(...
    [S.path,...
    '/Lightspeed'],...
    [S.path,...
    '/MatlabCentralFileExchange',...
    '/SplitVector'],...
    [S.path,...
    '/SMFiles',...
    '/Conversion'],...
    [S.path,...
    '/MatlabCentralFileExchange',...
    '/MultidimensionalHistograms',...
    '/nDimensionalHistogram'],...
    [S.path,...
    '/MatlabCentralFileExchange',...
    '/SmoothData'],...
    [S.path,...
    '/MatlabCentralFileExchange',...
    '/VisualTouchUps',...
    '/myaa'],...
    [S.path,...
    '/MatlabCentralFileExchange',...
    '/WaitbarAlternatives',...
    '/ImageAndProgressBar']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CleanUp:
function CleanUp()
%
%
%

S = what('Seneca');

% Remove path(s):
rmpath(...
    [S.path,...
    '/Lightspeed'],...
    [S.path,...
    '/MatlabCentralFileExchange',...
    '/SplitVector'],...
    [S.path,...
    '/SMFiles',...
    '/Conversion'],...
    [S.path,...
    '/MatlabCentralFileExchange',...
    '/MultidimensionalHistograms',...
    '/nDimensionalHistogram'],...
    [S.path,...
    '/MatlabCentralFileExchange',...
    '/SmoothData'],...
    [S.path,...
    '/MatlabCentralFileExchange',...
    '/VisualTouchUps',...
    '/myaa'],...
    [S.path,...
    '/MatlabCentralFileExchange',...
    '/WaitbarAlternatives',...
    '/ImageAndProgressBar']);


