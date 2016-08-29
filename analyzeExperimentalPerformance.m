% Copyright (C) 2016  Gregory L. Futia
% This work is licensed under a Creative Commons Attribution 4.0 International License.

% Description: Iterates throug the days of data. Data is split into 
% training-testing subsets for each days worth of data. 
% For each training testing subset, ROC curves of the individual features
% are created. Regressions are performed using the training side of each
% subset and operating point maximizing Ndet is found. Regressiosn are
% tested with the testing subset ath the previously found operationg point.
% Performance statistics are found for each of the Training Testing
% Subsets. 
% Function depends on :
% selectAndRegress2.m
% determineFeatureStat.m
% determineRegressionStat.m
% create_manuscript_histograms.m
% parsave.m

mainPath = [pwd, filesep, 'Processed Mats', filesep, 'Futia_Zenodo_2016'];

daysToProcess = {'day5', 'day7', 'day8', 'day9', 'day10', 'day12', 'day13', ...
    'day14', 'day15'};

daysToProcess = {'day5', 'day7',};

daysToProcess = {'day9', 'day10', 'day12', 'day13', ...
    'day14', 'day15'};


% this function also does the loading and merging

% for each day there are two folders. One of pure MCF7 data and one of pure
% WBC data. For each pairing of this data, compute cohen's D, and and ROC
% curve. Find the positon on the ROC curve that maximizes Ndet
% for each of those cut off positions ask what fp and sens is in remaining
% testing data. 
fsepchar = filesep;

DNACutoff = 65;

% load the regressions form files
% otherwise compute them
bLoad = 0; % 0 for compute; 1 for load. Running code in compute can take days.
bUsePooledRegs = 0; % the all regression instead of the training-testing pairings
bSwapTrainTest = 1; % 0 small train big test (normal), 1 big train small test 

Ndays = length(daysToProcess);


WBCDataAll = cell(1, Ndays);
MCF7DataAll = cell(1, Ndays);
MixDataAll = cell(1, Ndays);

% This will interate through each day's worth of data
for ii=1:length(daysToProcess)
    fprintf('Processing %s\n', daysToProcess{ii});    
    % Figure out the numbe of WBC and MCF7 file in these days
    
    wbcPath = [mainPath, fsepchar , daysToProcess{ii}, fsepchar, 'WBC', ...
        fsepchar];
    dirResult = dir([wbcPath '*.mat']);    
    
    Nwbcfiles = length(dirResult);
    wbcData = cell(1, Nwbcfiles);
    wbcFilenames = cell(1, Nwbcfiles);
    
    for jj=1:Nwbcfiles
        wbcData{jj} = load([wbcPath dirResult(jj).name]);
        wbcFilenames{jj} = dirResult(jj).name;
    end    
    
    ND = combineND(wbcData);
    % save the combined WBC data in the root day path with the day after
    % the file name
    parsave([mainPath, fsepchar , daysToProcess{ii}, fsepchar, 'wbc_', daysToProcess{ii}, '.mat'], ND);
    WBCDataAll{ii}.ND = ND;
    
    mcf7Path = [mainPath, fsepchar , daysToProcess{ii}, fsepchar, 'MCF7', ...
        fsepchar];
    dirResult = dir([mcf7Path, '*.mat']);       
    Nmcf7files = length(dirResult);
    mcf7Filenames = cell(1, Nmcf7files);
    
    mcf7Data = cell(1, Nmcf7files); 
   
    for jj=1:Nmcf7files
        mcf7Data{jj} = load([mcf7Path dirResult(jj).name]);
        mcf7Filenames{jj} = dirResult(jj).name;
    end
    
    ND = combineND(mcf7Data);
    
    % save the combined MCF7 data in the root day path with the after the
    % file name
    parsave([mainPath, fsepchar , daysToProcess{ii}, fsepchar, 'mcf7_', daysToProcess{ii}, '.mat'], ND);
    MCF7DataAll{ii}.ND = ND;
   
    mixPath = [mainPath, fsepchar , daysToProcess{ii}, fsepchar, 'Mix', ...
        fsepchar];
    dirResult = dir([mixPath, '*.mat']);
    Nmixfiles = length(dirResult);
    mixFilenames = cell(1, Nmixfiles);
    
    mixData = cell(1, Nmixfiles);
    
    for jj=1:Nmixfiles
        mixData{jj} = load([mixPath dirResult(jj).name]);
        mixFilenames{jj} = dirResult(jj).name;
    end    
    % group all the mixed data togeteher for one test set
    ND = combineND(mixData);
    mixTest = ND;
    
    % save the combined MCF7 data in the root day path with the after the
    % file name
    parsave([mainPath, fsepchar , daysToProcess{ii}, fsepchar, 'mix_', daysToProcess{ii}, '.mat'], ND);
    MixDataAll{ii}.ND = ND;
    
    pathToRegMats = [mainPath, fsepchar, daysToProcess{ii}, fsepchar, 'Reg Mats' fsepchar];
    pathToEqns = [mainPath, fsepchar, daysToProcess{ii}, fsepchar, 'Reg Equations' fsepchar];
    % if these directories don't exist create them
    if ~exist(pathToRegMats, 'dir')
        mkdir(pathToRegMats);
    end
    
    if ~exist(pathToEqns, 'dir')
        mkdir(pathToEqns);
    end
        
    featStats = cell(1, Nwbcfiles*Nmcf7files);
    regStats = cell(1, Nwbcfiles*Nmcf7files);
    indexFeatStats = 1;
    % iterate through each of the file name combinations
    for jj = 1:Nwbcfiles
        for kk = 1:Nmcf7files
            % this will be the training data
            wbcTrain = wbcData{jj};
            mcf7Train = mcf7Data{kk};
            wbcTrain = wbcTrain.ND;
            mcf7Train = mcf7Train.ND;
            
            % the testing data                    
            wbcTest = combineND(wbcData(1:Nwbcfiles~=jj));            
            mcf7Test = combineND(mcf7Data(1:Nmcf7files~=kk));
            
            if bSwapTrainTest == 1
                wbcTemp = wbcTrain;
                mcf7Temp = mcf7Train;
                
                wbcTrain = wbcTest;
                mcf7Train = mcf7Test;
                wbcTest = wbcTemp;
                mcf7Test = mcf7Temp;                
            end
                
            % featsStats is array of performance stats for the features
            % determineFeatureStats computes these metrics
            featStats{indexFeatStats} = determineFeatureStat(mcf7Train, wbcTrain, ...
                mcf7Test, wbcTest, DNACutoff);
           
            regMatsPath = pathToRegMats;
            eqnPath = pathToEqns;   
            
            if bUsePooledRegs == 0
               % compute the regressions on the training-testing pairs
               % in reg mats we need to create a folder called regvars
               % take the .mat off the wbc & mcf7 file names
               % this will be the name for this folder pairings
               regDirName = [wbcFilenames{jj}(1:end-4) '_' mcf7Filenames{kk}(1:end-4)];
               regMatsPath = [pathToRegMats, regDirName, filesep];
               eqnPath = [pathToEqns, regDirName, filesep];
               
               if ~exist(regMatsPath', 'dir')
                   mkdir(regMatsPath);
               end
                              
               if ~exist(eqnPath', 'dir')
                   mkdir(eqnPath);
               end
               
               % perform features selection and regression on the training
               % data
               if bLoad == 0
                % this will populate the RegMats directory if it isn't
                % already populated. 
                selectAndRegress2(wbcTrain, mcf7Train, regMatsPath, eqnPath);               
               end
            end
            
            % overloading in Matlab isn't so easy, changed to bring ND in
            % as 1x4 or 1x5 cell array so mixed data could be added
            % to the regStats
            regStats{indexFeatStats} = determineRegressionStat({mcf7Train, wbcTrain, ...
                mcf7Test, wbcTest, mixTest}, DNACutoff, regMatsPath);
            
%             featStatsAll{Npairmax*(ii-1) + Nwbcfiles*(jj-1) + kk} = featStats{indexFeatStats};
%             regStatsAll{Npairmax*(ii-1) + Nwbcfiles*(jj-1) + kk} = regStats{indexFeatStats};
%             featStatsAllName{Npairmax*(ii-1) + Nwbcfiles*(jj-1) + kk} = daysToProcess{ii};
            indexFeatStats = indexFeatStats + 1;
            %featStatsAllIndex = featStatsAllIndex + 1;
        end                
    end
    parsave([mainPath, fsepchar, daysToProcess{ii}, fsepchar, 'featStats.mat'], featStats);
    parsave([mainPath, fsepchar, daysToProcess{ii}, fsepchar, 'regStats.mat'], regStats);
        
    featFullData = createStatsTabAndPlots2(featStats, [mainPath, fsepchar, daysToProcess{ii}, fsepchar], 'feature');    
    createStatsTabAndPlots2(regStats, [mainPath, fsepchar, daysToProcess{ii}, fsepchar], 'regression');
    figure(6);
    clf;
    create_manuscript_histograms(MCF7DataAll{ii}.ND, WBCDataAll{ii}.ND, MixDataAll{ii}.ND, featFullData, DNACutoff, ...
        [mainPath, fsepchar, daysToProcess{ii}, fsepchar]);
end

% combine all the data and save the mats
ND = combineND(WBCDataAll);
save([mainPath, fsepchar , 'wbcAll', '.mat'], 'ND');
Dn = ND; % use for histograms

% combine all the data and save the mats
ND = combineND(MCF7DataAll);
save([mainPath, fsepchar , 'mcf7All', '.mat'], 'ND');
Dp = ND; % use for histograms

% combine all the data and save the mats
ND = combineND(MixDataAll);
save([mainPath, fsepchar , 'mixAll', '.mat'], 'ND');
Dmix = ND; % use for histograms


% this used to be in the above loop but is incompatible with parfor
Npairmax = 9;
% allocate an array for all the results
% 9 is more than we are expecting some days will have less than 9 pairings
featStatsAll = cell(1, Npairmax*length(daysToProcess));
featStatsAllName = cell(1, Npairmax*length(daysToProcess)); % this is 
% will be populated with the day of the featStats

% can reuse the feature name field and featureStatAllIndex
regStatsAll = cell(1, Npairmax*length(daysToProcess));

featStatsAllIndex = 1; % we will use the same index for both of these
for ii=1:length(daysToProcess)
    load([mainPath, fsepchar, daysToProcess{ii}, fsepchar, 'featStats.mat'], 'featStats');
    load([mainPath, fsepchar, daysToProcess{ii}, fsepchar, 'regStats.mat'], 'regStats');
    Nload = length(featStats);
    
    featStatsAll(featStatsAllIndex:featStatsAllIndex+Nload-1) = featStats;
    regStatsAll(featStatsAllIndex:featStatsAllIndex+Nload-1) = regStats;
    featStatsAllName(featStatsAllIndex:featStatsAllIndex+Nload-1) = repmat(daysToProcess(ii), 1, Nload);
    featStatsAllIndex = featStatsAllIndex + Nload;    
end


% remove empty elements
featStatsAll = featStatsAll(~cellfun('isempty', featStatsAll));
featStatsAllName = featStatsAllName(~cellfun('isempty',featStatsAllName));
regStatsAll = regStatsAll(~cellfun('isempty',regStatsAll));

% save the performance mats
save([mainPath, fsepchar, 'featStatsAll.mat'], 'featStatsAll', 'featStatsAllName');
save([mainPath, fsepchar, 'regStatsAll.mat'], 'regStatsAll', 'featStatsAllName');
featFullData = createStatsTabAndPlots2(featStatsAll, [mainPath, fsepchar], 'feature', featStatsAllName);
regFullData = createStatsTabAndPlots2(regStatsAll, [mainPath, fsepchar], 'regression', featStatsAllName);

figure(6);  
clf;
create_manuscript_histograms(Dp, Dn, Dmix, featFullData, DNACutoff, [mainPath, fsepchar]);


