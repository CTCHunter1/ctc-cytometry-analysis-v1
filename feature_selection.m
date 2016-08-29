% feature_selection.m

% performs feature selection to produce maximum descrimination between
% two classes. Uses discriminant analysis. 

% Copyright (C) 2016  Gregory L. Futia
% This work is licensed under a Creative Commons Attribution 4.0 International License.


function [ ] = feature_selection(selectedChannels, selectedFeatures, regName, filenames, saveFilePath)

if(~exist('filenames', 'var'))
    if ~exist('filenames.mat', 'file')
        [filenames{1}, pathnames{1}] = uigetfile('*.mat', 'Select WBC File', 'C:\Repository\Projects\Nuclear Second Moment\Code\Processed Mats');    

        if isequal(filenames{1},0) || isequal(pathnames{1},0)
            return;
        end

        [filenames{2}, pathnames{2}] = uigetfile('*.mat', 'Select Cancer Cell File', 'C:\Repository\Projects\Nuclear Second Moment\Code\Processed Mats');    

        if isequal(filenames{2},0) || isequal(pathnames{2},0)
            return;
        end
    end
else
    pathnames = cell(1, length(filenames));
    for ii = 1:2
        [pathnames{ii}, filenames{ii}, ext] = fileparts(filenames{ii});
        filenames{ii} = [filenames{ii}, ext];
    end
end

if(~exist('saveFilePath', 'var'))
    [saveFilePath, ~, ~] = fileparts(filenames{1});
end

% put a / at the end of save file path if it doesn't already have one
if saveFilePath(end) ~= filesep
    saveFilePath = [saveFilePath, filesep];
end
    
% For discriminant analysis we will attempt to regress against WBC=0,
% Cancer = 1
NumFiles = 2; % the last cell is empty


% determine number of points
Npts = 0;
data = cell(1, NumFiles);

% determine how many data points there are
for ii = 1:NumFiles;    
    data{ii} = load(fullfile(pathnames{ii}, filenames{ii}));
    Npts = Npts + data{ii}.ND.numROIs;
end

ND = combineND(data);
ND.Y = zeros(1, ND.numROIs, 'uint8');
% the second file is the MCF7 file set all those points equal to 1
ND.Y(data{1}.ND.numROIs + 1:end) = 1; 

% determine the number of features = num chnannels x metrics
numFeatures = 0;
channelNames = fieldnames(ND.Channels);
NchannelNames = length(channelNames);

if exist('selectedChannels', 'var') == 0 
    selectedChannels = {'DAPI', 'PanCK', 'Bodipy', 'CD45'};
end

if exist('selectedFeatures', 'var') == 0 
    selectedFeatures = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
end

features1 = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
varEqn = regexprep(selectedFeatures, features1, {'\Sigma', '<r>', '<r_f>', '<M>'});
unitsEqn = regexprep(selectedFeatures, features1, {' dBct ', ' \mu m ', ' \mu m^{-1} ', ''});

% selectedFeatures = {'totalSig_dBc'};
% varEqn = {'\Sigma'};
% unitsEqn = {' dBct '};

for ii = 1:NchannelNames
    featureNames = fieldnames(ND.Channels.(channelNames{ii}));
    % select the channels
    if sum(strcmp(channelNames{ii}, selectedChannels))
        % select the features
        for jj = 1:length(featureNames)
            if sum(strcmp(featureNames{jj}, selectedFeatures))
                [M, N] = size(ND.Channels.(channelNames{ii}).(featureNames{jj}));
                if (N > 1)
                    numFeatures = numFeatures + M; % compensates for 2D arrays
                end
            end
        end
    end
end

% set the DNA gate
dnaCutoff = 65;
gate = ND.Channels.DAPI.totalSig_dBc > dnaCutoff;


allData = zeros(numFeatures, ND.numROIs);
features = cell(1, numFeatures);
varnames = cell(1, numFeatures);
varnamesTable = cell(1, numFeatures);
units = cell(1, numFeatures);
d = zeros(1, numFeatures);
% now populate allData
featureindex = 1;
for ii = 1:NchannelNames
    featureNames = fieldnames(ND.Channels.(channelNames{ii}));
    % select the channels
    if sum(strcmp(channelNames{ii}, selectedChannels))
        for jj = 1:length(featureNames)
            [M, N] = size(ND.Channels.(channelNames{ii}).(featureNames{jj}));
            if (N > 1)
                % select the features
                % check the names matches the type we are looking for                
                varType = strcmp(featureNames{jj}, selectedFeatures);
                if sum(varType)
                    % features are selected here
                    allData(featureindex:featureindex+(M-1), :) = ND.Channels.(channelNames{ii}).(featureNames{jj});                
                   d(featureindex) = ... 
                        cohensD(ND.Channels.(channelNames{ii}).(featureNames{jj})(logical(ND.Y) & gate), ...
                        ND.Channels.(channelNames{ii}).(featureNames{jj})(~logical(ND.Y) & gate));


                    if M == 1
                    features{featureindex} = [channelNames{ii}, '.' , featureNames{jj}];
                    % build out the variable names list
                    varnames{featureindex} = [varEqn{varType}, '_{', channelNames{ii}, '}'];
                    % tex formated table var names
                    varnamesTable{featureindex} = [channelNames{ii}, ' $', varEqn{varType}, '$ ' ];

                    units{featureindex} = unitsEqn{varType};
                    else
                        for kk = 1:M
                        features{featureIndex + kk -1} = [channelNames{ii}, '.' , featureNames{jj}, 'OD', num2str(kk*2)];
                        end                
                    end
                    featureindex = featureindex + M;
                end
            end
        end
    end
end


Y = ND.Y(gate);

allData = allData(:, gate); % gate the data
wbcData = allData(:, Y==0);
mcf7Data = allData(:, Y==1);

% iterate through all combinations of all variables and perform feature selection
maxD = 0;
maxD_ind = 0;
maxD_numVar = 0;
Nfold = 10;
cvp = cvpartition(Y,'kfold', Nfold);
Y = Y;
% dReg holds the results from all the combinations of interations
% it is a cell array where each cell contains a jagged array resulting
% from the performance results from interations of all combations of that 
% many features
dReg = cell(1, numFeatures);
combinations = cell(1, numFeatures);
% we are going to keep the regression results from all of the final
% regressions
betaAll = cell(1, numFeatures);
muAll = cell(1, numFeatures);
sigAll = cell(1, numFeatures);
featuresAll = cell(1, numFeatures);
varnamesAll = cell(1, numFeatures);
unitsAll = cell(1, numFeatures);

regNameFile = strrep(regName, '\', ''); % remove slashes from save file 
% string
SNR = 2;

for ii = 1:numFeatures
    fprintf('%d feature regresssions \n', ii);
    C = nchoosek(1:numFeatures, ii);
    combinations{ii} = C;
    [Mc, Nc] = size(C);
    dReg{ii} = zeros(1,Mc);        
    betaAll{ii} = cell(1, Mc);
    muAll{ii} = cell(1, Mc);
    sigAll{ii} = cell(1, Mc);
    featuresAll{ii} = cell(1, Mc);
    varnamesAll{ii} = cell(1, Mc);
    unitsAll{ii} = cell(1, Mc);        

    for jj = 1:Mc
        dSum_Reg = 0;
        % select the data for the regression,
        % this selects the variables from the data in the combination
        % C(jj,:) 
        data = allData(C(jj,:),:).';
        fprintf('Regressing Feature Combination: %s \n', sprintf('%d ', C(jj,:)));

        for kk=1:Nfold
            %fprintf('On fold %d \n', kk);
            data_train = data(cvp.training(kk), :);
            [z, mu, sig] = zscore(data_train);                
            [len, ~] = size(z);
            beta_reg = regress(double(Y(cvp.training(kk))).',[ones(len, 1), z]);

            data_test = data(cvp.test(kk),:);
            Ytest = Y(cvp.test(kk));
            [Mt,Nt] = size(data_test);
            z_test = (data_test - ones(Mt,1)*mu)./(ones(Mt,1)*sig);
            scores_reg_test = [ones(Mt, 1), z_test]*beta_reg;

            % instead of optimizing on Cohen's D optimize on minimum
            % detectable threadholds
            dSum_Reg = dSum_Reg + abs(cohensD(scores_reg_test(logical(Ytest)), scores_reg_test(~logical(Ytest))));

            % this is slower than using cohen's d
%             [fp, sens, T, AUC] = compute_ROC(scores_reg_test(logical(Ytest)).', ...
%                 scores_reg_test(~logical(Ytest)).');
% 
%             % 0 false postive rate means we just didn't sample it
%             % set this to the sample size instead
%             fp(fp==0) = 1./length(scores_reg_test(~logical(Ytest)));
% 
%             [minDetect, minDetInd] = max((sens+SNR*fp)./(SNR*fp));
%             %FeatStat.fpMinDetect(ii) = fp(minDetInd);
%             %FeatStat.sensMinDetect(ii) = sens(minDetInd);
%             dSum_Reg = dSum_Reg + minDetect; % this is a hack for now
        end                     

        dReg{ii}(jj) =  dSum_Reg / Nfold;                     

        [d_max, ind] = max(abs(dReg{ii}));
        if(d_max > maxD)
             maxD = d_max;
             maxDInd = ind;
             maxDNumVar = ii;
        end                     

        varIndex = C(jj,:);
        featuresAll{ii}{jj} = features(varIndex);
        varnamesAll{ii}{jj} = varnames(varIndex);
        unitsAll{ii}{jj} = units(varIndex);
        [z_REG, muAll{ii}{jj}, sigAll{ii}{jj}] = zscore(allData(varIndex, :).');
        [len, ~] = size(z_REG);
        betaAll{ii}{jj} = regress(double(Y).', [ones(len, 1), z_REG]);
    end
end

varIndex = combinations{maxDNumVar}(maxDInd, :);
[z_REG, mu_REG, sigma_REG] = zscore(allData(varIndex, :).');
[len, ~] = size(z_REG);
betareg = regress(double(Y).', [ones(len, 1), z_REG]);
% features determines the order and structure of the matrix to build
if(~exist(saveFilePath, 'dir'))
    mkdir(saveFilePath);
end

save([saveFilePath, regNameFile, '.mat'], 'features', 'varIndex', 'z_REG', 'mu_REG', 'sigma_REG', 'betareg');
save([saveFilePath, regNameFile, '_all.mat'], 'featuresAll', 'muAll', 'sigAll', 'betaAll', 'dReg', 'varnamesAll', 'unitsAll');




