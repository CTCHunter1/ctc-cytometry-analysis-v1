function [FeatStat] = determineFeatureStat(Dp, Dn, Dptest, Dntest, DNACutOff)
% Copyright (C) 2016  Gregory L. Futia
% This work is licensed under a Creative Commons Attribution 4.0 International License.

% Description: Iterates throug the days of data. Data is split into 

if nargin == 0
    % this script is ment to be called by create_manuscript_figures
    % if this script is ran without parameters run that script instead
    estimateTechnicalVariances; 
    return;
end

% Set the DNA Gates
DNAp = Dn.Channels.DAPI.totalSig_dBc > DNACutOff;
DNApMCF7 = Dp.Channels.DAPI.totalSig_dBc > DNACutOff;
DNApDptest = Dptest.Channels.DAPI.totalSig_dBc > DNACutOff;
DNApDntest = Dntest.Channels.DAPI.totalSig_dBc > DNACutOff;

channelNames = fieldnames(Dp.Channels);
reorder = [3, 1, 4, 2]; % reorder the channel names

features = {'totalSig_dBc', 'radius_m', 'radius_invm', 'M2'};


% initialize the performance variables
% they aren't 'regressions' but feature nemes here
% this is for compatiblity with determineRegStat
regressionNames = cell(1, length(channelNames)*length(features));
Nregressions = length(regressionNames); 
FeatStat.regressionNames = regressionNames;
FeatStat.TMinDetect = zeros(1, Nregressions);
% these will be set using the test data
FeatStat.fpMinDetectTrain = zeros(1, Nregressions);
FeatStat.sensMinDetectTrain = zeros(1, Nregressions);
FeatStat.minDetectTrain = zeros(1, Nregressions);
% new to create ROC histograms
FeatStat.sensTrainROC = cell(1, Nregressions);
FeatStat.fpTrainROC = cell(1, Nregressions);


FeatStat.fpMinDetectTest = zeros(1, Nregressions);
FeatStat.sensMinDetectTest = zeros(1, Nregressions);
FeatStat.minDetectTest = zeros(1, Nregressions);
FeatStat.bFpIsZero = zeros(1, Nregressions); % 1 for is zero 
FeatStat.bFpTrainIsZero = zeros(1, Nregressions); % 1 for is zero 

FeatStat.AUC = zeros(1, length(regressionNames));
FeatStat.d_cohen = zeros(1, length(regressionNames));


% we'ere going to keep the raw data for anovas later
FeatStat.dpTrainScored = cell(1, Nregressions);
FeatStat.dnTrainScored = cell(1, Nregressions);
FeatStat.dpTestScored = cell(1, Nregressions);
FeatStat.dnTestScored = cell(1, Nregressions);

SNR = 2;
FeatStat.SNR = SNR;
NDnTest = sum(DNApDntest);
FeatStat.minDetectTestMax = (1+SNR*NDnTest.^-1)./(SNR*NDnTest.^-1);
FeatStat.fpTestMin = NDnTest.^-1;

zspace = linspace(-5, 5);
fp_no_use = normcdf(zspace);
sens_no_use = fp_no_use;


ii = 1;
for kk = 1:length(channelNames)
    jj = reorder(kk);
    %jj = kk;
    for ll = 1:4
    
    FeatStat.regressionNames = [channelNames{jj} '.' features{ll}];    
        
    [fp, sens, T, FeatStat.AUC(ii)] = compute_ROC(Dp.Channels.(channelNames{jj}).(features{ll})(DNApMCF7),...
        Dn.Channels.(channelNames{jj}).(features{ll})(DNAp));
    
    uptrain = mean(Dp.Channels.(channelNames{jj}).(features{ll})(DNApMCF7));
    untrain = mean(Dn.Channels.(channelNames{jj}).(features{ll})(DNAp));
    
    % if fp = 0 we didn't sample enough to measure it. set it equal
    % to thes sample size instead to prvent div 0 errors
    if sum(fp==0) > 0
        FeatStat.bFpTrainIsZero(ii) = 1;
    end
    
    fp(fp==0) = 1./length(Dn.Channels.(channelNames{jj}).(features{ll})(DNAp));      
    %sens = sens(fp~=0);
    %T = T(fp~=0);
    %fp = fp(fp~=0);    
    [FeatStat.minDetectTrain(ii), minDetInd] = max((sens+SNR*fp)./(SNR*fp));
    FeatStat.fpMinDetectTrain(ii) = fp(minDetInd);
    FeatStat.sensMinDetectTrain(ii) = sens(minDetInd);
    FeatStat.TMinDetect(ii) = T(minDetInd);
    FeatStat.d_cohen(ii) = abs(cohensD(Dp.Channels.(channelNames{jj}).(features{ll})(DNApMCF7), ...
        Dn.Channels.(channelNames{jj}).(features{ll})(DNAp)));
    FeatStat.regressionNames{ii} = [channelNames{jj}, '.', features{ll}];
    
    % new to create ROC historgrams :)
    FeatStat.sensTrainROC{ii} = sens;
    FeatStat.fpTrainROC{ii} = fp;
  
    %save data for later ANOVA analysis
    FeatStat.dpTrainScored{ii} = Dp.Channels.(channelNames{jj}).(features{ll})(DNApMCF7);
    FeatStat.dnTrainScored{ii} = Dn.Channels.(channelNames{jj}).(features{ll})(DNAp);
    FeatStat.dpTestScored{ii} = Dptest.Channels.(channelNames{jj}).(features{ll})(DNApDptest);
    FeatStat.dnTestScored{ii} = Dntest.Channels.(channelNames{jj}).(features{ll})(DNApDntest);
    
    % is the mean of disease positive more or less than disease negative
    if uptrain > untrain
        % diease positive is to right of disease negative
        sensTest = sum(Dptest.Channels.(channelNames{jj}).(features{ll})(DNApDptest) ...
            > FeatStat.TMinDetect(ii))./sum(DNApDptest);
        fpTest = sum(Dntest.Channels.(channelNames{jj}).(features{ll})(DNApDntest) ...
            > FeatStat.TMinDetect(ii))./sum(DNApDntest);
    else
        % diease positive is to left of disease negative
        sensTest = sum(Dptest.Channels.(channelNames{jj}).(features{ll})(DNApDptest) ...
            < FeatStat.TMinDetect(ii))./sum(DNApDptest);
        fpTest = sum(Dntest.Channels.(channelNames{jj}).(features{ll})(DNApDntest) ...
            < FeatStat.TMinDetect(ii))./sum(DNApDntest);
    end
        
    if fpTest == 0
        fpTest = 1/sum(DNApDntest);
        FeatStat.bFpIsZero(ii) = 1;
    end
    
    FeatStat.fpMinDetectTest(ii) = fpTest;
    FeatStat.sensMinDetectTest(ii) = sensTest;
    FeatStat.minDetectTest(ii) = (sensTest+SNR*fpTest)./(SNR*fpTest);    

    ii = ii + 1;
    end    
end

