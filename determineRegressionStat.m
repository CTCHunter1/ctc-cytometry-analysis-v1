function [FeatStat] = determineRegressionStat(NDAll, DNACutOff, regressionFolder)
% overloading in Matlab isn't so easy, changed to bring ND in
% as 1x4 or 1x5 cell array so mixed data could be added
% to the regStats
Ndata = length(NDAll); % this is a 1x4 or 1x5 cell array

switch Ndata
    case 4
        Dp = NDAll{1};
        Dn = NDAll{2};
        Dptest = NDAll{3};
        Dntest = NDAll{4};
    case 5
        Dp = NDAll{1};
        Dn = NDAll{2};
        Dptest = NDAll{3};
        Dntest = NDAll{4};
        Dmix = NDAll{5};
end


% Created By: G. Futia
% modified from determine FeatureStats

if nargin == 0
    % this script is ment to be called by create_manuscript_figures
    % if this script is ran without parameters run that script instead
    estimateTechnicalVariances; 
    return;
end

filesepchar = filesep();

% tack on a slash to the end of regression folder if it doesn't have it
if(regressionFolder(end) ~= filesepchar)
    regressionFolder = [regressionFolder, filesepchar];
end

% determine the number of regression_names in the regression folder
% this folder is assumed to be filled of .mat files with regression
% parameters in them
bGetNamesFolder = 0;

if(bGetNamesFolder)
    dirResult = dir([regressionFolder, '*.mat']);
    
    regressionNames = cell(1, length(dirResult));
    % one of these files probably has a _all.mat suffix that we want to exclude
    regIndex = 1;
    Nregressions = 0;
    for ii=1:length(dirResult)
        if isempty(strfind(dirResult(ii).name, '_all.mat'))
            regressionNames{regIndex} = dirResult(ii).name;
            regIndex = regIndex + 1;
            Nregressions = Nregressions + 1;
        end
    end

    regressionNames = regressionNames(1:Nregressions);
    % remove the .mat from end of regression names
    regressionNames = regexprep(regressionNames, '.mat', '');
else
    % this is just the easiest way to do it
    % the folder dir works but then we have to deal with ordering them
    regressionNames = {'Reg1_{All}', 'Reg2_{Sigma}', 'Reg3_{DAPI + CD45 + PanCK}', ...
        'Reg4_{DAPI + Bodipy + CD45}', 'Reg5_{DAPI + CD45}', 'Reg6_{DAPI + PanCK}', ...
        'Reg7_{DAPI + Bodipy}', 'Reg8_{DAPI}', 'Reg9_{Bodipy}', 'Reg10_{CD45}', 'Reg11_{PanCK}'};
    Nregressions = length(regressionNames);
end
    
    
% Set the DNA Gates
% DNAp stands for DNA positive
DNAp = Dn.Channels.DAPI.totalSig_dBc > DNACutOff;
DNApMCF7 = Dp.Channels.DAPI.totalSig_dBc > DNACutOff;
DNApDptest = Dptest.Channels.DAPI.totalSig_dBc > DNACutOff;
DNApDntest = Dntest.Channels.DAPI.totalSig_dBc > DNACutOff;


FeatStat.regressionNames = regressionNames;
FeatStat.TMinDetect = zeros(1, Nregressions);
% these will be set using the test data
FeatStat.fpMinDetectTrain = zeros(1, Nregressions);
FeatStat.sensMinDetectTrain = zeros(1, Nregressions);
FeatStat.minDetectTrain = zeros(1, Nregressions);
% to create ROC histograms
FeatStat.sensTrainROC = cell(1, Nregressions);
FeatStat.fpTrainROC = cell(1, Nregressions);

FeatStat.sensTestROC = cell(1, Nregressions);
FeatStat.fpTestROC = cell(1, Nregressions);
FeatStat.AUCTest = zeros(1, Nregressions);

FeatStat.fpMinDetectTest = zeros(1, Nregressions);
FeatStat.sensMinDetectTest = zeros(1, Nregressions);
FeatStat.minDetectTest = zeros(1, Nregressions);
FeatStat.bFpIsZero = zeros(1, length(regressionNames)); % 1 for is zero 
FeatStat.bFpTrainIsZero = zeros(1, length(regressionNames)); % 1 for is zero 


FeatStat.AUC = zeros(1, Nregressions);
FeatStat.d_cohen = zeros(1, Nregressions);
FeatStat.d_cohenTrain = zeros(1, Nregressions);

% we'ere going to keep the scored data to make a histogram with
FeatStat.dpTrainScored = cell(1, Nregressions);
FeatStat.dnTrainScored = cell(1, Nregressions);
FeatStat.dpTestScored = cell(1, Nregressions);
FeatStat.dnTestScored = cell(1, Nregressions);

if exist('Dmix', 'var')
   DNApMix = Dmix.Channels.DAPI.totalSig_dBc > DNACutOff;
   FeatStat.dMixScored = cell(1, Nregressions);   
end

SNR = 2;

zspace = linspace(-5, 5);
fp_no_use = normcdf(zspace);
sens_no_use = fp_no_use;

% iterate through each regression
for ii = 1:Nregressions
    regressionFullFile = [regressionFolder, regressionNames{ii}, '.mat'];
    dpScored = data_to_regressed_scored(Dp, regressionFullFile);
    dnScored = data_to_regressed_scored(Dn, regressionFullFile);
    dpTestScored = data_to_regressed_scored(Dptest, regressionFullFile);
    dnTestScored = data_to_regressed_scored(Dntest, regressionFullFile);
       
    if exist('Dmix', 'var')    
        dMixScored = data_to_regressed_scored(Dmix, regressionFullFile);
        dMixScored = dMixScored(DNApMix).';
        FeatStat.dMixScored{ii} = dMixScored;
    end
        
    % remove data that isn't DNA positive
    dpScored = dpScored(DNApMCF7).';
    dnScored = dnScored(DNAp).';
    dpTestScored = dpTestScored(DNApDptest).';
    dnTestScored = dnTestScored(DNApDntest).';
    
    % we could do this asignment directly. we orgionally weren't saving the
    % data
    FeatStat.dpTrainScored{ii} = dpScored;
    FeatStat.dnTrainScored{ii} = dnScored;
    FeatStat.dpTestScored{ii} = dpTestScored;
    FeatStat.dnTestScored{ii} = dnTestScored;
    
    [fp, sens, T, FeatStat.AUC(ii)] = compute_ROC(dpScored,...
        dnScored);
    % compute ROC corve on the test data set as well. Do this
    % because test data set is not what the regression is trained on
    % expecting this to be worse
    [fpTest, sensTest, ~, FeatStat.AUCTest(ii)] = compute_ROC(dpTestScored,...
        dnTestScored);
    
    
    uptrain = mean(dpScored);
    untrain = mean(dnScored);
    
    % if fp = 0 we didn't sample enough to measure it. set it equal
    % to these sample size instead to prvent div 0 errors
    if sum(fp==0) > 0
        FeatStat.bFpTrainIsZero(ii) = 1;
    end
        
    fp(fp==0) = 1./length(dnScored);      
    %sens = sens(fp~=0);
    %T = T(fp~=0);
    %fp = fp(fp~=0);    
    [FeatStat.minDetectTrain(ii), minDetInd] = max((sens+SNR*fp)./(SNR*fp));
    FeatStat.fpMinDetectTrain(ii) = fp(minDetInd);
    FeatStat.sensMinDetectTrain(ii) = sens(minDetInd);
    FeatStat.TMinDetect(ii) = T(minDetInd);
    FeatStat.d_cohen(ii) = abs(cohensD(dpTestScored, ...
        dnTestScored));
    FeatStat.d_cohenTrain(ii) = abs(cohensD(dpScored, ...
        dnScored));
    
    FeatStat.sensTrainROC{ii} = sens;
    FeatStat.fpTrainROC{ii} = fp;
    
    FeatStat.sensTestROC{ii} = sensTest;
    FeatStat.fpTestROC{ii} = fpTest;
    
    % is the mean of disease positive more or less than disease negative
    if uptrain > untrain
        % diease positive is to right of disease negative
        sensTest = sum(dpTestScored ...
            > FeatStat.TMinDetect(ii))./sum(DNApDptest);
        fpTest = sum(dnTestScored ...
            > FeatStat.TMinDetect(ii))./sum(DNApDptest);
    else
        % diease positive is to left of disease negative
        sensTest = sum(dpTestScored ...
            < FeatStat.TMinDetect(ii))./sum(DNApDptest);
        fpTest = sum(dnTestScored ...
            < FeatStat.TMinDetect(ii))./sum(DNApDptest);
    end
        
    if fpTest == 0
        fpTest = 1/sum(DNApDptest);
        FeatStat.bFpIsZero(ii) = 1;
    end
    
    FeatStat.fpMinDetectTest(ii) = fpTest;
    FeatStat.sensMinDetectTest(ii) = sensTest;
    FeatStat.minDetectTest(ii) = (sensTest+SNR*fpTest)./(SNR*fpTest);    

    ii = ii + 1;
    end    
end

