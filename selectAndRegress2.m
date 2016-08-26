% Created by: G. Futia
function [] = selectAndRegress2(wbcFileName, mcf7FileName, savepathmat, savepatheqn)

if ~exist('savepathmat', 'var')
    if(exist('wbcFileName', 'var'))
        savePath = fileparts(wbcFileName);
    else 
        savePath = '';
    end
    savepathmat = [savePath, filesep, 'Reg Mats', filesep];
    savepatheqn = [savePath, filesep, 'Reg Equations', filesep];
end

% make the directories if they don't already exist
if(~exist(savepathmat(1:end-1), 'dir'))
    mkdir(savepathmat(1:end-1));
end 

if(~exist(savepatheqn(1:end-1), 'dir'))
    mkdir(savepatheqn(1:end-1));
end 

bTest = 0;
if bTest
%changing regression names
selectedChannels{1} = {'DAPI'};
selectedFeatures{1} = {'totalSig_dBc', 'M2'};
regName{1} = 'Regtest_{All}';

selectedChannels{2} = {'DAPI'};
selectedFeatures{2} = {'totalSig_dBc'};
regName{2} = 'Reg2_{Sigma}';
else 
selectedChannels{1} = {'DAPI', 'Bodipy', 'PanCK', 'CD45'};
selectedFeatures{1} = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
regName{1} = 'Reg1_{All}';

% V.2
% instead of re running the regressiosn we will search the the results from
% the first regression and figure out what the best regression was for
% these other variables

selectedChannels{2} = {'DAPI', 'Bodipy', 'PanCK', 'CD45'};
selectedFeatures{2} = {'totalSig_dBc'};
regName{2} = 'Reg2_{Sigma}';

selectedChannels{3} = {'DAPI', 'PanCK', 'CD45'};
selectedFeatures{3} = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
regName{3} = 'Reg3_{DAPI + CD45 + PanCK}';

selectedChannels{4} = {'DAPI', 'Bodipy', 'CD45'};
selectedFeatures{4} = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
regName{4} = 'Reg4_{DAPI + Bodipy + CD45}';

selectedChannels{5} = {'DAPI', 'CD45'};
selectedFeatures{5} = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
regName{5} = 'Reg5_{DAPI + CD45}';

selectedChannels{6} = {'DAPI', 'PanCK'};
selectedFeatures{6} = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
regName{6} = 'Reg6_{DAPI + PanCK}';

selectedChannels{7} = {'DAPI', 'Bodipy'};
selectedFeatures{7} = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
regName{7} = 'Reg7_{DAPI + Bodipy}';

selectedChannels{8} = {'DAPI'};
selectedFeatures{8} = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
regName{8} = 'Reg8_{DAPI}';

selectedChannels{9} = {'Bodipy'};
selectedFeatures{9} = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
regName{9} = 'Reg9_{Bodipy}';

selectedChannels{10} = {'CD45'};
selectedFeatures{10} = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
regName{10} = 'Reg10_{CD45}';

selectedChannels{11} = {'PanCK'};
selectedFeatures{11} = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
regName{11} = 'Reg11_{PanCK}';
end

regressionAllFileName = [savepathmat strrep(regName{1}, '\', ''), '_all.mat']; 

numRegresions = length(regName);

fprintf('On Regression: %s\n', regName{1});

% this will create the regression All file name used 
% the for loop with sort though this file for the best regression
feature_selection(selectedChannels{1}, selectedFeatures{1}, regName{1}, {wbcFileName, mcf7FileName}, savepathmat);

for ii = 1:numRegresions
    regressionFile = [savepathmat, strrep(regName{ii}, '\', ''), '.mat'];
    
    fprintf('On Regression: %s\n', regName{ii});
    [betareg, mu_REG, sigma_REG, features, varnames, units, dAll] = findBestRegression(regressionAllFileName, selectedChannels{ii}, selectedFeatures{ii});
    % save the regression results
    save(regressionFile, 'betareg', 'mu_REG', 'sigma_REG', 'features', 'varnames', 'units', 'dAll');        
    % save the regression equation
    writePC1Eqn([savepatheqn, strrep(regName{ii}, '\', ''), '.tex'], ['\rm{', regName{ii},  '}'], varnames, units, mu_REG, sigma_REG, betareg(2:end));             
end



   