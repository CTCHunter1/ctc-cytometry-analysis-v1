function [scores] = data_to_regressed_scored(ND, regressionFilename)
% takes a data structure and scores it along a predone regression loaded
% from regression.mat

if nargin < 2
    regressionName = 'Reg1_{all}.mat';
end

% these constants are not needed
channelNames = fieldnames(ND.Channels);
NchannelNames = length(channelNames);
selectedFeatures = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
varEqn = {'\Sigma_', '<r>_', '<r_f>_', '<M>_'};
unitsEqn = {' dBct ', ' \mu m ', ' \mu m^{-1} ', ''};


% tack .mat onto regression file name if it isn't already there
if strcmp(regressionFilename(end-4:end), '.mat') ~= 0
    regressionFilename = [regressionFilename, '.mat'];
end
% the key to this code is building up the allData matrix 
% the rows of the all data matrix are in order of the string
% 'features' listed in the regressed mat
% other loaded parameters, z_REG, mu_REG, sigma_REG are of the length
% of the number of selected variables (varIndex)

load(regressionFilename);
% the updated version doesn't save varIndex just the feature names
% in this case just make this everything
if ~exist('varIndex', 'var')
    varIndex = 1:length(features);
end
    
numFeatures = length(features); % features is variable in the .mat file
allData = zeros(numFeatures, ND.numROIs);

% the front of the numFeatures string is the channel
% the second part is the feature example: DAPI.radius_m

% perform search and build
for ii = 1:numFeatures
    [channel, feat] = strtok(features{ii}, '.');
    feat = feat(2:end); % the period comes with it
    % I was going to do 2 for loops for the search
    % this is so much cleaner
    % if the feature doesn't exist this will throw an error; I think
    allData(ii, :) = ND.Channels.(channel).(feat);
end
 

[Nd, Md] = size(allData.');

scores = [(allData(varIndex,:) - mu_REG.'*ones(1,Nd))./(sigma_REG.'*ones(1,Nd))].'*betareg(2:end);

