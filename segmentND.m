function NDarr = setmentND(NDin, numSplits)
% numSlits - now many sizes to split ND into
% NDin ND to be segmented
% NDarr - array of ND objects
% Modified from combine ND

if(nargin < 1)
    [savefilename, saveFilePath] = uigetfile('*.mat', 'Pick a Mat file', 'Processed Mats/Futia_Zenodo_2016/');    
    index = 1;

    NDin = load(fullfile(saveFilePath, savefilename));
    numSplits = 10;    
    
    
elseif nargin == 1
    savefilename = 'Combined Data.mat';
end

if strcmp(fieldnames(NDin), 'ND')
    NDin = NDin.ND;
end

% generate random indexes for the ROIs to split the data into
splitInd = randi(numSplits, NDin.numROIs, 1);
% allocate ND arrar
NDarr = cell(1, numSplits); 

% create the heaade info int the ND structrues
for ii = 1:numSplits
    NDarr{ii}.fileName = sprintf('%s_split%i.mat', savefilename(1:end-4), ii);
    NDarr{ii}.numROIs = sum(splitInd == ii);
    NDarr{ii}.dx = NDin.dx;
    NDarr{ii}.dy = NDin.dy;
end
% determine number of points
Npts = 0;

channelNames = fieldnames(NDin.Channels);
NChannels = length(channelNames);

% interate though each channel
for ii = 1:NChannels
    featureNames = fieldnames(NDin.Channels.(channelNames{ii}));
    % and then each metric
    for jj = 1:length(featureNames)
        % scalars 
        if(length(NDin.Channels.(channelNames{ii}).(featureNames{jj})) == 1)
            for kk = 1:numSplits
                NDarr{kk}.Channels.(channelNames{ii}).(featureNames{jj}) = ...
                    NDin.Channels.(channelNames{ii}).(featureNames{jj});
            end
        % data
        else
            % populate the channels
            for kk = 1:numSplits            
                NDarr{kk}.Channels.(channelNames{ii}).(featureNames{jj}) = ...
                    NDin.Channels.(channelNames{ii}).(featureNames{jj})(splitInd == kk);         
            end
        end
    end            
end

% save the new ND file
if nargin < 1
    for ii = 1:numSplits
        ND = NDarr{ii};
        save(fullfile(saveFilePath, ND.fileName), 'ND');
    end
end
