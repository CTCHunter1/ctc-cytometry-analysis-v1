function ND = combineND(data)

if(nargin < 1)
    [filenames{1}, pathnames{1}] = uigetfile('*.mat', 'Pick a Mat file', 'C:\Repository\Projects\Nuclear Second Moment\Code\Processed Mats');    
    index = 1;

    while ~isequal(filenames{index},0) || ~isequal(pathnames{index},0)
        index = index + 1;
        [filenames{index}, pathnames{index}] = uigetfile('*.mat', 'Pick a Image file', pathnames{index-1});        
    end

    NumFiles = index - 1; % the last cell is empty

    if NumFiles == 0
        return;
    end

    % select file name to save as
    [savefilename, savepathname] = uiputfile('*.mat', 'Save .mat as', pathnames{index-1});    
    if isequal(savefilename, 0)
        return;
    end
    data = cell(1, NumFiles);
elseif nargin == 1
    savefilename = 'Combined Data.mat';
end

NumFiles = length(data);


% determine number of points
Npts = 0;

% determine how many data points there are
for ii = 1:NumFiles;    
    if nargin < 1
        data{ii} = load(fullfile(pathnames{ii}, filenames{ii}));
    end
    if isfield(data{ii}, 'ND')       
        Npts = Npts + data{ii}.ND.numROIs;
        NptsFile(ii) = data{ii}.ND.numROIs; 
    else
        Npts = Npts + data{ii}.numROIs;
        NptsFile(ii) = data{ii}.numROIs;
        data{ii}.ND = data{ii};
    end
end

ND.fileName = savefilename;
ND.numROIs = Npts;
ND.dx = data{ii}.ND.dx;
ND.dy = data{ii}.ND.dy;
channelNames = fieldnames(data{ii}.ND.Channels);
NChannels = length(channelNames);
pointsIndex = cumsum(NptsFile);
pointsIndex = [0, pointsIndex];

% interate though each channel
for ii = 1:NChannels
    featureNames = fieldnames(data{1}.ND.Channels.(channelNames{ii}));
    for jj = 1:length(featureNames)
        if(length(data{1}.ND.Channels.(channelNames{ii}).(featureNames{jj})) == 1)
            ND.Channels.(channelNames{ii}).(featureNames{jj}) = ...
                data{1}.ND.Channels.(channelNames{ii}).(featureNames{jj});
        else
            % populate the channels
            ND.Channels.(channelNames{ii}).(featureNames{jj}) = zeros(1, Npts);
            for kk = 1:NumFiles;
                ND.Channels.(channelNames{ii}).(featureNames{jj})(pointsIndex(kk)+[1:NptsFile(kk)]) = ...
                    data{kk}.ND.Channels.(channelNames{ii}).(featureNames{jj});
                    
            end
        end
    end            
end

% save the new ND file
if nargin < 1
save(fullfile(savepathname, savefilename), 'ND');
end
