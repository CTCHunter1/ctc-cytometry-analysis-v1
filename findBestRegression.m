function [beta, mu, sig, features, varnames, units, dAll] = findBestRegression(regressionAllFilename, selectedChannels, selectedFeatures)
%findBestRegression Looks through regression results for the best 
% subregression

    Nchannels = length(selectedChannels);
    Nfeatures = length(selectedFeatures);
    
    selectedElements = cell(1, Nchannels*Nfeatures);
    % build the features list
    index = 1;
    for ii = 1:Nchannels
        for jj = 1:Nfeatures
            selectedElements{index} = [selectedChannels{ii}, '.', selectedFeatures{jj}];
            index = index + 1;
        end
    end

    load(regressionAllFilename);
    
    % the regressionAllfile contains the results from all regressions
    % searched in that file. Search through them for ones that contain
    % elements from selectedElements cell list. Ones with elements outside
    % that list should be excluded
    
    dRegMax = 0;
    maxNumFeatures = 0;
    maxComb = 0;
    
    loopIterations = 0;
    %figure out how many loops
    for ii=1:length(featuresAll)
        loopIterations = loopIterations + length(featuresAll{ii});
    end
    dAll = zeros(1, loopIterations);
    iDAll = 1;
    
    % first loop is the Nth feature combination
    for ii = 1:length(featuresAll)
        % the second loop is combinations of those ii features
        for jj = 1:length(featuresAll{ii})
            regFeatures = featuresAll{ii}{jj};
            
            % regFeature needs to ba subset of selected elements
            % or we skip it, check this
            bcSum = 0;
            for kk =1:length(regFeatures);
                bC = strcmp(regFeatures{kk}, selectedElements);
                % there will be a 1 in that list if ther is a match
                bcSum = sum(bC) + bcSum;
            end
            
            % regFeatures only has features that are in selected elements
            if bcSum == length(regFeatures)
               dAll(iDAll) = dReg{ii}(jj);
               iDAll = iDAll + 1;
               if dReg{ii}(jj) > dRegMax 
                  maxNumFeatures = ii;
                  maxComb = jj;
                  dRegMax = dReg{ii}(jj);
               end
               
            end            
        end
    end
    beta = betaAll{maxNumFeatures}{maxComb};
    mu = muAll{maxNumFeatures}{maxComb};
    sig = sigAll{maxNumFeatures}{maxComb};
    features = featuresAll{maxNumFeatures}{maxComb};
    varnames = varnamesAll{maxNumFeatures}{maxComb};
    units = unitsAll{maxNumFeatures}{maxComb};    
    dAll = dAll(dAll ~= 0);
end

