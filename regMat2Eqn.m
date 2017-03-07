% Creates a .tex equation from each of the mat regression files in a
% directory

dirName = uigetdir('Processed Mats\Futia_Zenodo_2016\day9\Reg Mats', 'Select a directory');

D = dir([dirName, filesep, '*.mat']);

Nfiles = length(D);

for ii = 1:Nfiles
    regDat = load([dirName, filesep, D(ii).name]);
    regName = D(ii).name(1:end-4);
    
    % make sure this isn't some other kind of mat file
    if isfield(regDat, 'features')
        if ~isfield(regDat, 'varnames') && isfield(regDat, 'varIndex');
            channels = cell(1, length(regDat.varIndex));
            metric = cell(1, length(regDat.varIndex));
            regDat.varnames = cell(1, length(regDat.varIndex));
            
            for ii = 1:length(regDat.varIndex)
               split = strsplit(regDat.features{regDat.varIndex(ii)}, '.');
               channels{ii} = split{1};
               metric{ii} = split{2};           
            end

            features1 = {'totalSig_dBc', 'radius_m', 'radius_invm','M2'};
            % split up the feature strings
            
            regDat.units = regexprep(metric, features1, {' dBct ', ' \mu m ', ' \mu m^{-1} ', ''});
            metric = regexprep(metric, features1, {'\Sigma', '<r>', '<r_f>', '<M>'});            

            for ii = 1:length(regDat.varIndex)
                regDat.varnames{ii} = [metric{ii}, '_{', channels{ii}, '}'];
            end        
        end

        % make Sigma latex version
        regDat.varnames = strrep(regDat.varnames, 'Sigma', '\Sigma');
        regDat.units = strrep(regDat.units, 'mu', '\mu');
        regName = strrep(regName, 'Sigma', '\Sigma');

        fprintf('On Regression: %s\n', regName);   
        % save the regression equation
        writePC1Eqn([dirName, filesep, strrep(regName, '\', ''), '.tex'], ['\rm{', regName,  '}'], ...
            regDat.varnames, regDat.units, regDat.mu_REG, regDat.sigma_REG, regDat.betareg(2:end));             
    end
end



   