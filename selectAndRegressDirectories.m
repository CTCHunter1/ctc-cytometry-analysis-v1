

mainPath = uigetdir(pwd, 'Data Head Level');

daysToProcess = {'day5', 'day7', 'day8', 'day9', 'day10', 'day12', 'day13', ...
    'day14', 'day15'};

% interate through each day and process the files
parfor ii = 1:length(daysToProcess)
    wbcFileName = [mainPath, filesep, daysToProcess{ii}, filesep, 'wbc_', daysToProcess{ii}, '.mat'];
    mcf7FileName = [mainPath, filesep, daysToProcess{ii}, filesep, 'MCF7_', daysToProcess{ii}, '.mat'];
    mixFileName = [mainPath, filesep, daysToProcess{ii}, filesep, 'mix_', daysToProcess{ii}, '.mat'];
    
    if(~exist(wbcFileName, 'file'))
        fprintf('WBC File Name: %s does not exist', wbcFileName);
    end
    
    if(~exist(mcf7FileName, 'file'))
        fprintf('MCF7 File Name: %s does not exist', mcf7FileName);
    end
    
    if(~exist(mixFileName, 'file'))
        fprintf('Mix File Name: %s does not exist', mixFileName);   
    end
    
    fprintf('Processing %s\n', daysToProcess{ii});
    selectAndRegress2(wbcFileName, mcf7FileName, mixFileName)
end
