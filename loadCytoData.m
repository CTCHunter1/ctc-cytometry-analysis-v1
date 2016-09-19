function [NDp, NDn, NDmix, DNACutOff, pathNameRet] = loadCytoData(inputFileNames)

if(exist('inputFileNames', 'var') == 1)
    
    filename = cell(1, 3);
    pathname = cell(1, 3);
    extension = cell(1, 3);

    for ii = 1:3
        [pathname{ii}, filename{ii}, extension{ii}] = fileparts(inputFileNames{ii});
        % put the extension back on the file name
        filename{ii} = [filename{ii}, extension{ii}];
        pathname{ii} = [pathname{ii}, filesep];
    end
end

% filename = {'wbc_all_data.mat', 'mcf7_all_data.mat', ... 
%     'wbc_mcf7_all_data.mat'};
% 
% pathname = {'C:\Repository\Projects\Nuclear Second Moment\Code\Processed Mats\', ...
%     'C:\Repository\Projects\Nuclear Second Moment\Code\Processed Mats\', ...
%     'C:\Repository\Projects\Nuclear Second Moment\Code\Processed Mats\'};

%nargin = 0;
% no inputs get the files from the user
if(exist('inputFileNames', 'var') == 0)  
    
    [filenameUI, pathnameUI] = uigetfile('*.mat', 'Select WBC file (Cancel for Default)', ...
        'C:\Repository\Projects\Nuclear Second Moment\Code\Processed Mats\');    

    if isequal(filenameUI,0) || isequal(pathnameUI,0)
        % user pressed cancel
    else
        filename{1} = filenameUI;
        pathname{1} = pathnameUI;
        
        [filenameUI, pathnameUI] = uigetfile('*.mat', 'Select Cancer Cell File', pathname{1});    
    
        if isequal(filenameUI,0) || isequal(pathnameUI,0)
            % user pressed cancel
        else
            filename{2} = filenameUI;
            pathname{2} = pathnameUI;
        
            [filenameUI, pathnameUI] = uigetfile('*.mat', 'Mixed Data File', pathname{2});
            
            if ~isequal(filenameUI,0) && ~isequal(pathnameUI,0)
                filename{3} = filenameUI;
                pathname{3} = pathnameUI;        
            end           
        end        
    end
end


% we want the folder name too if possible
%filesepPos = strfind(pathname{1}, filesep);
%date_str = pathname(filesepPos(end-1):end);

% assume this is coming in in a ND file
% pathname 1 is WBC = Diesease Negative
DatM = load([pathname{1}, filename{1}]);
Dat = DatM.ND;
NDn = DatM.ND;

% Pathname 2 is MCF7 = Diesease Positive
DatMCF7 = load([pathname{2}, filename{2}]);
NDp = DatMCF7.ND;

DatMix = load([pathname{3}, filename{3}]);
NDmix = DatMix.ND;

DNACutOff = 65;

pathNameRet = pathname{1};
