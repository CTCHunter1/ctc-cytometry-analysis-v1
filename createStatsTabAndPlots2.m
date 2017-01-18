function [featuresFullData] = createStatsTabAndPlots2(featuresMinDet, pathname, filenameprefix, featLabels)
% featuresMinDet = cell array of performance statistics computed from determineFeatureStat
% pathname - path to where to save the figures
% filenameprefix - will put this prefix infront of the figure names
%                  useful to generate figures for different performance
%                  groups
% featlLabels - Seperate labels for the performance statistics. Used for
%               ANOVAs to determine source of experimental varation

% Copyright (C) 2016  Gregory L. Futia
% This work is licensed under a Creative Commons Attribution 4.0 International License.


% this will start by loading and restructuring the data, biggest part
% labels is expected to be the same length as featuresMinDet

% featuresMinDet is a  cell array
    
% field names are are AUC, d, sens. spec. Ndet
fieldNames = fieldnames(featuresMinDet{1});        
Nfields = length(fieldNames);
featuresFullData = struct;

% variables naames represent rows in the build out matrix

% go through each of the fields and build out 2D arrays for each of them
% one of the fields is the names of the field just copy that over
Nvariants = length(featuresMinDet);
Nregressions = length(featuresMinDet{1}.regressionNames);

regressionNames = featuresMinDet{1}.regressionNames;
featuresFullData.regressionNames = regressionNames;
% matlab field names can not have dots in them. We dot seperate the 
% channel and the feature
regressionNamesNoDot = regexprep(regressionNames, '\.', '_');
% matlab field names can not have {} in them. {} are used in the regressions
% they can't have + in them either % this will remove the _{} after
% regression names so it will work as a matlab field
regressionNamesNoDot = regexprep(regressionNamesNoDot, '_{[^\{]*$', '');

% we need to change the structure of the data for performing this analysis
% initalize the storage data set and copy over the data
for ii = 1:Nfields
    switch fieldNames{ii}
        case 'regressionNames'
            % there is nothing to do for this one
            % this case is to make an ROC histogram
        case {'sensTrainROC', 'fpTrainROC', 'sensTestROC', 'fpTestROC', ...
                'dpTrainScored', 'dnTrainScored', 'dpTestScored', 'dnTestScored', ...
                'dMixScored'}  
            bLabelVar = 0; % 0, don't create a label variable, 1 create
            
            % for these variables we also want to create a label variable
            % for the day the data came from and the technical varient
            % - these labels will allow us to use an ANOVA to determine if
            % between our within group variances are dominate
            switch fieldNames{ii}
                case { 'dpTrainScored', 'dnTrainScored', 'dpTestScored', 'dnTestScored'}
                    bLabelVar = 1; % 0, don't create a label variable, 1 create
            end
            
            % figure out how long the array will be
            featuresFullData.(fieldNames{ii}) = cell(1, Nregressions);
            for kk = 1:Nregressions
                % determine how many points we need to allociate for the
                % combined arrays
                Npts = 0;
                for jj=1:Nvariants
                    Npts = Npts + length(featuresMinDet{jj}.(fieldNames{ii}){kk});                          
                end

                % try to allocate the array here, it might not be square
                featuresFullData.(fieldNames{ii}){kk} = zeros(1, Npts);

                if bLabelVar == 1
                    featuresFullData.([fieldNames{ii} , 'TechLabel']){kk} = cell(1, Npts);
                    if exist('featLabels', 'var') == 1
                        featuresFullData.([fieldNames{ii} , 'DayLabel']){kk} = cell(1, Npts);                
                    end                
                end          

                                
                accum = 0;
                for jj = 1:Nvariants
                    NROC = length(featuresMinDet{jj}.(fieldNames{ii}){kk});                    
                    % one longfea array
                    featuresFullData.(fieldNames{ii}){kk}(accum+[1:NROC]) = ...
                        featuresMinDet{jj}.(fieldNames{ii}){kk};
                    
                    if bLabelVar == 1
                        for zz = 1:NROC
                            featuresFullData.([fieldNames{ii} , 'TechLabel']){kk}{accum+zz} = ...
                                int2str(jj);
                            if exist('featLabels', 'var') == 1
                                featuresFullData.([fieldNames{ii} , 'DayLabel']){kk}{accum+zz} = ...
                                    featLabels{jj};
                            end
                        end
                    end                    
                    accum = accum + NROC;
                end
            end
        case {'minDetectTestMax', 'fpTestMin', 'minDetectTrainMax', 'fpTrainMin'}
            featuresFullData.(fieldNames{ii}) = zeros(1, Nvariants);
            for kk=1:Nvariants
                featuresFullData.(fieldNames{ii})(kk) = featuresMinDet{kk}.(fieldNames{ii});
            end            
        otherwise         
            % this is for 1x1 scalar values           
            featuresFullData.(fieldNames{ii}) = zeros(Nvariants, Nregressions);
            % copy in the data
            for jj = 1:Nvariants
                featuresFullData.(fieldNames{ii})(jj,:) = featuresMinDet{jj}.(fieldNames{ii});        
            end

            % if labels exist perform anova
            if exist('featLabels', 'var') == 1
               for jj = 1:Nregressions            
                   [p, anovatab, stats] = anova1(featuresFullData.(fieldNames{ii})(:,jj), featLabels, 'off');
                   anovaData.(fieldNames{ii}).(regressionNamesNoDot{jj}).p = p;
                   anovaData.(fieldNames{ii}).(regressionNamesNoDot{jj}).anovatab = anovatab;
                   anovaData.(fieldNames{ii}).(regressionNamesNoDot{jj}).stats = stats;
               end
            end        
    end
end

% 1) save the anova results if features labels exist
% 2) average the ROC curves by day
% 3) group the ROC operating points by day
if exist('featLabels', 'var') == 1
    save([pathname, filenameprefix, 'anovaData.mat'], 'anovaData');
    % average the ROC curves by day
    identifiers = unique(featLabels); % get the unique identiifiers
    % if the identifiers are 'day#' reorder them by the #
    % this will work for all occurances of # after text for identifers
    numCellStr = regexp(identifiers, '\d*', 'match');
    % numcelLStr is cellstring array of length of the idenfiers
    if length(numCellStr) == length(identifiers)
        numCell = zeros(1, length(numCellStr));
        for ii = 1:length(numCellStr);
            numCell(ii) = str2num(numCellStr{ii}{1});
        end
        [~, newInd] = sort(numCell);
        identifiers = identifiers(newInd);
    end
    
    Nidentifiers = length(identifiers);
    NptsPerROC = 400;
    Npts = NptsPerROC*length(identifiers);
    % in the names strings
    for ii = 1:Nregressions
        % stitch together the averages from each of the days
        % iterate through the regression names and                       
        featuresFullData.sensTrainROCAvg{ii} = zeros(1, Npts);
        featuresFullData.fpTrainROCAvg{ii} = zeros(1, Npts);                     
        featuresFullData.sensTrainROCAvgMatrix{ii} = zeros(length(identifiers), NptsPerROC);        
        featuresFullData.fpTrainROCAvgMatrix{ii} = zeros(length(identifiers), NptsPerROC);
        featuresFullData.fpMinDetectTrainGrouped{ii} = cell(1, length(identifiers)); % this is jagged
        featuresFullData.sensMinDetectTrainGrouped{ii} = cell(1, length(identifiers));
        featuresFullData.fpMinDetectTrainGroupedAvg{ii} = cell(1, length(identifiers)); % this is jagged
        featuresFullData.sensMinDetectTrainGroupedAvg{ii} = cell(1, length(identifiers));
        
        
        if isfield(featuresMinDet{1}, 'sensTestROC')
        featuresFullData.sensTestROCAvg{ii} = zeros(1, Npts);
        featuresFullData.fpTestROCAvg{ii} = zeros(1, Npts);
        featuresFullData.sensTestROCAvgMatrix{ii} = zeros(length(identifiers), NptsPerROC);
        featuresFullData.fpTestROCAvgMatrix{ii} = zeros(length(identifiers), NptsPerROC);
        featuresFullData.fpMinDetectTestGrouped{ii} = cell(1, length(identifiers)); % this is jagged
        featuresFullData.sensMinDetectTestGrouped{ii} = cell(1, length(identifiers));
        featuresFullData.fpMinDetectTestGroupedAvg{ii} = cell(1, length(identifiers)); % this is jagged
        featuresFullData.sensMinDetectTestGroupedAvg{ii} = cell(1, length(identifiers));
        
        end
        
        for jj = 1:length(identifiers)
            indicies = strcmp(featLabels, identifiers{jj});
            % this is an array of featuresMinDet for each identifier
            featuresMinDetSub = featuresMinDet(indicies);
                        
            % first we need to put the sens and fp values into a 
            % cell array so they can be averaged
            senstrain_cell = cell(1, length(featuresMinDetSub));
            fptrain_cell = cell(1, length(featuresMinDetSub));
            featuresFullData.fpMinDetectTrainGrouped{ii}{jj} = zeros(1, length(featuresMinDetSub));
            featuresFullData.sensMinDetectTrainGrouped{ii}{jj} = zeros(1, length(featuresMinDetSub));
            
            if isfield(featuresMinDet{1}, 'sensTestROC')
                senstest_cell = cell(1, length(featuresMinDetSub));
                fptest_cell = cell(1, length(featuresMinDetSub));
                featuresFullData.fpMinDetectTestGrouped{ii}{jj} = zeros(1, length(featuresMinDetSub));
                featuresFullData.sensMinDetectTestGrouped{ii}{jj} = zeros(1, length(featuresMinDetSub));
            end
            
            for kk = 1:length(featuresMinDetSub);
                senstrain_cell{kk} = featuresMinDetSub{kk}.sensTrainROC{ii};
                fptrain_cell{kk} = featuresMinDetSub{kk}.fpTrainROC{ii};
                featuresFullData.fpMinDetectTrainGrouped{ii}{jj}(kk) = ...
                    featuresMinDetSub{kk}.fpMinDetectTrain(ii); 
                featuresFullData.sensMinDetectTrainGrouped{ii}{jj}(kk) = ...
                    featuresMinDetSub{kk}.sensMinDetectTrain(ii); 
                
                if isfield(featuresMinDet{1}, 'sensTestROC')
                    senstest_cell{kk} = featuresMinDetSub{kk}.sensTestROC{ii};
                    fptest_cell{kk} = featuresMinDetSub{kk}.fpTestROC{ii};
                    featuresFullData.fpMinDetectTestGrouped{ii}{jj}(kk) = ...
                    featuresMinDetSub{kk}.fpMinDetectTest(ii); 
                    featuresFullData.sensMinDetectTestGrouped{ii}{jj}(kk) = ...
                    featuresMinDetSub{kk}.sensMinDetectTest(ii); 
                end
            end
            [featuresFullData.fpTrainROCAvg{ii}(jj*NptsPerROC+[1:NptsPerROC]), ...
                featuresFullData.sensTrainROCAvg{ii}(jj*NptsPerROC+[1:NptsPerROC])] = ...
                averageROC(fptrain_cell, senstrain_cell, pi/4, NptsPerROC);  
            featuresFullData.sensTrainROCAvgMatrix{ii}(jj, :) = ...
                 featuresFullData.sensTrainROCAvg{ii}(jj*NptsPerROC+[1:NptsPerROC]);
            featuresFullData.fpTrainROCAvgMatrix{ii}(jj, :) = ...
                 featuresFullData.fpTrainROCAvg{ii}(jj*NptsPerROC+[1:NptsPerROC]);
            featuresFullData.fpMinDetectTrainGroupedAvg{ii}{jj} = ...
                mean(featuresFullData.fpMinDetectTrainGrouped{ii}{jj}); 
            featuresFullData.sensMinDetectTrainGroupedAvg{ii}{jj} = ...
                mean(featuresFullData.sensMinDetectTrainGrouped{ii}{jj}); 
             
             
            if isfield(featuresMinDet{1}, 'sensTestROC')
            [featuresFullData.fpTestROCAvg{ii}(jj*NptsPerROC+[1:NptsPerROC]), ...
                featuresFullData.sensTestROCAvg{ii}(jj*NptsPerROC+[1:NptsPerROC])] = ...
                averageROC(fptest_cell, senstest_cell, pi/4, NptsPerROC);         
            featuresFullData.sensTestROCAvgMatrix{ii}(jj, :) = ...
                featuresFullData.sensTestROCAvg{ii}(jj*NptsPerROC+[1:NptsPerROC]);
            featuresFullData.fpTestROCAvgMatrix{ii}(jj, :) =  ...
                featuresFullData.fpTestROCAvg{ii}(jj*NptsPerROC+[1:NptsPerROC]);
            featuresFullData.fpMinDetectTestGroupedAvg{ii}{jj} = ...
                mean(featuresFullData.fpMinDetectTestGrouped{ii}{jj}); 
            featuresFullData.sensMinDetectTestGroupedAvg{ii}{jj} = ...
                mean(featuresFullData.sensMinDetectTestGrouped{ii}{jj}); 
            end
        end
    end
end

% store the mean of the TCutParameter
% Op Cut for operating cut not optimal cut

% rename the labels
featuresFullData.regressionNames = regexprep(featuresFullData.regressionNames, '.totalSig_dBc', ' \\Sigma');
featuresFullData.regressionNames = regexprep(featuresFullData.regressionNames, '.M2', ' < M > ');
featuresFullData.regressionNames = regexprep(featuresFullData.regressionNames, '.radius_m', ' < r > ');
featuresFullData.regressionNames = regexprep(featuresFullData.regressionNames, '.radius_invm', ' < r_f > ');

%%
figSize = [.25, 2, 6.5/2, 3.75];

titlePos = [.26, .46];

figure(1);
clf;
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', figSize);
set(gcf, 'PaperUnits', 'Inches');
set(gcf, 'PaperPosition', figSize);

set(0, 'defaultTextFontSize', 7);
set(0, 'defaultAxesFontSize', 7);
fontSizeLabel = 8;
xshiftTick = .052;

switch filenameprefix
    case 'regression'
    fieldNamesForPlot = {'bFpIsZero', 'minDetectTest', 'fpMinDetectTest', ...
        'sensMinDetectTest', 'd_cohen', 'AUCTest', };
    case 'feature'
    fieldNamesForPlot = {'bFpIsZero', 'minDetectTest', 'fpMinDetectTest', ...
        'sensMinDetectTest', 'd_cohen', 'AUC', };
end

fieldNamesYLabel = {'FP Zero' 'N_{det}', 'Spec.', 'Sens.', ...
    'd', 'AUC', }; 

% TMinDetect needs to be scaled to the units of the feature
numFieldNamesForPlot = length(fieldNamesForPlot);

axh = zeros(1, numFieldNamesForPlot);

ao = .15;
aob = .28;
aot = .0;
aw = .75;
ah = .105;
ahscale = .73;


% need to rescale the thresholds on linear units
% this isn't a problem for the regressions
for ii = 1:numFieldNamesForPlot
axh(ii) = axes('position', [ao, aob+ ((ii-1)/numFieldNamesForPlot)*ahscale, aw, ah]);   

    switch fieldNamesForPlot{ii}
        case {'AUC', 'AUCTest', 'fpMinDetectTest', 'sensMinDetectTest', ...
                'fpMinDetectTrain', 'sensMinDetectTrain'}
        % put the AUC values on a z-scored axis
        switch fieldNamesForPlot{ii}
            case {'fpMinDetectTest', 'fpMinDetectTrain'}
            boxplot(axh(ii), norminv(1-featuresFullData.(fieldNamesForPlot{ii})));
            otherwise
            boxplot(axh(ii), norminv(featuresFullData.(fieldNamesForPlot{ii})));
        end
        yTicks = get(axh(ii), 'YTick');
        yTickLabels = sprintfc('%.3f', normcdf(yTicks));
       %yTickLabels = cellfun(@num2str, yTickLabels);
        xLim = get(axh(ii), 'Xlim');
        yLim = get(axh(ii), 'Ylim');
       
        set(axh(ii), 'YTickLabel', yTickLabels);
        axZ = axes('position', [ao, aob+ ((ii-1)/numFieldNamesForPlot)*ahscale, aw, ah]);
        h = plot(xLim, yLim);
        set(h, 'XData', []);
        set(h, 'YData', []);
        set(axZ, 'Color', 'none');
        set(axZ, 'Xlim', xLim);
        set(axZ, 'Ylim', yLim);
        set(axZ, 'YAxisLocation', 'right');
        set(axZ', 'XTick', []);
        h = ylabel(axZ, 'Z', 'Units', 'Normal', 'Rotation', 0);
        set(h, 'Position', [1.10 0.7 0]);
        axes(axh(ii));
                
        case 'minDetectTest'
        boxplot(axh(ii), featuresFullData.(fieldNamesForPlot{ii}));
        set(axh(ii), 'YScale', 'log');
        set(axh(ii), 'YTick', [10, 100, 1000]);
            %if kk ==1
            %set(axh(ii), 'Ylim', [1, 5000]);
            %else
            set(axh(ii), 'Ylim', [1, 5000]);   
            %end
        case {'bFpIsZero', 'bFpTrainIsZero'}
        plot(axh(ii), sum(featuresFullData.(fieldNamesForPlot{ii})), 'kx');  
        set(axh(ii), 'XLim', [.5, Nregressions+.5]);
        set(gca, 'XTick', [1:Nregressions]);
        otherwise
        boxplot(axh(ii), featuresFullData.(fieldNamesForPlot{ii}));
    end
    
    % label the variable names
    if(ii == 1)
    set(gca, 'XTickLabel', featuresFullData.regressionNames);
    %ylabel('Probability');
    set(gca, 'TickLabelInterpreter', 'Tex');
    set(gca, 'XTickLabelRotation', 90);
    % setting the x ticks resizes the axes, size it back
    else
    set(gca, 'XTickLabel', {});
    end

    ylabel(axh(ii), fieldNamesYLabel{ii}, 'Interpreter', 'Tex');
    
    annotation('line', [.03, .97], [.76, .76], 'Color', 'k', 'linewidth', .25);
    
% for soem unknown reason ploting in a axes resizes it
% lets just set it back to what it should be
set(gca, 'position', [ao, aob+ ((ii-1)/numFieldNamesForPlot)*ahscale, aw, ah]);

%save([matNames{kk}(1:end-4), 'CombinedData.mat'], 'featuresFullData');

end


% save the plots
saveas(gcf, [pathname, filenameprefix, 'TechVar', '.jpeg'], 'jpeg'); 
saveas(gcf, [pathname, filenameprefix, 'TechVar', '.eps'], 'epsc'); 
saveas(gcf, [pathname, filenameprefix, 'TechVar', '.png'], 'png'); 
saveas(gcf, [pathname, filenameprefix, 'TechVar', '.emf'], 'emf'); 



if exist('featLabels', 'var') == 1
    anovaFeat = {'p', 'f'};
    
    % anova features to plot
    for kk = 1:length(anovaFeat);    
        figure(1+kk);
        clf;
        clf;
        set(gcf, 'Units', 'Inches');
        set(gcf, 'Position', figSize);
        set(gcf, 'PaperUnits', 'Inches');
        set(gcf, 'PaperPosition', figSize);

        set(0, 'defaultTextFontSize', 7);
        set(0, 'defaultAxesFontSize', 7);
        fontSizeLabel = 8;
        xshiftTick = .052;


        axh = zeros(1, numFieldNamesForPlot);

        for ii = 1:numFieldNamesForPlot
        axh(ii) = axes('position', [ao, aob+ ((ii-1)/numFieldNamesForPlot)*ahscale, aw, ah]);   

            F = zeros(1, Nregressions);
            for jj = 1:Nregressions
                switch anovaFeat{kk}
                    case 'f'
                    F(jj) = anovaData.(fieldNamesForPlot{ii}).(regressionNamesNoDot{jj}).anovatab{2,5};
                    case 'p'
                    F(jj) = anovaData.(fieldNamesForPlot{ii}).(regressionNamesNoDot{jj}).p;
                end
            end

            plot(axh(ii), F, 'x');
            set(axh(ii), 'YScale', 'log');
            
            switch anovaFeat{kk}
                case 'p'
                hold(axh(ii), 'on');
                plot(axh(ii), [0, length(regressionNamesNoDot)+1], [.05, .05], 'k--')
                ylims = get(axh(ii), 'Ylim');
                set(axh(ii), 'YTick', [10^-10, 10^-5, 1]);                
                set(axh(ii), 'Ylim', [ylims(1), 2]);
            end
            
            set(axh(ii), 'XLim', [0, length(regressionNamesNoDot)+1]);

            % label the variable names
            if(ii == 1)
                set(gca, 'XTick', [1:Nregressions]);
                set(gca, 'XTickLabel', featuresFullData.regressionNames);
                %ylabel('Probability');
                set(gca, 'TickLabelInterpreter', 'Tex');
                set(gca, 'XTickLabelRotation', 90);
                % setting the x ticks resizes the axes, size it back
            else
            set(gca, 'XTickLabel', {});
            end

            ylabel(axh(ii), fieldNamesYLabel{ii}, 'Interpreter', 'Tex');

            annotation('line', [.03, .97], [.76, .76], 'Color', 'k', 'linewidth', .25);

        % for soem unknown reason ploting in a axes resizes it
        % lets just set it back to what it should be
        set(gca, 'position', [ao, aob+ ((ii-1)/numFieldNamesForPlot)*ahscale, aw, ah]);

        %save([matNames{kk}(1:end-4), 'CombinedData.mat'], 'featuresFullData');

        end

        % save the plots
        saveas(gcf, [pathname, filenameprefix, 'AnovaVar_', anovaFeat{kk}, '.jpeg'], 'jpeg'); 
        saveas(gcf, [pathname, filenameprefix, 'AnovaVar_', anovaFeat{kk}, '.eps'], 'epsc'); 
        saveas(gcf, [pathname, filenameprefix, 'AnovaVar_', anovaFeat{kk}, '.png'], 'png'); 
        saveas(gcf, [pathname, filenameprefix, 'AnovaVar_', anovaFeat{kk}, '.emf'], 'emf');     
    end
end

% this is a bifold size plot
% height is configured by number of regressions
% origionally designe for 16
% we need to build in a offset for the axis labels to 
% dynamically configure the height
Nvert = ceil(Nregressions/2);
figSize = [.25, 2, 6.5/2, (7.25/8)*Nvert+.75];

figure(4);
clf;
set(gcf, 'Units', 'Inches');
set(gcf, 'Position', figSize);
set(gcf, 'PaperUnits', 'Inches');
set(gcf, 'PaperPosition', figSize);

set(0, 'defaultTextFontSize', 7);
set(0, 'defaultAxesFontSize', 7);

axh = zeros(1, Nregressions);

% parameters for drawing the axes
ao = .12;
aob = .06;
aw = .4;
ah = .105;
ahscale = .87;
avscale = .83;
if Nregressions ==  16
    avspace = .007;
    ah = .102;
    aob = .09;
else
   avspace = 0;
   ah = .13;
   aob = .105;
end
    
aucZmin = -4;
aucZmax = 4;
zspace = linspace(aucZmin, aucZmax);
fp_no_use = normcdf(zspace);
sens_no_use = fp_no_use;

aucColormap = load('DensityPlotColorMapGray.mat');
thresholdMap = load('DensityPlotColorMapWhiteToRed.mat');
%8x2 figure
% two loops first loop is vertical, second loop is horizontal
% Nvert = NvertPan;
Nhorz = 2;
regIndex = 1;

ticks = [.00001, .0003, .005, .05, .2, .5, .8, .95, .995, .9997, .99999];
ticksZ = norminv(ticks);
tickLabel = cell(1, length(ticks));
for ii = 1:length(ticks);
    tickLabel{ii} = num2str(ticks(ii));
end

% use cloud plot or just overlay the plots
bCloudPlot = 0; % 1 use cloud plot, 0 overlay plots

switch filenameprefix
    case 'regression'
       fpFieldName = 'fpTestROC';
       sensFieldName = 'sensTestROC';
       threshFPFieldName = 'fpMinDetectTest';
       threshSensFieldName = 'sensMinDetectTest';
    case 'feature'
       fpFieldName = 'fpTrainROC';
       sensFieldName = 'sensTrainROC';
       threshFPFieldName = 'fpMinDetectTrain';
       threshSensFieldName = 'sensMinDetectTrain';       
end

% change between using averages or all data
if isfield(featuresFullData, 'fpTrainROCAvg')
    if bCloudPlot
        fpFieldName = [fpFieldName, 'Avg'];
        sensFieldName = [sensFieldName, 'Avg'];
    else
        fpFieldName = [fpFieldName, 'AvgMatrix'];
        sensFieldName = [sensFieldName, 'AvgMatrix'];
%         threshFPFieldName = [threshFPFieldName, 'Grouped'];
%         threshSensFieldName = [threshSensFieldName, 'Grouped'];
        threshFPFieldName = [threshFPFieldName, 'GroupedAvg'];
        threshSensFieldName = [threshSensFieldName, 'GroupedAvg'];
    end
end

% change loop direction to start at top instead of bottom
for ii=Nvert:-1:1
    for jj = 1:2
        % as we move between channels put a little bit of space
        % with the avspace parts of this
        axh(regIndex) = axes('position', [ao + ((jj-1)/Nhorz)*avscale, aob+ ((ii-1)/Nvert)*ahscale+avspace*(round((ii)/2)-1), aw, ah]);   
        %cloudPlot(norminv(featuresFullData.fpTrainROC{regIndex}), norminv(featuresFullData.sensTrainROC{regIndex}), [aucZmin aucZmax, aucZmin aucZmax], 0, [100, 100], aucColormap.map);
        
        if bCloudPlot || ~exist('identifiers', 'var')
            cloudPlot({norminv(featuresFullData.(fpFieldName){regIndex}), norminv(featuresFullData.(threshFPFieldName)(:, regIndex))}, {norminv(featuresFullData.(sensFieldName){regIndex}),norminv(featuresFullData.(threshSensFieldName)(:, regIndex))}, [aucZmin aucZmax, aucZmin aucZmax], 0, [200, 200], {aucColormap.map, thresholdMap.map}, {[0, 1], [0, 1]});
            hold on;
            plot(norminv(featuresFullData.(threshFPFieldName)(:, regIndex)).', ...
                 norminv(featuresFullData.(threshSensFieldName)(:, regIndex)).', 'r+', 'MarkerSize', 2);
            
        else
            plot(norminv(featuresFullData.(fpFieldName){regIndex}).', ...
                norminv(featuresFullData.(sensFieldName){regIndex}).');
            % // reset the color index
            hold on;
            ax = gca;
            ax.ColorOrderIndex = 1;
            % plot the operating points
            %plot(norminv(featuresFullData.(threshFPFieldName)(:, regIndex)).', ...
                %norminv(featuresFullData.(threshSensFieldName)(:, regIndex)).', '+', 'MarkerSize', 2);
            for kk = 1:Nidentifiers
                plot(norminv(featuresFullData.(threshFPFieldName){regIndex}{kk}), ...
                    norminv(featuresFullData.(threshSensFieldName){regIndex}{kk}), '+', 'MarkerSize', 4);
            end
            set(ax, 'Xlim', [aucZmin, aucZmax]);
            set(ax, 'Ylim', [aucZmin, aucZmax]);
        end
        
        plot(norminv(fp_no_use), norminv(sens_no_use), 'k--');
        set(axh(regIndex), 'XTick', ticksZ);
        set(axh(regIndex), 'YTick', ticksZ);
        set(axh(regIndex), 'XTickLabel', tickLabel);
        set(axh(regIndex), 'YTickLabel', tickLabel);
        
        
        if jj == 1
            hl = ylabel('Sens.'); 
            set(hl, 'Units', 'normalized');
            % hl is in data units
            % shift it up slightly
            pos = get(hl, 'Position');
            set(hl, 'Position', [pos(1)+.1, pos(2), pos(3)]);
        end
                
        % this will work for jagged Nregressiosn too
        if regIndex < Nregressions - 1
            set(axh(regIndex), 'XTickLabel', {});
        else
            set(axh(regIndex), 'XTickLabelRotation', 90);
            set(axh(regIndex), 'XTickLabelRotation', 90);
            hl = xlabel('1-Spec.');
            set(hl, 'Units', 'normalized');
            % hl is in data units
            % shift it up slightly
            pos = get(hl, 'Position');
            set(hl, 'Position', [pos(1), pos(2)+.1, pos(3)]);
        end
        % put z ticks on right axis
        if jj == 2
            set(axh(regIndex), 'YTickLabel', {});
            set(axh(regIndex), 'Box', 'off');
            axhtick = axes('position', get(axh(regIndex), 'position'), 'Nextplot', 'add', 'Box', 'off' );
            h = plot([aucZmin, aucZmax], [aucZmin, aucZmax]);
            set(h,'XData', []);
            set(h, 'YData', []);
            set(axhtick, 'Xlim', get(axh(regIndex), 'Xlim'));
            set(axhtick, 'Ylim', get(axh(regIndex), 'Ylim'));
            if ii == Nvert
                set(axhtick, 'XTick', [-3:3]);    
            else
                set(axhtick, 'XTick', get(axh(regIndex), 'XTick'));    
                set(axhtick, 'XTickLabel', {});
            end
            set(axhtick, 'YTick', [-3:3]);
            set(axhtick, 'YAxisLocation', 'right');
            set(axhtick, 'XAxisLocation', 'top');
            set(axhtick, 'Color', 'none');                            
        end
        
        % set the ticks on top of the top left axis
        if ii == Nvert && jj == 1
            set(axh(regIndex), 'Box', 'off');
            axhtick = axes('position', get(axh(regIndex), 'position'), 'Nextplot', 'add', 'Box', 'off' );
            h = plot([aucZmin, aucZmax], [aucZmin, aucZmax]);
            set(h,'XData', []);
            set(h, 'YData', []);
            set(axhtick, 'Xlim', get(axh(regIndex), 'Xlim'));
            set(axhtick, 'Ylim', get(axh(regIndex), 'Ylim'));
            set(axhtick, 'XTick', [-3:3]);    
            set(axhtick, 'YTick', get(axh(regIndex), 'YTick'));
            set(axhtick, 'YTickLabel', {});
            set(axhtick, 'YAxisLocation', 'right');
            set(axhtick, 'XAxisLocation', 'top');
            set(axhtick, 'Color', 'none');                            
        end
    
        % if identifiers exist 
        % add the legend 
        if jj==1 && ii==1 && exist('identifiers', 'var')
            h = columnlegend(3, identifiers, 'Orientation', 'Horizontal', 'Location', 'South');
            set(h, 'Position', [.15, -.075, .75, .12]);
        end
        
        %colormap(aucColormap.map);
        title(axh(regIndex), featuresFullData.regressionNames(regIndex), 'Units', 'Normal', 'Position', [.6, .05]); 
        if regIndex == Nregressions
            break;
        else            
            regIndex = regIndex + 1;       
        end
    end
end

annotation('textbox', [.92, .90 .1 .1], 'String', 'Z', 'EdgeColor', 'none');

% save the plots
saveas(gcf, [pathname, filenameprefix, 'ROCcloud', '.jpeg'], 'jpeg'); 
saveas(gcf, [pathname, filenameprefix, 'ROCcloud', '.eps'], 'epsc'); 
saveas(gcf, [pathname, filenameprefix, 'ROCcloud', '.png'], 'png'); 
saveas(gcf, [pathname, filenameprefix, 'ROCcloud', '.emf'], 'emf');  


if isfield(featuresFullData, 'dpTestScored') && ...
    isfield(featuresFullData, 'dpMixScored')
    figure(5);
    clf;
    CreateRegressionHistograms(featuresFullData.regressionNames, ...
        featuresFullData.dpTestScored, featuresFullData.dnTestScored, ...
        featuresFullData.dMixScored, featuresFullData.TMinDetect);

    % save the plots
    saveas(gcf, [pathname, filenameprefix, 'HistogramsTest', '.jpeg'], 'jpeg'); 
    saveas(gcf, [pathname, filenameprefix, 'HistogramsTest', '.eps'], 'epsc'); 
    saveas(gcf, [pathname, filenameprefix, 'HistogramsTest', '.png'], 'png'); 
    saveas(gcf, [pathname, filenameprefix, 'HistogramsTest', '.emf'], 'emf');      
end


bPeformStatTest = 1;

if strcmp(filenameprefix, 'regression') && bPeformStatTest
    PUT = {'minDetectTest', 'fpMinDetectTest', ...
    'sensMinDetectTest', 'd_cohen', 'AUCTest'};

fid = fopen([pathname, 'hypothesis_test_results.txt'], 'w');

% Test is Ndet is driven by sample size
p = ranksum(featuresFullData.minDetectTest(:,1), featuresFullData.minDetectTestMax);
if(p < .05)
    fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that Ndet of regressiosn 1 is the same as NdetMax\n', 100*p);
else
    fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that Ndet of regressiosn 1 is the same as NdetMax\n', 100*p);
end

% Test specificity is driven by sample size
p = ranksum(featuresFullData.fpMinDetectTest(:,1), featuresFullData.fpTestMin);
if(p < .05)
    fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that fpTest of regressiosn 1 is the same as fpMin\n', 100*p);
else
    fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that Ndet of fpTest of regressiosn 1 is the same as fpMin\n', 100*p);
end

% What if we look at the training data
% Ndet is driven by sample size
p = ranksum(featuresFullData.minDetectTrain(:,1), featuresFullData.minDetectTrainMax);
if(p < .05)
    fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that NdetTrain of regressiosn 1 is the same as NdetTrainMax\n', 100*p);
else
    fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that NdetTrain of regressiosn 1 is the same as NdetTrainMax\n', 100*p);
end

% Test specificity is driven by sample size
p = ranksum(featuresFullData.fpMinDetectTrain(:,1), featuresFullData.fpTrainMin);
if(p < .05)
    fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that fpTrain of regressiosn 1 is the same as fpTrainMin\n', 100*p);
else
    fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence fpTrain of regressiosn 1 is the same as fpTrainMin\n', 100*p);
end

% Ndet is driven by sample size
p = kruskalwallis([featuresFullData.minDetectTrain(:,1:4)- featuresFullData.minDetectTrainMax.'*ones(1,4)], [], 'off')
p = kruskalwallis([featuresFullData.minDetectTrain(:,1:4) featuresFullData.minDetectTrainMax.'], [], 'off');
if(p < .05)
    fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that NdetTrain of regressiosn 1-4 is the same as NdetTrainMax\n', 100*p);
else
    fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that NdetTrain of regressiosn 1-4 is the same as NdetTrainMax\n', 100*p);
end

fprintf(fid, '\n');

for ii=1:length(PUT)
    % test if Reg1_all and 3 feature Regs are equivalent
    switch PUT{ii}
        % try linearizing first
        case {'AUCTest', 'sensMinDetectTest', 'fpMinDetectTest'}
        %featuresFullData.(PUT{ii}) = norminv(featuresFullData.(PUT{ii}));
        otherwise        
    end
    
    [p, atab, stats] = kruskalwallis(featuresFullData.(PUT{ii})(:,1:4), [], 'off');

    if(p < .05)
        fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that the %s of regressiosn 1 through 4 are the same\n', 100*p, PUT{ii});
    else
        fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that the %s of regressiosn 1 through 4 are the same\n', 100*p,PUT{ii});
    end

    % test if 2 feature Regs are equivalent
    [p, atab, stats] = kruskalwallis(featuresFullData.(PUT{ii})(:,5:7), [], 'off');
    if(p < .05)
        fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that the %s of regressions with two channels are the same\n', 100*p, PUT{ii});
    else
        fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that the  %sof regressions with two channels are the same\n', 100*p, PUT{ii});
    end

    
    % test if 2 &1 feature Regs are equivalent
    [p, atab, stats] = kruskalwallis(featuresFullData.(PUT{ii})(:,5:11), [], 'off');
    if(p < .05)
        fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that the %s of regressions with two and one channels are the same\n', 100*p, PUT{ii});
    else
        fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that the % sof regressions with two and one channels are the same\n', 100*p, PUT{ii});
    end

    % test if 1 feature Regs are equivalent
    [p, atab, stats] = kruskalwallis(featuresFullData.(PUT{ii})(:,8:11), [], 'off');
    if(p < .05)
        fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that the %s of regressions with a single channel are the same\n', 100*p, PUT{ii});
    else
        fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that the %s of regressions with a single channel are the same\n', 100*p, PUT{ii});
    end

    % test if spatial features perform equal sum signal
    p = ranksum(featuresFullData.(PUT{ii})(:,1), featuresFullData.(PUT{ii})(:,2));
    if(p < .05)
        fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that the %s a Reg1_all and Reg2_sigma are the same\n', 100*p, PUT{ii});
    else
        fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that the %s a Reg1_all and Reg2_sigma are the same\n', 100*p, PUT{ii});
    end

    % test if Reg1_all beats Reg3_DAPI+Bodipy+CD45
    p = ranksum(featuresFullData.(PUT{ii})(:,1), featuresFullData.(PUT{ii})(:,3));
    if(p < .05)
        fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that the %s Reg1_all and Reg3_DAPI+CD45+Panck are the same\n', 100*p, PUT{ii});
    else
        fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that the %s Reg1_all and Reg3_DAPI+CD45+Panck are the same\n', 100*p, PUT{ii});
    end
    
    % test if DAPI+Bodipy+CD45 is better than DAPI+CD45+PanCK
    p = ranksum(featuresFullData.(PUT{ii})(:,3), featuresFullData.(PUT{ii})(:,4));
    if(p < .05)
        fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that the %s of Reg3_DAPI+PanCK+CD45 is the same as Reg4_DAPI+Bodipy+CD45\n', 100*p,PUT{ii});
    else
        fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that the  %s Reg3_DAPI+PanCK+CD45 is the same as Reg4_DAPI+Bodipy+CD45\n', 100*p, PUT{ii});
    end


    [h, p] = vartest2(featuresFullData.(PUT{ii})(:,3), featuresFullData.(PUT{ii})(:,4));
    if(p < .05)
        fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that variances of the %s for the Reg3_DAPI+PanCK+CD45 is the same as Reg4_DAPI+Bodipy+CD45\n', 100*p, PUT{ii});
    else
        fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that variances of the %s for a Reg3_DAPI+PanCK+CD45 is the same as Reg4_DAPI+Bodipy+CD45\n', 100*p, PUT{ii});
    end

    
    fprintf(fid, '\n');    
end
    fclose(fid);
end

fid = fopen([pathname, 'stats_', filenameprefix ,'_ndet.txt'], 'w');

for ii = 1:Nregressions
    fprintf(fid, 'Biomarker: %s, Ndet = %.2f (median), %.2f (min), %.2f (max) \n', ...
        featuresFullData.regressionNames{ii}, ...
        median(featuresFullData.minDetectTest(:, ii)), ...
        min(featuresFullData.minDetectTest(:, ii)), max(featuresFullData.minDetectTest(:, ii)));
    
    
end

fclose(fid);

fid = fopen([pathname, 'stats_', filenameprefix ,'_featurevarations.txt'], 'w');   

fprintf(fid, 'ANOVA results comparing technical and interday varation between groups\n');

for ii = 1:Nregressions
    fprintf('Regression %d \n', ii);
    [pDp] = anova1(featuresFullData.dpTrainScored{ii}, featuresFullData.dpTrainScoredTechLabel{ii}, 'off');   
    [pDn] = anova1(featuresFullData.dnTrainScored{ii}, featuresFullData.dnTrainScoredTechLabel{ii}, 'off');
    
    fprintf(fid, '%s D+ Tech: ', featuresFullData.regressionNames{ii});
    if(pDp < .05)
        fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that mean positions of D+ technical varients are the same\n', 100*pDp);
    else
        fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that mean positions of D+ technical varients are the same\n', 100*pDp);    
    end
    fprintf('\n');
    
    fprintf(fid, '%s D- Tech: ', featuresFullData.regressionNames{ii});
    if(pDn < .05)
        fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that mean positions of D- technical varients are the same\n', 100*pDn);
    else
        fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that mean positions of D- technical varients are the same\n', 100*pDn);    
    end
    
    if exist('featLabels', 'var') == 1
        [pDp] = anova1(featuresFullData.dpTrainScored{ii}, featuresFullData.dpTrainScoredDayLabel{ii}, 'off');   
        [pDn] = anova1(featuresFullData.dnTrainScored{ii}, featuresFullData.dnTrainScoredDayLabel{ii}, 'off');

        fprintf(fid, '%s D+ Day (1-way): ', featuresFullData.regressionNames{ii});
        if(pDp < .05)
            fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that mean positions of D+ days are the same\n', 100*pDp);
        else
            fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that mean positions of D+ days varients are the same\n', 100*pDp);    
        end
        
        fprintf(fid, '%s D- Day (1-way): ', featuresFullData.regressionNames{ii});
        if(pDn < .05)
            fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that mean positions of D- technical varients are the samem\n', 100*pDn);
        else
            fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence that mean positions of D- technical varients are the same\n', 100*pDn);    
        end 
        
        % do for 2 way anova
        [pDp] = anovan(featuresFullData.dpTrainScored{ii}, {featuresFullData.dpTrainScoredDayLabel{ii}, featuresFullData.dpTrainScoredTechLabel{ii}}, 'model', 2, 'varnames', {'days', 'tech'}, 'display', 'off');  
        [pDn] = anovan(featuresFullData.dnTrainScored{ii}, {featuresFullData.dnTrainScoredDayLabel{ii}, featuresFullData.dnTrainScoredTechLabel{ii}}, 'model', 2, 'varnames', {'days', 'tech'}, 'display', 'off'); 

        fprintf(fid, '%s D+ Day (2-way) Day: ', featuresFullData.regressionNames{ii});
        if(pDp(1) < .05)
            fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that mean positions of D+ Days are the same\n', 100*pDp(1));
        else
            fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence mean positions of D+ Days are the same\n', 100*pDp(1));    
        end
        
        fprintf(fid, '%s D+ Day (2-way) Tech: ', featuresFullData.regressionNames{ii});
        if(pDp(2) < .05)
            fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that mean positions of D+ Tech are the same\n', 100*pDp(2));
        else
            fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence mean positions of D+ Tech are the same\n', 100*pDp(2));    
        end
        
        fprintf(fid, '%s D+ Day (2-way) Day*Tech: ', featuresFullData.regressionNames{ii});
        if(pDp(3) < .05)
            fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence there is no interaction between day and tech in D+ \n', 100*pDp(3));
        else
            fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence there is no interaction between day and tech in D+\n', 100*pDp(3));    
        end
        
        fprintf(fid, '%s D- Day (2-way) Day: ', featuresFullData.regressionNames{ii});
        if(pDn(1) < .05)
            fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that mean positions of D- Days are the same\n', 100*pDn(1));
        else
            fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence mean positions of D- Days are the same\n', 100*pDn(1));    
        end
        
        fprintf(fid, '%s D- Day (2-way) Tech: ', featuresFullData.regressionNames{ii});
        if(pDn(2) < .05)
            fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence that mean positions of D- Tech are the same\n', 100*pDn(2));
        else
            fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence mean positions of D- Tech are the same\n', 100*pDn(2));    
        end
        
        fprintf(fid, '%s D- Day (2-way) Day*Tech: ', featuresFullData.regressionNames{ii});
        if(pDn(3) < .05)
            fprintf(fid, 'Reject null hypothesis with p=%.2f%% confidence there is no interaction between day and tech in D- \n', 100*pDn(3));
        else
            fprintf(fid, 'Fail to reject null hypothesis with p=%.2f%% confidence there is no interaction between day and tech in D- \n', 100*pDn(3));    
        end
               
    end
    fprintf('\n');
end

fclose(fid);
