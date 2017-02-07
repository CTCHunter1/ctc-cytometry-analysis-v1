function [] = CreateRegressionHistograms(regressionNames, scores_Dp, scores_Dn, scores_mixed, thresholds )
%CreateRegressionHistograms Creates Histograms of the regression names
% regressionNames is size 1xNregressions
% scoresDp is cell array size 1xNregressions
% scoresDn is cell array size 1xNregreesions 
% scoresMixed is cell array size 1XNregressions
% Assume Figure was called before entering this function. We don't call
% call saveas after to save the figure

numBins = 100;

ao = .12;
aob = .05;
aot = .0;
aw = .84;
ah = .065;
ahscale = .95;

figSize = [.25, 2, 6.5/2, 5.5];

titlePos = [.26, .50];

set(gcf, 'Units', 'Inches');
set(gcf, 'Position', figSize);
set(gcf, 'PaperUnits', 'Inches');
set(gcf, 'PaperPosition', figSize);

set(0, 'defaultTextFontSize', 7);
set(0, 'defaultAxesFontSize', 7);
fontSizeLabel = 8;
xshiftTick = .052;


boxPlotYoffset = .38;

% histograms for regressions, 
% this was for a bug fix to compare histograms
numRegresions = length(regressionNames);
reorder = numRegresions:-1:1;
for ii=1:length(regressionNames)    
    % reindex off of dsort. Plot in order of separation
    jj = reorder(ii);    
    
    % set the position of the axis
    axes('position', [ao, aob+ ((ii-1)/numRegresions)*ahscale, aw, ah]);   
    % switch for the different data types
    [n1, xout] = hist(scores_Dp{jj}, numBins);
    [n1, xout] = makejagged(n1, xout); 
    plot(xout, n1, 'r'); hold on;    
    [n2, xout2] = hist(scores_Dn{jj}, numBins);
    [n2, xout2] = makejagged(n2, xout2); 
    n2 = max(n1)*n2/max(n2);
    plot(xout2, n2, 'b'); hold on;    
    [n3, xout3] = hist(scores_mixed{jj}, numBins);
    [n3, xout3] = makejagged(n3, xout3); 
    n3 = max(n1)*n3/max(n3);
    plot(xout3, n3, 'k--'); hold on;

    %set(gca, 'YTick', []);                
    xlim = [-1.5, 1.5];
    set(gca, 'Xlim', xlim);
    a = title(regressionNames{jj}, 'Units', 'Normal', 'Position', titlePos, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'left', ...
            'FontSize', fontSizeLabel);        
    shiftXTicks(xshiftTick);    
    yLim = get(gca, 'YLim');
    
    % this sets the ylimit, 
    ylimMax = 500;
    if(yLim(2) < ylimMax)
        set(gca, 'Ylim', [0, ylimMax]);
    else
        set(gca, 'Ylim', [0, yLim(2)*1.3]);
    end
    
    yTicks = get(gca, 'YTick');
    yTickLabels = sprintfc('%.0f', yTicks);
    set(gca, 'YTickLabel', yTickLabels);
    
    if exist('thresholds', 'var') == 1
        currentPost = get(gca, 'position');
        axes('position', currentPost);    
        boxplot(thresholds(:, jj), 'Orientation', 'horizontal', 'colors', 'k');
        set(gca, 'Ylim', get(gca', 'Ylim') - boxPlotYoffset)
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', {});
        set(gca, 'Xlim', xlim);
        set(gca, 'Color', 'none');
        set(gca, 'position', currentPost);
    end
end

