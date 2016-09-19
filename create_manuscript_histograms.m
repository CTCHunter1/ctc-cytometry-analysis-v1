function [] = create_manuscript_histograms(Dp, Dn, Mixed, featuresFullData, DNACutOff, pathname)
% create_manuscript_histograms: 
% DP is disease poisitve = MCF7
% Dn is diesease negative = WBCs

% Copyright (C) 2016  Gregory L. Futia
% This work is licensed under a Creative Commons Attribution 4.0 International License.

if nargin == 0
   [Dp, Dn, Mixed, DNACutOff, pathname] = loadCytoData();   
end


channelNames = fieldnames(Dp.Channels);

% Emily's renaming of DAPI->DNA and Bodipy->Lipids
if(1)
    channelNamesLabels = strrep(channelNames, 'Bodipy', 'Lipids (Bodipy)');
    channelNamesLabels = strrep(channelNamesLabels, 'DAPI', 'DNA (DAPI)');
else
    channelNamesLabels = channelNames;
end


% Set the DNA Gates
DNApWBC = Dn.Channels.DAPI.totalSig_dBc > DNACutOff; % WBCs are diesease negative
DNApMCF7 = Dp.Channels.DAPI.totalSig_dBc > DNACutOff;% MCF7s are disease positive
DNApMix = Mixed.Channels.DAPI.totalSig_dBc > DNACutOff;

colorNames = {'g', 'y', 'b', 'r'};
reorder = [4, 2, 1, 3]; % order is from bottom to top

numBins = 100;

ao = .12;
aob = .0148;
aot = .0;
aw = .795;
ah = .0475;

size = [.25, 2, 6.5/2, 7.5];

set(gcf, 'Units', 'Inches');
set(gcf, 'Position', size);
set(gcf, 'PaperUnits', 'Inches');
set(gcf, 'PaperPosition', size);

set(0, 'defaultTextFontSize', 7);
set(0, 'defaultAxesFontSize', 7);
labelFontSize = 8;

boxPlotYoffset = .35;
titlePos = [.78, .45];
xshiftTick = .05;
linewidth = 1.0;  % default line width = .5

for kk = 1:length(channelNames)
    jj = reorder(kk); % reorder is defined above
    %jj = kk;
    
    axes('position', [ao, 4*aob+3*ah + (kk-1)*4*(aob+ah), aw, ah]);
    
    [n1, xout] = hist(Dp.Channels.(channelNames{jj}).totalSig_dBc(DNApMCF7), numBins);
    [n1, xout] = makejagged(n1, xout); 
    plot(xout, n1, 'r', 'LineWidth', linewidth); hold on;    
    [n2, xout2] = hist(Dn.Channels.(channelNames{jj}).totalSig_dBc(DNApWBC), numBins);
    [n2, xout2] = makejagged(n2, xout2); 
    n2 = max(n1)*n2/max(n2);
    h = plot(xout2, n2, 'b', 'LineWidth', linewidth); hold on; 
    
    p = ranksum(Dp.Channels.(channelNames{jj}).totalSig_dBc(DNApMCF7), Dn.Channels.(channelNames{jj}).totalSig_dBc(DNApWBC));
    fprintf('P=%e for %s total signal', p, channelNames{jj});
    if p == 0
        fprintf(' to machine prcession');
    end
    fprintf('\n');
    % clean up -INFs in CD45 dB channel
    %Mixed.Channels.(channelNames{jj}).totalSig_dBc(DNApMix)( ...
    %    Mixed.Channels.(channelNames{jj}).totalSig_dBc(DNApMix)== -INF) = -1;
    % -1 makes sense because log(1 count) = 0 so log(0) - can be right next
    % to it.   
    Mixed.Channels.(channelNames{jj}).totalSig_dBc( ...
        Mixed.Channels.(channelNames{jj}).totalSig_dBc == -Inf) = -1;
    [n3, xout3] = hist(Mixed.Channels.(channelNames{jj}).totalSig_dBc(DNApMix), numBins);
    [n3, xout3] = makejagged(n3, xout3); 
    n3 = max(n1)*n3/max(n3);
    plot(xout3, n3, 'k--', 'LineWidth', linewidth); hold on;    
      
    %set(gca, 'YTick', []);    
    yticks = get(gca, 'Ytick');
    set(gca, 'YTick', yticks(1:end-1));
    xlim = [56, 89];
    ylim = get(gca, 'YLim');
    set(gca, 'YLim', [0, ylim(2)*1.2]);
    set(gca, 'Xlim', xlim);    
    xlabel('dBct', 'Units', 'Normal', 'Position', [1.045, .0], 'FontSize', labelFontSize);    
    title([channelNamesLabels{jj}, ' \Sigma'], 'Units', 'Normal', 'Position', titlePos, ...
        'FontWeight', 'bold', 'FontSize', labelFontSize);
    shiftXTicks(xshiftTick);
    
    if exist('featuresFullData', 'var') == 1
        currentPost = get(gca, 'position');
        axes('position', currentPost);    
        ind = strcmp(featuresFullData.regressionNames, [channelNames{jj}, ' \Sigma']);
        boxplot(abs(featuresFullData.TMinDetect(:, ind)), 'Orientation', 'horizontal', 'colors', 'k');
        set(gca, 'Ylim', get(gca', 'Ylim') - boxPlotYoffset)
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', {});
        set(gca, 'Xlim', xlim);
        set(gca, 'Color', 'none');
        set(gca, 'position', currentPost);
    end
    
    % 2nd Moment
    axes('position', [ao, 3*aob+2*ah+ (kk-1)*4*(aob+ah), aw, ah]);
    [n1, xout] = hist(Dp.Channels.(channelNames{jj}).radius_m(DNApMCF7)*10^6, numBins);
    [n1, xout] = makejagged(n1, xout); 
    plot(xout, n1, 'r', 'LineWidth', linewidth); hold on;    
   
    [n2, xout2] = hist(Dn.Channels.(channelNames{jj}).radius_m(DNApWBC)*10^6, numBins);
    [n2, xout2] = makejagged(n2, xout2); 
    n2 = max(n1)*n2/max(n2);
    plot(xout2, n2, 'b', 'LineWidth', linewidth); hold on;    
        
    [n3, xout3] = hist(Mixed.Channels.(channelNames{jj}).radius_m(DNApMix)*10^6, numBins);
    [n3, xout3] = makejagged(n3, xout3); 
    n3 = max(n1)*n3/max(n3);
    plot(xout3, n3, 'k--', 'LineWidth', linewidth); hold on;
    
    p = ranksum(Dp.Channels.(channelNames{jj}).radius_m(DNApMCF7), Dn.Channels.(channelNames{jj}).radius_m(DNApWBC));
    fprintf('P=%e for %s radius', p, channelNames{jj});
    if p == 0
        fprintf(' to machine prcession');
    end
    fprintf('\n');
    
    yticks = get(gca, 'Ytick');
    set(gca, 'YTick', yticks(1:end-1));
    xlim = [.1, 9.9];
    set(gca, 'Xlim', xlim);   
    ylim = get(gca, 'YLim');
    set(gca, 'YLim', [0, ylim(2)*1.4]);   
    xlabel('\mum', 'Units', 'Normal', 'Position', [1.045, .06], 'FontSize', labelFontSize);
    title([channelNamesLabels{jj}, ' <r>'], 'Units', 'Normal', 'Position', titlePos, ...
        'FontWeight', 'bold', 'FontSize', labelFontSize);
    shiftXTicks(xshiftTick);
    %xlabel('Radius (m)');
    
    if exist('featuresFullData', 'var') == 1
        currentPost = get(gca, 'position');
        axes('position', currentPost);    
        ind = strcmp(featuresFullData.regressionNames, [channelNames{jj}, ' < r > ']);
        boxplot(abs(featuresFullData.TMinDetect(:, ind))*10^6, 'Orientation', 'horizontal', 'colors', 'k');
        set(gca, 'Ylim', get(gca', 'Ylim') - boxPlotYoffset)
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', {});
        set(gca, 'Xlim', xlim);
        set(gca, 'Color', 'none');
        set(gca, 'position', currentPost);
    end
    
    % Second spatial frequency moment    
    axes('position', [ao, aot+2*aob+ah+(kk-1)*4*(aob+ah), aw, ah]);
    [n1, xout] = hist(Dp.Channels.(channelNames{jj}).radius_invm(DNApMCF7)*10^-6, numBins);
    [n1, xout] = makejagged(n1, xout); 
    plot(xout, n1, 'r', 'LineWidth', linewidth); hold on;
    
    [n2, xout2] = hist(Dn.Channels.(channelNames{jj}).radius_invm(DNApWBC)*10^-6, numBins);
    [n2, xout2] = makejagged(n2, xout2); 
    n2 = max(n1)*n2/max(n2);
    plot(xout2, n2, 'b', 'LineWidth', linewidth); hold on;    
        
    [n3, xout3] = hist(Mixed.Channels.(channelNames{jj}).radius_invm(DNApMix)*10^-6, numBins);
    [n3, xout3] = makejagged(n3, xout3); 
    n3 = max(n1)*n3/max(n3);
    plot(xout3, n3, 'k--', 'LineWidth', linewidth); hold on;
    
    p = ranksum(Dp.Channels.(channelNames{jj}).radius_invm(DNApMCF7), Dn.Channels.(channelNames{jj}).radius_invm(DNApWBC));
    fprintf('P=%e for %s spat. freq. radius \n', p, channelNames{jj});
    if p == 0
        fprintf(' to machine prcession');
    end
    fprintf('\n');
    
    %set(gca, 'YTick', []);
    yticks = get(gca, 'Ytick');
    set(gca, 'YTick', yticks(1:end-1));
    xlim = [.66, 1.12];
    set(gca, 'Xlim', xlim);
    ylim = get(gca, 'YLim');
    set(gca, 'YLim', [0, ylim(2)*1.4]);   
    
    xlabel('\mum^{-1}', 'Units', 'Normal', 'Position', [1.060, .06], 'FontSize', labelFontSize);
    title([channelNamesLabels{jj}, ' <r_f>'], 'Units', 'Normal', ...
        'Position', [titlePos(1), titlePos(2)-.12], 'FontWeight', 'bold', ...
        'FontSize', labelFontSize);
    %xlabel('Spatial Freq. Radius (m^{-1})');
    shiftXTicks(xshiftTick);
    
    if exist('featuresFullData', 'var') == 1
        currentPost = get(gca, 'position');
        axes('position', currentPost);    
        ind = strcmp(featuresFullData.regressionNames, [channelNames{jj}, ' < r_f > ']);
        boxplot(abs(featuresFullData.TMinDetect(:, ind))*10^-6, 'Orientation', 'horizontal', 'colors', 'k');
        set(gca, 'Ylim', get(gca', 'Ylim') - boxPlotYoffset)
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', {});
        set(gca, 'Xlim', xlim);
        set(gca, 'Color', 'none');
        set(gca, 'position', currentPost);
    end
    
    % M2 product
    axes('position', [ao, aob + (kk-1)*4*(aob+ah), aw, ah]);
    [n1, xout] = hist(Dp.Channels.(channelNames{jj}).M2(DNApMCF7), numBins);
    [n1, xout] = makejagged(n1, xout); 
    plot(xout, n1, 'r', 'LineWidth', linewidth); hold on;
    
    [n2, xout2] = hist(Dn.Channels.(channelNames{jj}).M2(DNApWBC), numBins);
    [n2, xout2] = makejagged(n2, xout2); 
    n2 = max(n1)*n2/max(n2);
    plot(xout2, n2, 'b', 'LineWidth', linewidth); hold on;    
        
    [n3, xout3] = hist(Mixed.Channels.(channelNames{jj}).M2(DNApMix), numBins);
    [n3, xout3] = makejagged(n3, xout3); 
    n3 = max(n1)*n3/max(n3);    
    plot(xout3, n3, 'k--', 'LineWidth', linewidth); hold on;
    
    p = ranksum(Dp.Channels.(channelNames{jj}).M2(DNApMCF7), Dn.Channels.(channelNames{jj}).M2(DNApWBC));
    fprintf('P=%e for %s <M> \n', p, channelNames{jj});
    if p == 0
        fprintf(' to machine prcession');
    end
    fprintf('\n');
    
    %set(gca, 'Ytick', []);
    yticks = get(gca, 'Ytick');
    set(gca, 'YTick', yticks(1:end-1));
    xlim = [1.2, 8];
    set(gca, 'Xlim', xlim);
    ylim = get(gca, 'YLim');
    set(gca, 'YLim', [0, ylim(2)*1.4]);   
    
    title([channelNamesLabels{jj}, '<M>'], 'Units', 'Normal', 'Position', titlePos, ...
    'FontWeight', 'bold', 'FontSize', labelFontSize);
    %xlabel('M2');   
    shiftXTicks(xshiftTick);
    
     if exist('featuresFullData', 'var') == 1
        currentPost = get(gca, 'position');
        axes('position', currentPost);    
        ind = strcmp(featuresFullData.regressionNames, [channelNames{jj}, ' < M > ']);
        boxplot(abs(featuresFullData.TMinDetect(:, ind)), 'Orientation', 'horizontal', 'colors', 'k');
        set(gca, 'Ylim', get(gca', 'Ylim') - boxPlotYoffset)
        set(gca, 'YTick', []);
        set(gca, 'XTickLabel', {});
        set(gca, 'Xlim', xlim);
        set(gca, 'Color', 'none');
        set(gca, 'position', currentPost);
    end
    
end

h = text(0, 0, 'Frequency', 'rotation', 90, 'FontSize', labelFontSize);
set(h, 'Units', 'Normal');
set(h, 'Position', [-.13, -6, 0]);
%h = text(-.05, 5.8, filename, 'Units', 'Normal', 'Interpreter', 'None');
%h2 = text(.25, 6.0, date_str, 'Units', 'Normal', 'Interpreter', 'None');


set(gcf, 'PaperUnits', 'Inches');
set(gcf, 'PaperPosition', size); 

saveas(gcf, [pathname, 'FeatureHistograms', '.emf'], 'emf'); 
saveas(gcf, [pathname, 'FeatureHistograms', '.jpeg'], 'jpg'); 
saveas(gcf, [pathname, 'FeatureHistograms', '.png'], 'jpg'); 
saveas(gcf, [pathname, 'FeatureHistograms', '.eps'], 'epsc'); 

end

% function shiftXTicks(dy)
% xtickpos = get(gca, 'Xtick');
% ylimvals = get(gca, 'Ylim');
% set(gca, 'XTickLabel', []);
% 
%     for ii = 1:length(xtickpos);
%     t = text(xtickpos(ii), ylimvals(1), num2str(xtickpos(ii)));
%     set(t, 'Units', 'Normal');
%     set(t, 'Position', get(t, 'Position') + [0, dy, 0] -[0, .18, 0] );
%     %set(t, 'Units', 'pixels');
%     t1 = get (t, 'Position');    
%     t1ex = get(t, 'Extent');
%     set(t, 'Position', [t1(1) - t1ex(3)/4, t1(2), t1(3)]);
%     end
%     
% end
