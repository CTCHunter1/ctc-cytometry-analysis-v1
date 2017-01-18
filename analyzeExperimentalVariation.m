% Copyright (C) 2017  Gregory L. Futia
% This work is licensed under a Creative Commons Attribution 4.0 International License.

% Description: 

mainPath = [pwd, filesep, 'Processed Mats', filesep, 'Futia_Zenodo_2016'];

daysToProcess = {'day5', 'day7', 'day8', 'day9', 'day10', 'day12', 'day13', ...
    'day14', 'day15'};

daysToProcess = {'day5', 'day7',};

daysToProcess = {'day9', 'day10', 'day12', 'day13', ...
    'day14', 'day15'};



% for each day there are two folders. One of pure MCF7 data and one of pure
% WBC data. For each pairing of this data, compute cohen's D, and and ROC
% curve. Find the positon on the ROC curve that maximizes Ndet
% for each of those cut off positions ask what fp and sens is in remaining
% testing data. 
fsepchar = filesep;

DNACutoff = 65;

% load the regressions form files
% otherwise compute them
bLoad = 1; % 0 for compute; 1 for load. Running code in compute can take days.
bUsePooledRegs = 0; % the all regression instead of the training-testing pairings
bSwapTrainTest = 0; % 0 small train big test (normal), 1 big train small test 

Ndays = length(daysToProcess);

% one cell for each day of data
WBCDataAll = cell(1, Ndays);
MCF7DataAll = cell(1, Ndays);
MixDataAll = cell(1, Ndays);

% This will interate through each day's worth of data
for ii=1:length(daysToProcess)
    fprintf('Processing %s\n', daysToProcess{ii});    
    % Figure out the numbe of WBC and MCF7 file in these days
    
    wbcPath = [mainPath, fsepchar , daysToProcess{ii}, fsepchar, 'WBC', ...
        fsepchar];
    dirResult = dir([wbcPath '*.mat']);    
    
    Nwbcfiles = length(dirResult);
    wbcData = cell(1, Nwbcfiles);
    wbcFilenames = cell(1, Nwbcfiles);
    
    for jj=1:Nwbcfiles
        wbcData{jj} = load([wbcPath dirResult(jj).name]);
        wbcData{jj} = wbcData{jj}.ND;
        wbcFilenames{jj} = dirResult(jj).name;
    end    
        
    mcf7Path = [mainPath, fsepchar , daysToProcess{ii}, fsepchar, 'MCF7', ...
        fsepchar];
    dirResult = dir([mcf7Path, '*.mat']);       
    Nmcf7files = length(dirResult);
    mcf7Filenames = cell(1, Nmcf7files);
    
    mcf7Data = cell(1, Nmcf7files); 
   
    for jj=1:Nmcf7files
        mcf7Data{jj} = load([mcf7Path dirResult(jj).name]);
        mcf7Data{jj} = mcf7Data{jj}.ND;
        mcf7Filenames{jj} = dirResult(jj).name;
    end
    
       
    mixPath = [mainPath, fsepchar , daysToProcess{ii}, fsepchar, 'Mix', ...
        fsepchar];
    dirResult = dir([mixPath, '*.mat']);
    Nmixfiles = length(dirResult);
    mixFilenames = cell(1, Nmixfiles);
    
    mixData = cell(1, Nmixfiles);
    
    for jj=1:Nmixfiles
        mixData{jj} = load([mixPath dirResult(jj).name]);
        mixData{jj} = mixData{jj}.ND;
        mixFilenames{jj} = dirResult(jj).name;
    end
    
    % put the technical varients in to full day structures
    WBCDataAll{ii} = wbcData;
    MCF7DataAll{ii} = mcf7Data;
    MixDataAll{ii} = mixData;
end


channelNames = fieldnames(WBCDataAll{1}{1}.Channels);
features = {'totalSig_dBc', 'radius_m', 'radius_invm', 'M2'};
numBins = 50;

figure (1);

numBins = 100;

ao = .12;
aob = .10;
aot = .0;
aw = .83;
ah = .11;
ahscale = .9;

figSize = [.25, 2, 6.4/2, 5.5];

titlePos = [.26, .50];

set(gcf, 'Units', 'Inches');
set(gcf, 'Position', figSize);
set(gcf, 'PaperUnits', 'Inches');
set(gcf, 'PaperPosition', figSize);

set(0, 'defaultTextFontSize', 8);
set(0, 'defaultAxesFontSize', 8);
fontSizeLabel = 7;
xshiftTick = .052;

numAxis = length(daysToProcess);
numHistPts = 50;

featureLabels = features;
featureUnits = cell(1,4);

figure(1);
for ii = 1:length(channelNames)
    for jj = 1:length(features)
        % change for different feature types
        scale = 1;
        xlim = [0, 20];
        
        % plot the Dp data
        switch features{jj}
            case 'M2'
                xlim = [0, 9.5];
                featureLabels{jj} = '<M2>';
                featureUnits{jj} = '';
            case 'totalSig_dBc'
                xlim = [60 84];
                featureLabels{jj} = '\Sigma';
                featureUnits{jj} = 'dBc';
            case 'radius_m'
                scale = 10^6;
                xlim = [1 9.5];                
                featureLabels{jj} = '<r>';
                featureUnits{jj} = '\mu m';
            case 'radius_invm'
                scale = 10^-6;
                xlim = [.65, 1.5] ;
                featureLabels{jj} = '<r_f>';
                featureUnits{jj} = '\mu m^{-1}';
        end       
        
        wbc_merge = cell(1, length(daysToProcess));
        wbc_merge_offset = cell(1, length(daysToProcess)); % renormalize to 
        % offset dtermined by MCF7
        Nwbc_arr = zeros(1,  length(daysToProcess)); % so we don't have to 
        % go back and recount
        anova_ss_wbc = zeros(1, length(daysToProcess));
        anova_df_wbc = zeros(1, length(daysToProcess));
        anova_bs_wbc = zeros(1, length(daysToProcess));
        anova_bsdf_wbc = zeros(1, length(daysToProcess));
                
        mcf7_merge = cell(1, length(daysToProcess));
        mcf7_merge_offset = cell(1, length(daysToProcess)); % renormalize to 
        % offset dtermined by WBC
        Nmcf7_arr = zeros(1,  length(daysToProcess)); % so we don't have to 
        % go back and recount
        anova_ss_mcf7 = zeros(1, length(daysToProcess));
        anova_df_mcf7 = zeros(1, length(daysToProcess));
        anova_bs_mcf7 = zeros(1, length(daysToProcess));
        anova_bsdf_mcf7 = zeros(1, length(daysToProcess));
        
        clf;        
        for kk = 1:length(daysToProcess)
            axes('position', [ao, aob+ ((kk-1)/numAxis)*ahscale, aw, ah]);   
            
            % we get to loop through twice because growing arrays in matlab
            % isn't efficent
            Ndaywbc = 0;
            n1 = 200;
            for ll = 1:length(WBCDataAll{kk})
                wbcDNDp = WBCDataAll{kk}{ll}.Channels.DAPI.totalSig_dBc > DNACutoff;
                Ndaywbc = Ndaywbc + sum(wbcDNDp);
                [n2, xout2] = hist(scale*WBCDataAll{kk}{ll}.Channels.(channelNames{ii}).(features{jj})(wbcDNDp), numHistPts);                
                [n2, xout2] = makejagged(n2, xout2); 
                n2 = n2/sum(wbcDNDp);
                plot(xout2, n2, 'b'); hold on;    
            end
            % allocate arrays for the full day merged
            wbcFeatMerg = zeros(1, Ndaywbc);
            wbcFeatLab = zeros(1, Ndaywbc);
            lastInd = 0;
            for ll = 1:length(WBCDataAll{kk})
                wbcDNDp = WBCDataAll{kk}{ll}.Channels.DAPI.totalSig_dBc > DNACutoff;
                Nptsin = length(WBCDataAll{kk}{ll}.Channels.(channelNames{ii}).(features{jj})(wbcDNDp));
                wbcFeatMerg(lastInd+[1:Nptsin]) = WBCDataAll{kk}{ll}.Channels.(channelNames{ii}).(features{jj})(wbcDNDp); % the data
                wbcFeatLab(lastInd+[1:Nptsin]) = ll;
                lastInd = lastInd + Nptsin;
            end
            
            
            
            Ndaymcf7 = 0;
            for ll = 1:length(MCF7DataAll{kk})
                mcf7DNDp = MCF7DataAll{kk}{ll}.Channels.DAPI.totalSig_dBc > DNACutoff;
                Ndaymcf7 = Ndaymcf7 + sum(mcf7DNDp);
                [n2, xout2] = hist(scale*MCF7DataAll{kk}{ll}.Channels.(channelNames{ii}).(features{jj})(mcf7DNDp), numHistPts);
                [n2, xout2] = makejagged(n2, xout2); 
                n2 = n2/sum(mcf7DNDp);
                plot(xout2, n2, 'r'); hold on;    
            end
            % allocate arrays for the full day merged
            mcf7FeatMerg = zeros(1, Ndaymcf7);
            mcf7FeatLab = zeros(1, Ndaymcf7);
            lastInd = 0;
            for ll = 1:length(MCF7DataAll{kk})
                mcf7DNDp = MCF7DataAll{kk}{ll}.Channels.DAPI.totalSig_dBc > DNACutoff;
                Nptsin = length(MCF7DataAll{kk}{ll}.Channels.(channelNames{ii}).(features{jj})(mcf7DNDp));
                mcf7FeatMerg(lastInd+[1:Nptsin]) = MCF7DataAll{kk}{ll}.Channels.(channelNames{ii}).(features{jj})(mcf7DNDp); % the data
                mcf7FeatLab(lastInd+[1:Nptsin]) = ll;
                lastInd = lastInd + Nptsin;
            end
            
            set(gca, 'Xlim', xlim);
                        
            [p, astats] = anova1(wbcFeatMerg, wbcFeatLab, 'off');
            tobj = text(2, 100, sprintf('P_{wbc} = %2.2f', p*100));
            set(tobj, 'Units', 'normal');
            set(tobj, 'Position', [.05, .8, 0]);
            set(tobj, 'FontSize', 7);
            
            anova_ss_wbc(kk) = astats{2,2};
            anova_df_wbc(kk) = astats{2,3};
            anova_bs_wbc(kk) = astats{3,2};
            anova_bsdf_wbc(kk) = astats{3,3};
            wbc_merge{kk} = wbcFeatMerg;      
            % reoffst the data against the mcf7 position
            if(kk > 1) 
                wbc_merge_offset{kk} = wbcFeatMerg + -mean(mcf7FeatMerg) + mean(mcf7_merge{1});
                % controlt to first sample, same data should be p=1
                % wbc_merge_offset{kk} = wbcFeatMerg + -mean(wbcFeatMerg) + mean(wbc_merge{1});
            else
                wbc_merge_offset{kk} = wbcFeatMerg;
            end
            Nwbc_arr(kk) = Ndaywbc;
                        
            [p, astats] = anova1(mcf7FeatMerg, mcf7FeatLab, 'off');
            tobj =text(10, 100, sprintf('P_{mcf7} = %2.2f', p*100));
            set(tobj, 'Units', 'normal');
            set(tobj, 'Position', [.7, .8, 0]);
            set(tobj, 'FontSize', 7);
            
            anova_ss_mcf7(kk) = astats{2,2};
            anova_df_mcf7(kk) = astats{2,3};
            anova_bs_mcf7(kk) = astats{3,2};
            anova_bsdf_mcf7(kk) = astats{3,3};
            mcf7_merge{kk} = mcf7FeatMerg;
             % reoffst the data against the wbc position
            if(kk > 1) 
                mcf7_merge_offset{kk} = mcf7FeatMerg + -mean(wbcFeatMerg) + mean(wbc_merge{1});
                % control to first sample, same data should be p=1
                % mcf7_merge_offset{kk} = mcf7FeatMerg + -mean(mcf7FeatMerg) + mean(mcf7_merge{1});
            else
                mcf7_merge_offset{kk} = mcf7FeatMerg;
            end
            Nmcf7_arr(kk) = Ndaymcf7;
            
            
            tobj =text(10, 100, daysToProcess{kk});
            set(tobj, 'Units', 'normal');
            set(tobj, 'Position', [.4, .9, 0]);
            
            tobj =text(10, 100, featureUnits{jj});
            set(tobj, 'Units', 'normal');
            set(tobj, 'Position', [.95, -.15, 0]);
            
        end
        
        wbc_days_merg = zeros(1, sum(Nwbc_arr));
        wbc_day_merg_offset = zeros(1, sum(Nwbc_arr));
        wbc_day_merg_lab = zeros(1, sum(Nwbc_arr));
        lastIndWbc = 0;
        mcf7_days_merg = zeros(1, sum(Nmcf7_arr));
        mcf7_day_merg_offset = zeros(1, sum(Nwbc_arr));
        mcf7_day_merg_lab = zeros(1, sum(Nmcf7_arr));
        lastIndMcf7 = 0;
        % merge the days data and perform final anova
        for kk=1:length(daysToProcess)
            wbc_days_merg(lastIndWbc + [1:Nwbc_arr(kk)]) = wbc_merge{kk};
            wbc_day_merg_offset(lastIndWbc + [1:Nwbc_arr(kk)]) = wbc_merge_offset{kk};
            wbc_day_merg_lab(lastIndWbc + [1:Nwbc_arr(kk)]) = kk;
            lastIndWbc = lastIndWbc + Nwbc_arr(kk); 
            mcf7_days_merg(lastIndMcf7 + [1:Nmcf7_arr(kk)]) = mcf7_merge{kk};
            mcf7_day_merg_offset(lastIndMcf7 + [1:Nmcf7_arr(kk)]) = mcf7_merge_offset{kk};
            mcf7_day_merg_lab(lastIndMcf7 + [1:Nmcf7_arr(kk)]) = kk;
            lastIndMcf7 = lastIndMcf7 + Nmcf7_arr(kk);            
        end
        
        [p_day_wbc, astats_wbc_day] = anova1(wbc_days_merg, wbc_day_merg_lab, 'off');
        [p_day_mcf7, astats_mcf7_day] = anova1(mcf7_days_merg, mcf7_day_merg_lab, 'off');
        
        [p_day_wbc_offset, astats_wbc_day_offset] = anova1(wbc_day_merg_offset, wbc_day_merg_lab, 'off');
        [p_day_mcf7_offset, astats_mcf7_day_offset] = anova1(mcf7_day_merg_offset, mcf7_day_merg_lab, 'off');
        
        % null hypothesis group error no different than that of individual
        % anovas
        % now wish to determine if the F number is different        
        MSdenom = sum(anova_ss_wbc)./sum(anova_df_wbc);
        BSdenom = sum(anova_bs_wbc)/sum(anova_bsdf_wbc);
        F = astats_wbc_day{2,4}./MSdenom;
        Foffset = astats_wbc_day_offset{2,4}./MSdenom;
        Ftechwbc = MSdenom/BSdenom;
        p_wbc = 1 - fcdf(F, astats_wbc_day{2,3}, sum(anova_df_wbc));
        p_wbc_offset = 1 - fcdf(Foffset, astats_wbc_day_offset{2,3}, sum(anova_df_wbc));
        
        a = annotation('textbox',[.01, .025, .98,.05],'String', ...
            sprintf('Day v. Tech WBC: P = %2.2f%%, P_{off} = %2.2f%%, f= %2.1f, f_{tech}=%2.1f', 100*p_wbc, 100*p_wbc_offset, F, Ftechwbc));
        a.EdgeColor = 'none';
        a.FontSize = 7;
        
        % now wish to determine if the F number is different
        MSdenom = sum(anova_ss_mcf7)./sum(anova_df_mcf7);
        BSdenom = sum(anova_bs_mcf7)/sum(anova_bsdf_mcf7);
        F = astats_wbc_day{2,4}./MSdenom;
        Foffset = astats_wbc_day_offset{2,4}./MSdenom;
        Ftechmcf7 = MSdenom/BSdenom;        
        p_mcf7 = fpdf(F, astats_mcf7_day{2,3}, sum(anova_df_mcf7));
        p_mcf7_offset = 1 - fcdf(Foffset, astats_mcf7_day_offset{2,3}, sum(anova_df_wbc));
        
        a = annotation('textbox', [.01, .002, .98,.05],'String', ...
            sprintf('Day v. Tech MCF7: P = %2.2f%%, P_{off} = %2.2f%%, f= %2.1f, f_{tech}=%2.1f', 100*p_mcf7, 100*p_mcf7_offset, F, Ftechmcf7));
        a.EdgeColor = 'none';
        a.FontSize = 7;
        
        a = annotation('textbox',[.15,.96,.4,.05],'String',[channelNames{ii}, ' ', featureLabels{jj}]);
        
        a.FontSize = 11;
        a.EdgeColor = 'none';
        
        saveas(gcf, [mainPath, filesep, 'VariationHist', filesep, channelNames{ii}, '_', features{jj} '.png'], 'png');
    end    
end
