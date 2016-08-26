% returns sensitivity and fp estimate for D+, D- inputs
function [fp, sens, T, AUC] = compute_ROC(Dp, Dn, method, bGenerateFigures, pathname)

if nargin == 0
    create_manuscript_figures; 
    return;    
elseif nargin <= 2
    method = 'perfcurve';
elseif nargin == 3
    bGenerateFigures = 0;    
end

uDp = mean(Dp);
sigDp = std(Dp);
uDn = mean(Dn);
sigDn = std(Dn);

bNeg = 0; % code below not ran
% ROC curves typically assume the mean of the Dp is greater than
% the mean of Dn if this is not true flip them so we get normal looking 
% ROC curves
if uDn > uDp
    % old code is swap them
    %temp = Dp;
    %Dp = Dn;
    %Dn = temp; 
    % new code is take the negative
    Dp = - Dp;
    Dn = - Dn;
    bNeg = 1;
end

uDp = mean(Dp);
sigDp = std(Dp);
uDn = mean(Dn);
sigDn = std(Dn);


Y = zeros(1, length(Dp) + length(Dn), 'uint8');
Y(1:length(Dp)) = 1; % set everything at the end equal to 1

X = [Dp, Dn];

% perform using perfcurve
switch method
    case {'GregsFitLorentz', 'GregsFitGaussLin', 'GregsFitGaussLog'}
        Xmin = min(X)*1.5;
        Xmax = max(X)*1.5;
        numBins = 50;
        xCuts = linspace(Xmin, Xmax, numBins);
        
        
        %DpBin = xCuts;
        %DnBin = xCuts;
        
        [nDp, DpBin] = hist(Dp, numBins);
        % window nDp
        cut = -.19;
        nDp = nDp(DpBin > cut);
        DpBin = DpBin(DpBin > cut);
        
        cut = -cut;
        [nDn, DnBin] = hist(Dn, numBins);
        nDn = nDn(DnBin < cut);
        DnBin = DnBin(DnBin < cut);
        
        lorFun1 = @(P, xdata) P(1)./(abs((xdata - P(2))./P(3)).^2+ 1).^1.8 ;
        lorFun2 = @(P, xdata) P(1)./(abs((xdata - P(2))./P(3)).^2 + 1).^1.8 ;
        
        %gaussfun1 = @(P, xdata) P(1).*exp(-abs((xdata-P(2))./P(3)).^2);
        
        [nDpMax, IMax] = max(nDp);
        x0Dp = [nDpMax, uDp, sigDp];
        parmDp = lsqcurvefit(lorFun1, x0Dp, DpBin, nDp);       
      
        [nDnMax, IMax] = max(nDn);
        x0Dn = [nDnMax, uDn, sigDn];
        parmDn = lsqcurvefit(lorFun2, x0Dn, DnBin, nDn);       
      
        xInterp = linspace(Xmin, Xmax, numBins*4);        
        if strcmp(method, 'GregsFitGaussLin')
            fDp = fit(DpBin.', nDp.', 'gauss2');
            fDn = fit(DnBin.', nDn.', 'gauss2');
        
            nDp_fit = fDp(xInterp);
            nDn_fit = fDn(xInterp);       
        
        elseif strcmp(method, 'GregsFitGaussLog')        
            % remove nDp = 0 points
            DpBin = DpBin(nDp > 0);
            nDp = nDp(nDp > 0);

            DnBin = DnBin(nDn > 0);
            nDn = nDn(nDn > 0);

            % fit the gauss in the log domain, the gauss in a log
            % domain is a polynomial 
            fDp = fit(DpBin.', log10(nDp.'), 'poly2');
            %ParmPolyDp = polyfit(DpBin, log10(nDp), 2);
            fDn = fit(DnBin.', log10(nDn.'), 'poly2');
            %ParmPolyDn = polyfit(DnBin, log10(nDn), 2);
            nDp_fit = 10.^fDp(xInterp);
            nDp_fit = nDp_fit/max(nDp_fit)*nDpMax;

            nDn_fit = 10.^fDn(xInterp);
            nDn_fit = nDn_fit/max(nDn_fit)*nDnMax;         
        end
        
        
        
        %nDp_fit = 10.^polyval(ParmPolyDp, xInterp);        
        %nDn_fit = 10.^polyval(ParmPolyDn, xInterp); 
        
        if strcmp(method, 'GregsFitLorentz')
            nDp_fit = lorFun1(parmDp, xInterp);
            nDn_fit = lorFun2(parmDn, xInterp);
        end
        
        NExtrap = length(nDp_fit);
               
              
%        nDpEx = interp1(DpBin, nDp, xCuts, 'pchip', 'extrap');
%        nDnEx = interp1(DnBin, nDn, xCuts, 'pchip', 'extrap');
        
        NDpFit = sum(nDp_fit);
        NDnFit = sum(nDn_fit);
        
        if(bGenerateFigures)
            figure(7);
            clf;
            set(gcf, 'Units', 'Inches');
            set(gcf, 'Position', [4, 4, 2.5, 2.5]);
            set(gcf, 'PaperUnits', 'Inches');
            set(gcf, 'PaperPosition', [4, 4, 2.5, 2.5]);
            
            clf;
            set(0, 'defaultTextFontSize', 16);           
            set(0, 'defaultAxesFontSize', 16);
            
            ax1 = axes('position', [.18, .22, .78, .3]);
            plot(DpBin, log10(nDp), 'r'); hold on;
            plot(DnBin, log10(nDn), 'b');
            set(gca, 'Ylim', [0, 6]);
            set(gca, 'Xlim', [-1.5, 1.5]);
            ylabel('log_{10}(N)');
            xlabel('Regressed Cord.');     
            set(gca, 'YTick', []);
            
            %figure(8);
            plot(xInterp, log10(nDp_fit), 'r--'); hold on;
            plot(xInterp, log10(nDn_fit), 'b--');

            ax2 = axes('position', [.18, .65, .78, .3]);
            plot(DpBin, nDp, 'r'); hold on;
            plot(DnBin, nDn, 'b');
            
            %figure(8);
            plot(xInterp, (nDp_fit), 'r--'); hold on;
            plot(xInterp, (nDn_fit), 'b--');
            set(gca, 'Xlim', [-1.5, 1.5]);
            ylabel('N');
            set(gca, 'YTick', []);

        
            if nargin >= 5               
                saveas(gcf, [pathname, 'DistFit', '.png'], 'png'); 
                saveas(gcf, [pathname, 'DistFit', '.eps'], 'epsc');                 
            end
        end
        
        %set(gca, 'Ylim', [-3, 3]);
        sens = zeros(1, NExtrap);
        fp = zeros(1, NExtrap);
        
        for ii = 1:NExtrap
            sens(ii) = sum(nDp_fit(ii:end))./NDpFit;
            fp(ii) = sum(nDn_fit(ii:end))./NDnFit;            
        end
        
     case 'perfcurve'
        [fp, sens, T, AUC] = perfcurve(Y, X, 1);
     
     case  'GregsNormal' % crate histograms but don't fit/ extrapolate
        numBins = 200; % this should be input parameter
        Xmin = min(X);
        Xmax = max(X);
        xCuts = linspace(Xmin, Xmax, numBins);
        
        [nDp] = hist(Dp, xCuts);
        [nDn] = hist(Dn, xCuts);
         
        NDp = length(Dp);
        NDn = length(Dn);
        
        for ii = 1:numBins
            sens(ii) = sum(nDp(ii:end))./NDp;
            fp(ii) = sum(nDn(ii:end))./NDn;            
        end
        
        %[sens, fp] = makejagged(sens, fp);
end

if(bNeg == 1)
    %undo the negative thing we did
    T = -T;
end
