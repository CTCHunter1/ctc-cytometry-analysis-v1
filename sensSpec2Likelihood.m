function [LRTp, LRTn] = sensSpec2Likelihood(sens, spec, bCleanUp)
% likelihood ratio definitions
% LR+ = sens/(1-spec) - Undifined if FP = 0, also will limit
% sens > 0
% LR- = (1-sens)./spec - Undefined if SENS = 1, also will limit
% 1-sens > 0

% cleanup flag excludes the INFs for those conditions
LRTp = sens./(1-spec);
LRTn = (1-sens)./spec;

% Remove INF/NAN/ and 0s
if exist('bCleanUp', 'var')
    if bCleanUp == 1
        indRem = LRTp == Inf | isnan(LRTp) | LRTp == 0 ...
            | LRTn == Inf | isnan(LRTn) | LRTn == 0;
        % a useless biomarker has LRTp=LRTn=1
        LRTp(indRem) = 1;
        LRTn(indRem) = 1;
    end
end
   
