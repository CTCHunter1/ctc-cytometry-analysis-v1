function [fp_avg, sens_avg] = averageSensSpec2(fp, sens, theta)
% fp and sens are assumed to be jagged cell arrays possibly jagged
% theta is the angle of rotation prior to averaging, in radians

Npts = length(fp);

fpCell = cell(1, Npts);
sensCell = cell(1, Npts);

for ii = 1:Npts
   fpCell{ii} = fp(ii);
   sensCell{ii} = sens(ii);
end

[fp_avg, sens_avg] = averageROC(fpCell, sensCell, theta, Npts);