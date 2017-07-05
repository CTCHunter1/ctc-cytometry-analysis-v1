function [fp_avg, sens_avg] = averageSensSpec3(fp, sens)
% fp and sens are assumed to be jagged cell arrays possibly jagged
% theta is the angle of rotation prior to averaging, in radians

LRtp = sens./fp;
LRtn = (1-sens)./(1-fp);

LRtp_avg = mean(LRtp);
LRtn_avg = mean(LRtn);

fp_avg = (1-LRtn_avg)./(LRTp_avg-1);
sens_avg = (1-LRtn_avg)./(1-(LRtn_avg./LRtp_avg));

 
