function [fp_avg, sens_avg] = averageSensSpec(fp, sens, theta)
% fp and sens are assumed to be jagged cell arrays possibly jagged
% theta is the angle of rotation prior to averaging, in radians

% create the rotation matrix
Mrot = [cos(theta), sin(theta); -sin(theta), cos(theta)];
MrotInv = inv(Mrot);

% fp and sens are 1xN length objects
N = length(fp);

udata = zeros(1, N);
vdata = zeros(1, N);
   
for jj = 1:N
    wdata = Mrot*[fp(jj); sens(jj)];
    udata(jj) = wdata(1);
    vdata(jj) = wdata(2);
end


% average the v_data and rotate it back
udata_avg = mean(udata);
vdata_avg = mean(vdata);

wp = MrotInv*[udata_avg; vdata_avg];
fp_avg = wp(1);
sens_avg = wp(2); 
