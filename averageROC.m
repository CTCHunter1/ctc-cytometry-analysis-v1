function [fp_avg, sens_avg] = averageROC(fp, sens, theta, Ninterp)
% fp and sens are assumed to be jagged cell arrays possibly jagged
% theta is the angle of rotation prior to averaging, in radians

% we will average on to the longest of the curves. start by figuring out 
% the data ranges
if length(fp) ~= length(sens)
    return;
end

%remove values where fp,sens=0,1, thes are +/- INF on a zscored scale
for ii = 1:length(fp)
   Nfp = length(fp{ii}); % fp and sens are same length
%    bool = fp{ii} == 0 | fp{ii} ==1 | sens{ii} == 0 | sens{ii} ==1;
%    fp{ii} = (fp{ii}(~bool));
%    sens{ii} = (sens{ii}(~bool));    
   % instead change fp,sens=0,1 to be at 1/N
   fp{ii}(fp{ii}==0) = 1/Nfp;
   sens{ii}(sens{ii}==1) = 1-1/Nfp;
   %fp{ii}(fp{ii}==1) = 1-1/Nfp;
   %sens{ii}(sens{ii}==0) = 1/Nfp;
   % for these cases remove the data
   bool = fp{ii} == 1 | sens{ii} == 0;
   fp{ii} = fp{ii}(~bool);
   sens{ii} = sens{ii}(~bool);
end

Naverage = length(fp);

fp_min = min(fp{1});
fp_max = max(fp{1});
sens_min = min(sens{1});
sens_max = max(sens{1});

for ii = 2:Naverage
    fp_min_p = min(fp{ii});
    fp_max_p = max(fp{ii});
    sens_min_p = min(sens{ii});
    sens_max_p = max(sens{ii});    
    
    if fp_min_p < fp_min
        fp_min = fp_min_p;
    end
    
    if fp_max_p > fp_max;
        fp_max = fp_max_p;
    end
    
    if sens_min_p < sens_min
        sens_min = sens_min_p;
    end
    
    if sens_max_p > sens_max;
        sens_max = sens_max_p;
    end    
end

% we are going to do the interpolation in the rotated cordinate system

% create the rotation matrix
Mrot = [cos(theta), sin(theta); -sin(theta), cos(theta)];
MrotInv = inv(Mrot);

% make z-linear spaces on the two axis
fp_space = normcdf(linspace(norminv(fp_min), norminv(fp_max), Ninterp));
sens_space = normcdf(linspace(norminv(sens_min), norminv(sens_max), Ninterp));
v_space = zeros(1, Ninterp);
u_interp = zeros(1, Ninterp);

% create the rotated interpolated space
for ii = 1:Ninterp
    w = Mrot*[fp_space(ii); sens_space(ii)];
    u_interp(ii) = w(1);
    v_space(ii) = w(2);
end


% w_min = Mrot*norminv([fp_min; sens_min]);
% w_max = Mrot*norminv([fp_max; sens_max]);
% 
% % axis to interpolate onto equal spacing in z
% u_interp = linspace((w_min(1)), (w_max(1)), Ninterp);
%u_interp = (u_interp_z);

udata = cell(1, Naverage);
vdata = cell(1, Naverage);

%udata_interp = repmat(u_interp, Naverage, 1);
vdata_interp = zeros(Naverage, Ninterp);

% project data into rotated coordinate space
for ii = 1:Naverage
    udata{ii} = zeros(1, length(fp{ii}));
    vdata{ii} = zeros(1, length(fp{ii}));
    
    
    for jj = 1:length(fp{ii})
        wdata = Mrot*[fp{ii}(jj); sens{ii}(jj)];
        udata{ii}(jj) = wdata(1);
        vdata{ii}(jj) = wdata(2);
    end
    [udata{ii}, index] = unique(udata{ii});
    
    % interpolate on to this axis for averaging fill in extrapolated 
    % values with zero
    vdata_interp(ii, :) = interp1(udata{ii}, vdata{ii}(index), u_interp, 'linear',0);  
end

% average the v_data and rotate it back
vdata_interp_avg = mean(vdata_interp);
fp_avg = zeros(1, Ninterp);
sens_avg = zeros(1, Ninterp);

for ii = 1:Ninterp
    wp = MrotInv*[u_interp(ii); vdata_interp_avg(ii)];
    fp_avg(ii) = (wp(1));
    sens_avg(ii) =(wp(2)); 
end
