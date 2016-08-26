function [] = writePC1Eqn(fileName, eqnName, varNames, units, mu, sigma, coeff)

fid = fopen(fileName, 'w');

s = cell(1, length(coeff));
for ii = 1:length(coeff)
    if coeff(ii) >= 0
        s{ii} = '+';
    else
        s{ii} = '-';
    end
end

N = length(mu);

[~, idx] = sort(abs(coeff), 'descend');


for ii = 1:N
    sort_ind = idx(ii);
    % change the units if needed
    % print these units to strings here
    if strcmp(units{sort_ind}, ' \mu m ') 
        mu(sort_ind) = mu(sort_ind)*10^6;
        mu_str = sprintf('%.2f', mu(sort_ind));
        sigma(sort_ind) = sigma(sort_ind)*10^6;
        sigma_str = sprintf('%.2f', sigma(sort_ind));
    elseif strcmp(units{sort_ind}, ' \mu m^{-1} ') 
        mu(sort_ind) = mu(sort_ind)*10^-6;
        mu_str = sprintf('%.2f', mu(sort_ind));
        sigma(sort_ind) = sigma(sort_ind)*10^-6;
        sigma_str = sprintf('%.2f', sigma(sort_ind));
    elseif strcmp(units{sort_ind}, ' Counts ') 
        muLog = log10(mu(sort_ind));
        muLog = floor(muLog);
        mu_str = sprintf('%.2f \\rm{x}10^{%d}', mu(sort_ind)/10^muLog, muLog);        
        sigma(sort_ind)
        sigmaLog = log10(mu(sort_ind));
        sigmaLog = floor(muLog);
        sigma_str = sprintf('%.2f \\rm{x}10^{%d}', sigma(sort_ind)/10^sigmaLog, sigmaLog);        
    else
        mu_str = sprintf('%.2f', mu(sort_ind));
        sigma_str = sprintf('%.2f', sigma(sort_ind));
    end
    
       
    % write the eqn. tex file. 
    if ii == 1
    fprintf(fid, '%s \\; = \\; & %.2f * \\frac{%s - %s\\; %s}{%s \\; %s} \\\\', ...
        eqnName, coeff(sort_ind), varNames{sort_ind}, mu_str, units{sort_ind}, sigma_str, units{sort_ind});
    elseif ii < N
    fprintf(fid, '%s & \\; %.2f * \\frac{%s - %s \\; %s}{%s \\; %s } \\\\', ...
        s{sort_ind}, abs(coeff(sort_ind)), varNames{sort_ind}, mu_str, units{sort_ind}, sigma_str, units{sort_ind});
    else
    fprintf(fid, '%s & \\; %.2f * \\frac{%s - %s \\; %s }{%s \\; %s }', ...
        s{sort_ind}, abs(coeff(sort_ind)), varNames{sort_ind}, mu_str, units{sort_ind}, sigma_str, units{sort_ind});
    end    
end

fclose(fid);