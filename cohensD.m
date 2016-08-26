function [d] = cohensD(x1, x2)
numerator = (length(x1)-1)*var(x1) + (length(x2)-1)*var(x2);
denominator = length(x1) + length(x2) - 1;
s = sqrt(numerator/denominator);
d = (mean(x1) - mean(x2))/s;

end
