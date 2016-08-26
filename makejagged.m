function [nout, xout] = makejagged(nin, xin)
    
    N = length(xin);
    
    xout = zeros(1, N*2-1);
    nout = zeros(1, N*2-1);
    
    xout(1:2:N*2)= xin;
    xout(2:2:N*2-1) = xin(2:end);
    
    nout(1:2:N*2) = nin;
    nout(2:2:N*2-1) = nin(1:end-1);
end
