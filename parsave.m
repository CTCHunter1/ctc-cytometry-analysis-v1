function parsave(filename, varargin)
narginchk(2, Inf);

nargoutchk(0, 0);

for I = 2:nargin

    varname = genvarname(inputname(I));
    eval([varname ' = varargin{' num2str(I-1) '};'])

    if (I == 2)
        save(filename, varname)
    else
        save(filename, varname, '-append')
    end

end
