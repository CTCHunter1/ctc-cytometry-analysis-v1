function [ varargout ] = cloudPlot(X, Y, axisLimits, useLogScale, bins, cmaps, clims)
%CLOUDPLOT Does a cloud plot of the data in X and Y.
% CLOUDPLOT(X,Y) draws a cloudplot of Y versus X. A cloudplot is in essence
% a 2 dimensional histogram showing the density distribution of the data
% described by X and Y. As the plot displays density information, the
% dimensionality of X and Y are ignored. The only requirement is that X and
% Y have the same number of elements. Cloudplot is written to visualize
% large quantities of data and is most appropriate for datasets of 10000
% elements or more.
%
% CLOUDPLOT is fully compatible with SUBPLOT, COLORBAR and COLORMAP, and
% hold states.
%
% CLOUTPLOT(X,Y,axisLimits) plots the values of X and Y that are within the
% limits specified by axisLimits. axisLimits is a 4-element vector with the
% same meaning as "axis( axisLimits )" would have on a regular plot.
% axisLimits can be [] in which case it will be ignored.
%
% CLOUTPLOT(X,Y,axisLimits,useLogScale) If useLogScale == true, plots the
% base10 logarithm of the density at each point. Useful for distributions
% where the density varies greatly. Setting it empty will be itnerprated as
% linaer scale.
% 
% CLOUDPLOT(X,Y,axisLimits,useLogScale,bins) By setting bins to a 2-element
% vector, the number of bins used in the x and y direction can be adjusted.
% If you for example have a smaller data set, setting a lower number of
% bins may be clearer. The default value is to have one bin per pixel of
% the screen
%
% h = CLOUTPLOT(...) returns the handle to the cloud image in h.
%
% Example:
%   cloudPlot( randn(10000,1), randn(10000,1) );
%    Will plot a Gaussian scattered distribution.
%
%   cloudPlot( randn(1000), randn(1000), [0 3 -1 4]);
%    Will plot a Gaussian scattered distribution of X and Y values within
%    the limits.
%    
% Contributions:
% Thanks to Eric for the addition of the xscale/yscale output option.
% Thanks to Romesh who provided a bugfix, the bins input option, and the
% ability to resize the plot and a fix to remove the need to return scaling
% factors.

% GLF - Modified to allow for cell imput containing various lengths arrays
%       each cell will be ploted in a different colormap 
if nargin == 0
    CreateImage3PanelComboFigure;
    return;
end

% double array input
if ~iscell(X)
    Xp = X;
    Yp = Y;
    clear X Y
    X{1} = Xp;
    Y{1} = Yp;
end
    
% Check the data size
assert ( numel(X) == numel(Y), ...
    'The number of arrays in X and Y must be the same.' );

for ii = 1:numel(X)    
    assert ( numel(X{ii}) == numel(Y{ii}), ...
        'The number of elements in X{ii} and Y{ii} must be the same.' );
end


if nargin < 6 || isempty(cmaps)
    % Could make this more exciting
    % just make all colormaps the same
    for ii = 1:numel(X)
        cmaps{ii} = colormap;
    end
else
    if ~iscell(cmaps)
    cmapp = cmaps;
    clear cmaps
    cmaps{1} = cmapp;        
    end
    
    if(numel(X) > 1)
    assert ( numel(cmaps) == numel(X), ...
        'The number of colormaps must be the same as number of arrays X.');
    end
end


if nargin < 7 || isempty(clims)
    % Could make this more exciting
    % clims is undefined
    clims = {};
else
    if iscell(clims)
        assert ( numel(clims) == numel(X), ...
            'Passed as a cell the # of clims must be the same as number of arrays X.');
    else
        climsp = clims;
        clear clims;
        for ii = 1:numel(X)
            clims{ii} = climsp;
        end
    end
end

% If there is axis limits discard the points outside of the axis limits
if ( nargin >= 3 && ~isempty(axisLimits) )      
    for ii = 1:numel(X)
        pointSelect = X{ii}<=axisLimits(2) & X{ii}>=axisLimits(1) & ...
            Y{ii}<=axisLimits(4) & Y{ii}>=axisLimits(3);
        X{ii} = X{ii}(pointSelect);
        Y{ii} = Y{ii}(pointSelect);
    end
    axisLimitsSet = true;
else
    axisLimitsSet = false;
end

if ( nargin < 4 || isempty(useLogScale) )
    useLogScale = false;
end

if ( nargin < 5 )
    bins = [];
end


  
%Remove any nans or Infs in the data as they have no meaning in this
%context.
for ii = 1:numel(X)
    pointSelect = ~(isinf(X{ii}) | isnan(X{ii}) | isinf(Y{ii}) | isnan(Y{ii}));
    X{ii} = X{ii}(pointSelect);
    Y{ii} = Y{ii}(pointSelect);
end    

% Plot to get appropriate limits
h = [];
if ( axisLimitsSet )    
    for ii = 1:numel(X)
        h = plot( X{ii}(:), Y{ii}(:), '.' );        
        g = get( h, 'Parent' );    
        set (g, 'Xlim', [axisLimits(1) axisLimits(2)] );
        set (g, 'Ylim', [axisLimits(3) axisLimits(4)] );
    end
else
    for ii = 1:numel(X)
        h = plot( X{ii}(:), Y{ii}(:), '.' );        
        g = get( h, 'Parent' );  
        hold on;        
    end
end
xLim = get(g, 'Xlim' );
yLim = get(g, 'Ylim' );

%Get the bin size.
unitType = get(g,'Unit');
set(g,'Unit','Pixels')
axesPos = get(g,'Position');
nHorizontalBins = axesPos(3);
nVerticalBins = axesPos(4);
set(g,'Unit', unitType );

% Clear the data, as we actually don't want to see it.
if ( ~isempty(h) )
    set ( h, 'XData', [] );
    set ( h, 'YData', [] );
end

% Allocate an area to draw on
if ( isempty(bins) )
    bins = ceil([nHorizontalBins nVerticalBins ]);
else
    assert ( isnumeric(bins) && isreal(bins) );
    assert ( numel(bins) == 2 );
end
binSize(2) = diff(yLim)./(bins(2));
binSize(1) = diff(xLim)./(bins(1));

rgbImageNeg = zeros([bins, 3]);

% Draw in the canvas
for ii = 1:numel(X)
    canvas = zeros(bins);
    xBinIndex = floor((X{ii} - xLim(1))/binSize(1))+1;
    yBinIndex = floor((Y{ii} - yLim(1))/binSize(2))+1;

    % Added security: Make sure indexes are never outside canvas. May not be
    % possible. For the axislimits set case was taken care of
    pointsSelect = xBinIndex > 0 & xBinIndex <= bins(1) & ...
        yBinIndex > 0 & yBinIndex <= bins(2);
    xBinIndex = xBinIndex(pointsSelect);
    yBinIndex = yBinIndex(pointsSelect);
    
    % now we prepare multiple canvases, changed from previous
    for i = 1:numel(xBinIndex);
        canvas(xBinIndex(i),yBinIndex(i)) = ...
            canvas(xBinIndex(i),yBinIndex(i)) + 1;
    end     
    % to merege the colormaps we will store the negative image then 
    % uninvert it
    if(isempty(clims))
        clim = [min(min(canvas)), max(max(canvas))];        
    else
        clim = clims{ii};
    end
    
    %rescale the canvas to clims - remember color indexing starts at 1
    CmapLen = size(cmaps{ii},1);
    canvas = round((canvas/(clim(2) - clim(1)))*(CmapLen -1)) + 1;
    
    canvas(canvas < 1) = 1;
    canvas(canvas > CmapLen) = CmapLen;
        
    % do the transpose here
    if (useLogScale)
        rgbImageNeg = rgbImageNeg + (1 - ind2rgb(log(canvas.'),cmaps{ii}));
    else
        rgbImageNeg = rgbImageNeg + (1 - ind2rgb(canvas.',cmaps{ii}));    
    end
end

rgbImage = 1 - rgbImageNeg;
rgbImage(rgbImage > 1) = 1;
rgbImage(rgbImage < 0 ) = 0; 

% if there is only one color we can plot the old way. This will preserve
% existing interface functionality that assumes cData is NxNx1. 
if(numel(X) == 1)
    if ( useLogScale )
        h = imagesc(xLim, yLim, log10(canvas)');    
    else
        h = imagesc(xLim, yLim, canvas'); 
    end
else
    % Show the canvas and adjust the grids.
    h = image(xLim, yLim, rgbImage);
end

axis ( 'xy' );
axis ( 'tight' );
set ( g, 'units', 'normalized' ); % Enable resizing 

% Optionally return a handle.
if ( nargout == 1)
    varargout{1} = h;
end



