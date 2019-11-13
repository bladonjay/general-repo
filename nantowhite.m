function [im, clim] = nantowhite(cvals, clim, cmap)

% Get the default color map.
if(nargin < 3); cmap = get(0,'DefaultFigureColormap');
elseif(ischar(cmap)); cmap = feval(cmap,64);
end
ncolors = size(cmap,1);

% Determine the range of the color values provided.
finites = isfinite(cvals(:));
if(~any(finites))
    cvalmin = 0;
    cvalmax = 0;
else
    cvalmin = min(cvals(finites));
    cvalmax = max(cvals(finites));
end

% Process the CLim input, if one is provided.
% If clim isn't provided, default to [-inf inf], which will will be set
% later to match the range of cavls.
if(nargin < 2 || isempty(clim)); clim = [-inf inf];
    
% If clim is a single value, treat it as a maximum (which cannot be -inf).
elseif(isscalar(clim) && clim~=-inf); clim = [-inf clim];
elseif(isscalar(clim));
    error('nantowhite:InvalidCLim','CLim(2) cannot be -inf.');
    
% Otherwise, make sure it is a two value vector.
elseif(numel(clim)~=2);
    error('nantowhite:InvalidCLim','CLim must be a scalar or 2 element vector.');

% Special case, the clim(1) and clim(2) are equal to each other, and equal
% to all values in cvals.
elseif(clim(1)==clim(2) && clim(1)==cvalmin && cvalmin==cvalmax); clim = [-inf inf];
    
% Otherwise, clim(2) must be greater than clim(1).
elseif(any(isnan(clim)) || clim(2)<=clim(1));
    error('nantowhite:InvalidCLim','CLim values must be increasing and non-NaN.');
end

% If clim(1) is -inf, replace with actual minimum.
% Three possibilities:
%   clim(2) == inf, we want clim(1) = cvalmin
%   clim(2) > cvalmin, we want clim(1) = cvalmin
%   clim(2) <= cvalmin, we want clim(1) < clim(2)
if(isinf(clim(1))); clim(1) = min([clim(2)-eps(clim(2))*64, cvalmin]); end

% If clim(2) is inf, replace with actual maximum.
% Three possibilities:
%   clim(1) == -inf, set above to be cvalmin, we want clim(2) > clim(1)
%   clim(1) < cvalmax, we want clim(2) = cvalmax
%   clim(1) >= cvalmax, we want clim(2) > clim(1)
if(isinf(clim(2))); clim(2) = max([clim(1)+eps(clim(1))*64, cvalmax]); end

% Distribute the bins equally among the CLim.
cbins = linspace(clim(1), clim(2), ncolors+1);

% Make the first and last bin '-inf' and 'inf' to include all values.
cbins(1) = -inf;
cbins(end) = inf;

% Bin the color values into the appropriate color bin.
[~,cvals] = histc(cvals,cbins);

% Move 'inf' values to the last bin.
cvals(cvals==ncolors+1) = ncolors;

% Initialize an entirely white image.
im = ones([size(cvals) 3]);

% Copy the appropriate values into the red, green, and blue layers.
for ii = 1:3;
    imlayer = im(:,:,ii);
    imlayer(cvals~=0) = cmap(cvals(cvals~=0),ii);
    im(:,:,ii) = imlayer;
end

end