function cmap = getcolormap(n, varargin)
%GETMUCOLORMAP Summary of this function goes here
%   Detailed explanation goes here
P = inputParser;
addRequired(P, 'nUnit', @isnumeric)
addParameter(P,'mapFcn','brewermap',@(x) ischar(x) && ismember(x,{'brewermap','colorcet'}))
addParameter(P,'mapName','spectral',@ischar)
addParameter(P,'brightness',-0.1,@isscalar)
parse(P, n, varargin{:})

% get map
switch P.Results.mapFcn
    case 'brewermap'
        cmap = brewermap(n,P.Results.mapName);
    case 'colorcet'
        cmap = colorcet(P.Results.mapName,'N',n);
end

% adjust brightness
cmap = max([0 0 0],cmap+P.Results.brightness);

end

