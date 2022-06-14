function plotforce(forces, conditions, experiment, sessionIndex, varargin)

P = inputParser;
addRequired(P, 'forces', @istable)
addRequired(P, 'conditions', @istable)
addRequired(P, 'experiment', @(x) ischar(x) && ismember(x, {'stim','dynamic','posture'}))
addRequired(P, 'sessionIndex', @isnumeric)
addOptional(P, 'layers', {'trials','mean','count'}, @(x) iscell(x) && all(ismember(x,{'trials','mean','sem','target','count'})))
% key selection
addParameter(P, 'conditionBlock', 1, @isscalar)
addParameter(P, 'conditionIndex', [], @isnumeric)
% figure setup
addParameter(P, 'fig', 'figure', @(x) ischar(x) && ismember(x,{'figure','clf'}))
addParameter(P, 'orient', 'square', @(x) ischar(x) && ismember(x,{'square','horizontal','vertical'}))
addParameter(P, 'axis', [], @(x) ischar(x) && ismember(x,{'square','equal'}))
addParameter(P, 'xScale', 1, @(x) isscalar(x))
addParameter(P, 'yScale', 8, @(x) isscalar(x))
addParameter(P, 'yInt', 4, @isscalar)
addParameter(P, 'pad', 0.05, @(x) isscalar(x))
addParameter(P, 'style', 'plain', @(x) ischar(x) && ismember(x,{'plain','fancy'}))
parse(P, forces, conditions, experiment, sessionIndex, varargin{:})

% get selection keys
if isempty(P.Results.conditionIndex)
    key = struct('experiment', experiment, ...
        'session_index',sessionIndex, ...
        'condition_block',P.Results.conditionBlock);
    condSel = selectrows(conditions, key);
    keys = getkeys(conditions, table2struct(condSel(:,1:4)));
else
    keys = getkeys(conditions, arrayfun(@(cIdx) struct(...
        'experiment', experiment, ...
        'session_index',sessionIndex, ...
        'condition_block', P.Results.conditionBlock,...
        'condition_index', cIdx), P.Results.conditionIndex));
end

% sort keys by condition index
keys = table2struct(sortrows(struct2table(keys), 'condition_index'));

% number of plot conditions
nCond = length(keys);

% subplot settings
eval(P.Results.fig)
switch P.Results.orient
    case 'square'
        nRow = ceil(sqrt(nCond));
        nCol = ceil(nCond/nRow);
        
    case 'horizontal'
        nRow = 1;
        nCol = nCond;
        
    case 'vertical'
        nRow = nCond;
        nCol = nCond;
end

% force range
frcMax = -Inf;
for ii = 1:length(keys)
    forceSel = selectrows(forces, keys(ii));
    frcMax = max(frcMax, max(forceSel.force{1}));
end

forceRange = [-P.Results.yInt, P.Results.yInt+(frcMax-rem(frcMax,P.Results.yInt))];

for iCo = 1:nCond
    
    subplot(nRow,nCol,iCo)
    hold on
    
    condSel = selectrows(conditions, keys(iCo));
    
    % get time vector and target force profile
    t = condSel.t_behavioral{1};
    targFrc = condSel.target_force{1};
    
    % get single-trial forces
    forceSel = selectrows(forces, keys(iCo));
    X = cell2mat(forceSel.force);
    
    if ~isempty(X)
        if ismember('trials',P.Results.layers)
            plot(t,X','color',0.8*ones(1,3));
        end
        if ismember('sem',P.Results.layers)
            mu = mean(X,1);
            ste = std(X,0,1)/sqrt(size(X,1));
            patch([t,fliplr(t)],[mu-ste,fliplr(mu+ste)],'k','facealpha',0.125,'edgealpha',0)
        end
        if ismember('mean',P.Results.layers)
            plot(t,mean(X,1),'k','LineWidth',2)
        end
    end
    if ismember('target',P.Results.layers)
        plot(t,targFrc,'c--','linewidth',2)
    end
    if ismember('count',P.Results.layers) && ~strcmp(P.Results.style,'fancy')
        text(t(1)+0.025*range(t), mean([forceRange(1),0]), sprintf('n = %i',size(X,1)))
    end
    plot(t([1 end]),[0 0],'color',.8*[1 1 1])
    if keys(iCo).stim_id == 0
        title(sprintf('condition %i\n(ID: %i)',keys(iCo).condition_index,keys(iCo).targ_id))
    else
        title(sprintf('condition %i\n(Elec. %i, Current %i)\n(ID: %i)',...
            keys(iCo).condition_index,keys(iCo).stim_electrode,...
            keys(iCo).stim_current,keys(iCo).targ_id))
    end
    xlim(t([1 end]))
    ylim(forceRange);
    box off
    if ~isempty(P.Results.axis)
        axis(P.Results.axis)
    end
    if strcmp(P.Results.style,'fancy')
        axis('square')
        title('')
        axis off
        % format x-axis
        if (strcmp(P.Results.orient,'square') && iCo*nRow>nCond)...
                || (strcmp(P.Results.orient,'horizontal'))...
                || (strcmp(P.Results.orient,'vertical') && iCo==nCond)
            
            if strcmp(keys(iCo).target_type,'STA')
                xScaleCoord = t(1)+[0,P.Results.xScale];
            else
                xScaleCoord = t(find(t>=0,1))+[0,P.Results.xScale];
            end
            plot(xScaleCoord,forceRange(1)*[1 1],'k','LineWidth',2)
            text(xScaleCoord(2)+P.Results.pad*range(t),forceRange(1),sprintf('%i s',P.Results.xScale))
        end
        % format y-axis
        if (strcmp(P.Results.orient,'square') && mod(iCo,nCol)==1)...
                || (strcmp(P.Results.orient,'horizontal') && iCo==1)...
                || strcmp(P.Results.orient,'vertical')
            
            plot((-P.Results.pad*range(t)+t(1))*[1 1],round([0,P.Results.yScale]),'k','LineWidth',2)
            text(-P.Results.pad*range(t)+t(1),0,sprintf('%i N',round(P.Results.yScale)),'Rotation',90,'HorizontalAlignment','Left','VerticalAlignment','Bottom')
        end
        set(gcf,'Color',[1 1 1])
        xlim([-P.Results.pad*range(t)+t(1), t(end)])
    end
    drawnow
end
subplot(nRow,nCol,1+(nRow-1)*nCol);

end