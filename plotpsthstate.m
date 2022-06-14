function plotpsthstate(psths, conditions, keys, varargin)
%PLOTMUSVSTIME Summary of this function goes here
% Detailed explanation goes here

P = inputParser;
addRequired(P, 'psths', @istable)
addRequired(P, 'conditions', @istable)
addRequired(P, 'keys', @iscell)
addParameter(P, 'linkAxes', 'all', @(x) iscell(x) || iscode(x,{'none','all'}))
addParameter(P, 'orient', 'square', @(x) iscode(x,{'square','horizontal','vertical'}))
addParameter(P, 'padDuration', [], @(x) isempty(x) || isnumeric(x))
addParameter(P, 'pad', 0.02, @isnumeric)
parse(P, psths, conditions, keys, varargin{:})

% SETUP -------------------------------------------------------

% ensure proper inputs
assert(all(cellfun(@(x) all(cellfun(@length,x)==2), keys)),...
    'All condition key sets must be of length 2')

% figure orientation
nPlt = length(keys);
switch P.Results.orient
    case 'square'
        nCol = ceil(sqrt(nPlt));
        nRow = ceil(nPlt/nCol);
        
    case 'horizontal'
        nCol = nPlt;
        nRow = 1;
        
    case 'vertical'
        nCol = 1;
        nRow = nPlt;
end

% number of conditions
nCond = cellfun(@length,keys);

% PLOT --------------------------------------------------------
figure
ax = [];
for iPl = 1:nPlt
    ax(iPl) = subplot(nRow,nCol,iPl);
    hold on
    
    % condition colormap
    cmap = getcolormap(nCond(iPl),'mapName','Set2');
    
    % figure handle and legend text
    fh = zeros(1,nCond(iPl));
    legTxt = cell(1,nCond(iPl));
    
    % layer subplot components
    rMax = 0;
    for iCo = 1:nCond(iPl)
        
        % fetch psths
        psth = cell(1,2);
        for ii = 1:2
            psthSel = selectrows(psths, keys{iPl}{iCo}(ii));
            psth(ii) = psthSel.motor_unit_psth;
        end
        
        % condition selection (assuming the same for both keys)
        condSel = selectrows(conditions, keys{iPl}{iCo}(1));
        
        % truncate
        if ~isempty(P.Results.padDuration)
            t = condSel.t_behavioral{1};
            targDur = condSel.target_duration;
            
            padDur = P.Results.padDuration(sum(nCond(1:iPl-1))+iCo, :);
            tIdx = t>=-padDur(1) & t<=(targDur+padDur(2));
            psth = cellfun(@(x) x(tIdx), psth, 'uni',false);
        end
        
        fh(iCo) = plot(psth{1},psth{2},'color',cmap(iCo,:),'linewidth',2);
        
        rMax = max(1, max(rMax, ceil(max(cellfun(@max,psth))/5)*5));
        
        switch condSel.target_type{1}
            case 'STA'
                legTxt{iCo} = sprintf('Cond %i (STA %i N), block %i',...
                    condSel.condition_index,...
                    condSel.force_max * condSel.target_offset,...
                    condSel.condition_block);
                
                if condSel.stim_id > 0
                    legTxt{iCo} = [legTxt{iCo}, sprintf(', elec %i (%i uA)',...
                        condSel.stim_electrode,...
                        condSel.stim_current)];
                end
                
            case 'RMP'
                legTxt{iCo} = sprintf('Cond %i (RMP %i N/s), block %i',...
                    condSel.condition_index,...
                    (condSel.force_max * condSel.target_amplitude) / condSel.target_duration,...
                    condSel.condition_block);
                
            case 'SIN'
                legTxt{iCo} = sprintf('Cond %i (SIN %.2f Hz), block %i',...
                    condSel.condition_index,...
                    condSel.target_frequency1,...
                    condSel.condition_block);
                
            case 'CHP'
                legTxt{iCo} = sprintf('Cond %i (CHP %.2f - %.2f Hz), block %i',...
                    condSel.condition_index,...
                    condSel.target_frequency1,...
                    condSel.target_frequency2,...
                    condSel.condition_block);
                
        end
    end
    
    % firing rate scale bar
    rScale = ceil(rMax*0.25/10)*10;
    plot([0,rScale], -P.Results.pad*rMax*[1,1], 'k', 'LineWidth',2)
    plot(-P.Results.pad*rMax*[1,1], [0,rScale], 'k', 'LineWidth',2)
    text(0, -P.Results.pad*rMax, sprintf('%i spikes/s',rScale), 'HorizontalAlignment','Left', 'VerticalAlignment','Top')
    text(-P.Results.pad*rMax, 0, sprintf('%i spikes/s',rScale), 'HorizontalAlignment','Left', 'VerticalAlignment','Bottom', 'Rotation',90)
    
    % condition legend
    legend(fh,legTxt{:})
    
    % title
    figTxt = sprintf('%s %s\nMU %i vs %i',...
        condSel.muscle_head{1},...
        condSel.muscle_name{1},...
        keys{iPl}{1}(2).motor_unit_id,...
        keys{iPl}{1}(1).motor_unit_id);
    text(0,1,figTxt,'Units','normalized','HorizontalAlignment','left','VerticalAlignment','top')
    
    % format
    axis([-P.Results.pad*rMax,rMax,-P.Results.pad*rMax,rMax])
    axis square
    axis off
end

if ischar(P.Results.linkAxes)
    if strcmp(P.Results.linkAxes,'all')
        linkaxes(ax,'xy')
    end
else
    for ii = 1:length(P.Results.linkAxes)
        linkaxes(ax(P.Results.linkAxes{ii}),'xy')
    end
end

set(gcf,'Color',[1,1,1])

end

