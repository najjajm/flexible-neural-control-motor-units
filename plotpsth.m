function plotpsth(psths, conditions, keys, varargin)
%PLOTMUSVSTIME Summary of this function goes here
% Detailed explanation goes here
FIGURE_COMPONENTS = {'Spike','Psth','Force'};

P = inputParser;
addRequired(P, 'psths', @istable)
addRequired(P, 'conditions', @istable)
addRequired(P, 'keys', @iscell)
addOptional(P, 'forces', table(), @istable)
addOptional(P, 'spikes', table(), @istable)
addParameter(P, 'HeightRatios', struct('Spike',0.2, 'Psth', 0.5, 'Force',0.2),...
    @(x) isstruct(x) && all(ismember(fieldnames(x),FIGURE_COMPONENTS))...
    && sum(cell2mat(struct2cell(x)))<=1)
addParameter(P, 'forceLim', [], @(x) isempty(x) || isnumeric(x))
addParameter(P, 'forceScale', 8, @isnumeric)
addParameter(P, 'colorSet', 'plot', @(x) ischar(x) && ismember(x,{'plot','index'}))
addParameter(P, 'xScale', 1, @isnumeric)
addParameter(P, 'psthScale', 20, @isnumeric) % FR
addParameter(P, 'padDuration', [], @(x) isempty(x) || isnumeric(x))
addParameter(P, 'layerPad', 0.03, @isnumeric)
addParameter(P, 'linky', true, @islogical)
addParameter(P, 'orient', 'square', @(x) ischar(x) && ismember(x,{'square','horizontal','vertical'}))
addParameter(P, 'fontsize', 8, @isnumeric)
parse(P, psths, conditions, keys, varargin{:})

%% Figure setup

forces = P.Results.forces;
spikes = P.Results.spikes;

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

% set heights to zero if missing data
HeightRatios = P.Results.HeightRatios;
if isempty(forces) && isfield(HeightRatios, 'Force')
    HeightRatios.Force = 0;
end
if isempty(spikes) && isfield(HeightRatios, 'Spike')
    HeightRatios.Spike = 0;
end
totHeight = sum(cell2mat(struct2cell(HeightRatios)));

fields = fieldnames(HeightRatios);
for ii = 1:length(fields)
    HeightRatios.(fields{ii}) = HeightRatios.(fields{ii}) / totHeight;
end

% setup figure struct
Fig = struct();
for ii = 1:length(fields)
    if HeightRatios.(fields{ii}) > 0
        Fig.(fields{ii}).height = HeightRatios.(fields{ii});
    end
end
compFields = fieldnames(Fig);

% vertical spacing between components
vSpace = (1-sum(cell2mat(struct2cell(HeightRatios))))/max(1,length(compFields)-1);

% psth limits
if isfield(Fig, 'Psth')
    
    psthLim = [zeros(nRow,1),P.Results.psthScale*ones(nRow,1)];
    
    for iRow = 1:nRow
        keyIdx = intersect(1:nPlt,(1:nCol)+nCol*(iRow-1));
        keySet = cell2mat(keys(keyIdx));
        
        for ii = 1:length(keySet)
            
            psthSel = selectrows(psths, keySet(ii));
            psth = psthSel.motor_unit_psth{1};
            psthSem = psthSel.motor_unit_psth_sem{1};
            
            psthLim(iRow,1) = min(psthLim(iRow,1), min(psth-psthSem));
            psthLim(iRow,2) = max(psthLim(iRow,2), max(psth+psthSem));
        end
    end
end

% force limits
if isfield(Fig, 'Force')
    
    frcLim = [zeros(nRow,1),P.Results.forceScale*ones(nRow,1)];
    
    for iRow = 1:nRow
        keyIdx = intersect(1:nPlt,(1:nCol)+nCol*(iRow-1));
        keySet = cell2mat(keys(keyIdx));
        
        for ii = 1:length(keySet)
            
            forceSel = selectrows(forces, keySet(ii));
            frc = forceSel.force{1};
            
            frcMean = mean(frc,1);
            frcSem = std(frcMean,0,1)/sqrt(size(frc,1));
            frcLim(iRow,1) = min(frcLim(iRow,1), min(frcMean-frcSem));
            frcLim(iRow,2) = max(frcLim(iRow,2), max(frcMean+frcSem));
        end
    end
end

%% Plot
figure
for iPl = 1:nPlt
    subplot(nRow,nCol,iPl)
    hold on
    
    % motor unit indexes
    psthSel = selectrows(psths, struct('experiment',keys{iPl}(1).experiment,'session_index',keys{iPl}(1).session_index));
    muIndexes = unique(table2array(psthSel(:,'motor_unit_index')));  
    
    nLayer = length(keys{iPl});
    
    cmap = getcolormap(nLayer, 'mapName', 'Set1');
    
    % spike block height
    if isfield(Fig, 'Spike')
        spkPad = 0.05 * Fig.Spike.height;
        hSpk = (Fig.Spike.height-(spkPad*(nLayer-1)))/nLayer;
    end
    
    % layer subplot components
    for iLa = 1:length(keys{iPl})
        
        % condition time
        condSel = selectrows(conditions, keys{iPl}(iLa));
        tBeh = condSel.t_behavioral{1};
        tNeu = condSel.t_neural{1};
        targDur = condSel.target_duration;
        
        if isempty(P.Results.padDuration)
            tContIdx = 1:length(tNeu);
            tSgIdx = 1:length(tBeh);
        else
            tContIdx = tNeu>=-P.Results.padDuration(iPl,1) & tNeu<=(targDur+P.Results.padDuration(iPl,2));
            tSgIdx = tBeh>=-P.Results.padDuration(iPl,1) & tBeh<=(targDur+P.Results.padDuration(iPl,2));
        end
        tNeu = tNeu(tContIdx);
        tBeh = tBeh(tSgIdx);
        
        tContRange = range(tNeu);
        tSgRange = range(tBeh);
        
        tNeu = tNeu/tContRange;
        tBeh = tBeh/tSgRange;
        
        axis([tBeh(1)-P.Results.layerPad,tBeh(end),-2*P.Results.layerPad,1])
        
        offset = 0;
        for iCf = 1:length(compFields)
            
            switch compFields{iCf}
                
                case 'Spike'
                    
                    spikesSel = selectrows(spikes, keys{iPl}(iLa));
                    spikeTrains = spikesSel.motor_unit_spikes;
                    
                    spikeTrains = cellfun(@(x) x(tContIdx), spikeTrains, 'uni',false);
                    
                    nTrial = length(spikeTrains);
                    baseline = (iLa-1)*(hSpk+spkPad);
                    for iTr = 1:nTrial
                        if any(spikeTrains{iTr})
                            plot(tNeu(spikeTrains{iTr}),...
                                baseline+hSpk*(iTr-1)/(nTrial-1),...
                                'o', 'color',cmap(iLa,:),...
                                'MarkerFaceColor',cmap(iLa,:),...
                                'MarkerSize',1)
                        end
                    end
                    
                case 'Psth'
                    
                    psthSel = selectrows(psths, keys{iPl}(iLa));
                    psth = psthSel.motor_unit_psth{1};
                    psthSem = psthSel.motor_unit_psth_sem{1};
                    
                    psth = psth(tSgIdx);
                    psthSem = psthSem(tSgIdx);
                    
                    rows = repmat((1:nRow)',1,nCol)';
                    scaleFactor = Fig.Psth.height/range(psthLim(rows(iPl),:));
                    
                    patch([tBeh,fliplr(tBeh)],offset+scaleFactor*[psth-psthSem,fliplr(psth+psthSem)],...
                        cmap(iLa,:),'EdgeAlpha',0,'FaceAlpha',0.125)
                    plot(tBeh,offset+scaleFactor*psth,'color',cmap(iLa,:),'LineWidth',1.5)
                    
                    % vertical scale bar
                    if iLa == 1
                        plot(tBeh(1)-P.Results.layerPad*[1,1],offset+scaleFactor*[0,P.Results.psthScale],'k','LineWidth',2)
                        text(tBeh(1)-2*P.Results.layerPad, offset, [num2str(P.Results.psthScale) ' spikes/s'], 'Rotation', 90, 'fontsize',P.Results.fontsize);
                    end
                    
                    % unit labels
                    unitIdx = find(keys{iPl}(iLa).motor_unit_index == muIndexes);
                    text(tBeh(1), offset+0.1+0.04*(iLa-1), sprintf('MU %i (ID %i)',unitIdx, psthSel.motor_unit_id), 'color', cmap(iLa,:),...
                        'HorizontalAlignment','Left', 'VerticalAlignment','Bottom', 'fontsize',P.Results.fontsize)
                    
                case 'Force'
                    
                    if iLa > 1
                        continue
                    end
                    
                    forceSel = selectrows(forces, keys{iPl}(iLa));
                    frc = forceSel.force{1};
                    
                    frc = frc(:,tSgIdx);
                    
                    frcMean = mean(frc,1);
                    frcSem = std(frc,0,1)/sqrt(size(frc,1));
                    
                    rows = repmat((1:nRow)',1,nCol)';
                    scaleFactor = Fig.Force.height/range(frcLim(rows(iPl),:));
                    
                    patch([tBeh,fliplr(tBeh)],offset+scaleFactor*[frcMean-frcSem,fliplr(frcMean+frcSem)],...
                        'k','EdgeAlpha',0,'FaceAlpha',0.125)
                    plot(tBeh,offset+scaleFactor*frcMean,'k','linewidth',1.5)
                    
                    % vertical scale bar
                    if iLa == 1
                        plot(tBeh(1)-P.Results.layerPad*[1,1],offset+scaleFactor*[0,P.Results.forceScale],'k','LineWidth',2)
                        text(tBeh(1)-2*P.Results.layerPad, offset, [num2str(P.Results.forceScale) ' N'], 'Rotation', 90, 'fontsize',P.Results.fontsize);
                    end
            end
            
            % increment offset
            offset = offset + Fig.(compFields{iCf}).height + vSpace;
        end
    end
    
    % format
    axis square
    
    % time scale bar
    plot(tBeh(1)+[0,P.Results.xScale]/tSgRange,-P.Results.layerPad*[1,1],'k','LineWidth',2)
    text(tBeh(1),-2*P.Results.layerPad,[num2str(P.Results.xScale) ' s'],'HorizontalAlignment','Left', 'VerticalAlignment','Top', 'fontsize',P.Results.fontsize)
    axis off
    set(gcf,'color',[1 1 1])
    
    % stim scale bar
    if condSel.stim_id > 0
        nPulse = condSel.stim_pulses;
        pulseFreq = condSel.stim_frequency;
        pulseTrainDur = nPulse/pulseFreq;
        plot([0,pulseTrainDur]/tSgRange, -2*P.Results.layerPad*[1,1],'c','LineWidth',2)
    end
    
    % title
    muscleHead = condSel.muscle_head{1};
    muscleName = condSel.muscle_name{1};
    posture = condSel.posture_id;
    if condSel.stim_id > 0
        stimElec = condSel.stim_electrode;
        stimCurr = condSel.stim_current;
        figTxt = sprintf('%s %s\nposture %i\ncondition %i (elec %i, %i uA)',...
            muscleHead, muscleName, posture, keys{iPl}(1).condition_index, stimElec, stimCurr);
    else
        figTxt = sprintf('%s %s\nposture %i\ncondition %i',...
            muscleHead, muscleName, posture, keys{iPl}(1).condition_index);
    end
    text(tBeh(1),1,figTxt,'HorizontalAlignment','Left','VerticalAlignment','Top','fontsize',P.Results.fontsize)
end
end