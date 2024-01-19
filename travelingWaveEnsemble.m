function [lateSpikeEnsemble, nolateSpikeEnsemble] = travelingWaveEnsemble(spikeTrials,lateSpikeTrials,nolateSpikeTrials,ROIcentroid,AverageImage)
% Parse spike trials into late/no late spike and then calculate the
% ensemble space that it occupies. Using this we can then reproject into
% the broad ensemble space and (hopefully) cluster the data

%% FOV 
count = 1;
for i = lateSpikeTrials
    lateSpikes{count} = spikeTrials{i};
    count = count+1;
end

count = 1;
for i = nolateSpikeTrials
    nolateSpikes{count} = spikeTrials{i};
    count = count+1;
end

if iscell(lateSpikes)
    lateSpikes = horzcat(lateSpikes{:});
end

if iscell(nolateSpikes)
    nolateSpikes = horzcat(nolateSpikes{:});
end

% Now extract representative ensembles
% Late Spike
FactorCorrection = 5*floor(size(lateSpikes,2)/5); % Correct for frame size aquisition
lateSpikeEnsemble = ensembleAnalysis(lateSpikes(:,1:FactorCorrection),ROIcentroid);
[~,I] = sort(cellfun(@length,lateSpikeEnsemble.NodeList),'descend'); %sort max node size
rankEdges = lateSpikeEnsemble.NumEdges(:,I);
rankEnsembles = lateSpikeEnsemble.NodeList(:,I); 
[grad,~]=colorGradient([1 0 0],[0 0 0],6);
lateSpikeEnsemble.rankEnsembles = rankEnsembles;
lateSpikeEnsemble.rankEdges = rankEdges;
lateSpikeEnsemble.Spikes = lateSpikes;
try
    figure,
    for i = 2:5
        axis off
        color = jet(3);
        EnsembleMap(AverageImage,ROIcentroid,lateSpikeEnsemble.rankEnsembles{i},5,grad(i,:))
        set(gcf,'Position',[100 100 500 500])
        drawnow
        hold on
    end
catch
end
% No Late Spike
FactorCorrection = 5*floor(size(nolateSpikes,2)/5); % Correct for frame size aquisition
nolateSpikeEnsemble = ensembleAnalysis(nolateSpikes(:,1:FactorCorrection),ROIcentroid);
[~,I] = sort(cellfun(@length,nolateSpikeEnsemble.NodeList),'descend'); %sort max node size
rankEdges = nolateSpikeEnsemble.NumEdges(:,I);
rankEnsembles = nolateSpikeEnsemble.NodeList(:,I); 
[grad,~]=colorGradient([0 0 1],[0 0 0],6);
nolateSpikeEnsemble.rankEnsembles = rankEnsembles;
nolateSpikeEnsemble.rankEdges = rankEdges;
nolateSpikeEnsemble.Spikes = nolateSpikes;
try
    figure,
    for i = 1:5
        axis off
        color = jet(3);
        EnsembleMap(AverageImage,ROIcentroid,nolateSpikeEnsemble.rankEnsembles{i},4,grad(i,:))
        set(gcf,'Position',[100 100 500 500])
        drawnow
        hold on
    end
catch
end

% Combine Maps
figure,imagesc(interp2(nolateSpikeEnsemble.sim_index,2)),colormap(jet),caxis([0.13 .8])
K = (1/25)*ones(5);
figure,imagesc(interp2(conv2(nolateSpikeEnsemble.sim_index,K,'same'),2)),colormap(jet),caxis([.08 .3])
