%% Remove ROIs
if exist('badComponents','var') && ~exist('badComFlag','var')
[DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A] = ...
    removeROI(DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A,unique(badComponents));
badComFlag = 1;
end
%% Fix centroids
ROIcentroid = [];
for i = 1:length(ROI)
    blah = vertcat(ROI{i}{:});
    ROIcentroid(i,:) = floor(mean(blah,1));
end
%% Analysis
addpath(genpath('main'));
std_threshold = 6;
static_threshold = .01;
Spikes = Spike_Detector_Single(dDeltaFoverF,std_threshold,static_threshold);
%Excude inactive cells
% numSpikes = sum(Spikes,2);
% keepSpikes = find(numSpikes>(.01*mean(numSpikes)));
% Spikes = Spikes(keepSpikes,:);
[coactive_cells,detected_spikes] = coactive_index(Spikes,length(Spikes));
cell_count = length(ROI);
time = time_adjust(size(DeltaFoverF,2),30);
for i = 1:size(DeltaFoverF,1)
    calcium_avg{i} = STA(DeltaFoverF(i,:),2,120);
end

% Perform shuffling and pairwise if data is small enough
if size(DeltaFoverF,2)<2000    
    factorCorrection = 10*floor(size(Spikes,2)/10); % Correct for frame size aquisition
    [vectorized,sim_index] = cosine_similarity(Spikes(:,1:factorCorrection),10);
    corr = correlation_dice(Spikes);
    Connected_ROI = Connectivity_dice(corr, ROI,0.3);
    [NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
        ActivityCentroid,ActivityCentroidVariance]...
        = Network_Analysis(ROIcentroid,Connected_ROI);
end
% Pairwise Velocity Analysis
% velocityPairwise(VR_data,Spikes)
% Ensemble Analysis
% figure,[Coor,json_file] = plot_contours(A,C,ops,0); % contour plot of spatial footprints
factorCorrection = 5*floor(size(Spikes,2)/5); % Correct for frame size aquisition
Ensemble = ensembleAnalysis(Spikes(:,1:factorCorrection),ROI,ROIcentroid);
% Ensemble = ensembleNetworks(Ensemble);
% Plot Ensemble
% ensembleVid(Ensemble,AverageImage,ROIcentroid,files);
% Displays ensemble overlay
[~,I] = sort(cellfun(@length,Ensemble.NodeList),'descend'); %sort max node size
rankEdges = Ensemble.NumEdges(:,I);
rankEnsembles = Ensemble.NodeList(:,I); 
[grad,~]=colorGradient([1 0 0],[0 0 0],3);
Ensemble.rankEnsembles = rankEnsembles;
figure,
for i = 1
    axis off
    color = jet(3);
    EnsembleMap(AverageImage,ROIcentroid,Ensemble.rankEnsembles{i},4,grad(i,:))
    set(gcf,'Position',[100 100 500 500])
    drawnow
    hold on
end
% Combine Maps
figure,imagesc(interp2(Ensemble.sim_index,2)),colormap(jet),caxis([0.13 .4])
K = (1/25)*ones(5);
figure,imagesc(interp2(conv2(Ensemble.sim_index,K,'same'),2)),colormap(jet),caxis([.08 .3])

%
sizeE = cellfun(@size,Ensemble.rankEnsembles,'UniformOutput',false);
sizeE = cell2mat(sizeE);
sizeE(sizeE==1) = [];
sizeE = sizeE/max(sizeE);
figure,plot(sizeE),title('Ranked Ensembles')

sizeEdge = cell2mat(rankEdges);
sizeEdge = sizeEdge/max(sizeEdge);
figure,plot(sizeEdge),title('Ranked Connections')
% Ensemble Stats
% EnsembleStats
% Information Entropy
informationEntropy = shannonEntropy(Ensemble.rankEnsembles);
% Plot Centroid Boundary
figure,hold on
%Rank by number of Activity points
[~,I] = sort(cellfun(@length,Ensemble.ActivityCoords),'descend'); %sort max node size
rankedActivityCoords = Ensemble.ActivityCoords(:,I);
Ensemble.rankedActivityCoords = rankedActivityCoords;
% Now look at trial specific changes
% Standard trial length is 360 frames ** old data set had 240 frames
%%
spikeTrials = [];
trialLength = 360;
for i = 1:size(Spikes,2)/trialLength
    spikeTrials{i} = Spikes(:,((i-1)*trialLength+1):i*trialLength);
    DeltaTrials(:,:,i) = DeltaFoverF(:,((i-1)*trialLength+1):i*trialLength);
end

% check for edge case when calcium and trial numbers don't match
if length(spikeTrials)<trialData.responsiveTrials.trialNum(end)
    trialData.responsiveTrials.lateSpikeTrialsOld = trialData.responsiveTrials.lateSpikeTrials;
    trialData.responsiveTrials.lateSpikeTrials = trialData.responsiveTrials.lateSpikeTrials...
        (trialData.responsiveTrials.lateSpikeTrials<=length(spikeTrials));
    
    trialData.responsiveTrials.noLateSpikeTrialsOld = trialData.responsiveTrials.noLateSpikeTrials;
    trialData.responsiveTrials.noLateSpikeTrials = trialData.responsiveTrials.noLateSpikeTrials...
        (trialData.responsiveTrials.noLateSpikeTrials<=length(spikeTrials));
end

figure,
 for i = 1:length(spikeTrials)
   subplot(6,6,i),Show_Spikes(spikeTrials{i})
 end
% Late vs no Late spike ensembles
[lateSpikeEnsemble, nolateSpikeEnsemble] =...
    travelingWaveEnsemble(spikeTrials,trialData.responsiveTrials.lateSpikeTrials,trialData.responsiveTrials.noLateSpikeTrials,ROIcentroid,AverageImage);

% manifold analysis and entropy
lateSpikeEnsemble = ensembleMetric(lateSpikeEnsemble,AverageImage,ROIcentroid);
nolateSpikeEnsemble = ensembleMetric(nolateSpikeEnsemble,AverageImage,ROIcentroid);

% some statistics about these ensembles
lateSpikeEnsemble = ensembleStat(lateSpikeEnsemble);
nolateSpikeEnsemble = ensembleStat(nolateSpikeEnsemble);
% Quantify Node Reactivation
lateSpikeEnsemble = findCritNode(trialData.responsiveTrials.lateSpikeTrials,ROI,spikeTrials,lateSpikeEnsemble); % trials, ROI, spikeTrials, Spikes, Ensemble
nolateSpikeEnsemble = findCritNode(trialData.responsiveTrials.noLateSpikeTrials,ROI,spikeTrials,nolateSpikeEnsemble);

%% Trial specific connectivity
simM = [];
Connected_ROI = [];
Ls = [];
% figure,
count = 1;
for i = trialData.responsiveTrials.lateSpikeTrials  %t1
    [~,sim_index] = cosine_similarity(spikeTrials{i},10);
%     subplot(5,4,count),imagesc(sim_index),colormap(jet)
    Ls(count) = mean(sim_index(sim_index),'all'); %mean center
    count = count+1;
end
% t1 = vertcat(simValue{:});
%  figure,boxplot(Ls);ylim([0.0 0.5])
 lateSpikeEnsemble.LsSim = Ls;

Connected_ROI = {};
simM = [];
simValue = [];
nLs = [];
count = 1;
% % figure,
for i = trialData.responsiveTrials.noLateSpikeTrials%t2
    [~,sim_index] = cosine_similarity(spikeTrials{i},10);
%     subplot(5,4,count),imagesc(sim_index),colormap(jet)
    nLs(count) = mean(sim_index(sim_index>0),'all');
    count = count+1;
end
% % t2 = vertcat(simValue{:});
% figure,boxplot(nLs);ylim([0.0 0.5])
nolateSpikeEnsemble.nLsSim = nLs;

close all
%% Overlapping ensemble metric
% Check number of shared nodes by adding each interative ensemble rank
checkIteration = length(lateSpikeEnsemble.rankEnsembles)- length(nolateSpikeEnsemble.rankEnsembles);
if checkIteration>0 %ie more late spike ensembles
    ensembleIteration = length(nolateSpikeEnsemble.rankEnsembles); %only index across smaller one (and we'll negate the remaining number
else %ie more no late spike ensembles
    ensembleIteration = length(lateSpikeEnsemble.rankEnsembles);
end

sharedEnsemble = [];
count = 1;
for i = 1:40
    lsEnsemble = lateSpikeEnsemble.rankEnsembles{i};
    nlsEnsemble = nolateSpikeEnsemble.rankEnsembles{i};
    checkIdx = length(lsEnsemble)-length(nlsEnsemble); %find which one is larger and check existing indices
    if ~isempty(checkIdx)
        if checkIdx>0
            idx = ismember(lsEnsemble,nlsEnsemble);
        else
            idx = ismember(lsEnsemble,nlsEnsemble);
        end
        if count>1
            sharedEnsemble(count) = sum(idx)+sharedEnsemble(count-1);
        else
            sharedEnsemble(count) = sum(idx);
        end
        count = count+1;
    end
end

%%
lateSpikeConnectedROI = vertcat(lateSpikeEnsemble.Connected_ROI{:});
[r,c] = find(lateSpikeConnectedROI==0);
lateSpikeConnectedROI(r,:) = [];
nolateSpikeConnectedROI = vertcat(nolateSpikeEnsemble.Connected_ROI{:});
[r,c] = find(nolateSpikeConnectedROI==0);
nolateSpikeConnectedROI(r,:) = [];
figure('Name','Late Spike Network Map'); NodeSize = 3;EdgeSize = 1;Cell_Map_Dice(AverageImage,lateSpikeConnectedROI,ROIcentroid,NodeSize,EdgeSize)
figure('Name','No Late Spike Network Map'); NodeSize = 2;EdgeSize = 2;Cell_Map_Dice(AverageImage,nolateSpikeConnectedROI,ROIcentroid,NodeSize,EdgeSize)
%%
W12_10Entropy.informationEntropy = informationEntropy;
W12_10Entropy.rankedEnsembles = rankEnsembles;
W12_10Entropy.rankedEdges = rankEdges;
W12_10Entropy.Ensemble = Ensemble;

%% SVD/PCA of Ensembles
[vectorizedL,sim_indexL] = cosine_similarity(lateSpikeEnsemble.Spikes,10);
[vectorizedNL,sim_indexNL] = cosine_similarity(nolateSpikeEnsemble.Spikes,40);

comVect = [vectorizedL vectorizedNL];
[tProjq1, tProjq2, uProjq1, uProjq2] = featureProject(comVect,length(vectorizedL),0);
legend('Late Spike','No Late Spike')

comSim = [sim_indexL sim_indexL];
svd_analysis(comSim)
%% Trial by Trial analysis ##Only use with batch processed files##
addpath(genpath('Figures'));
[batchSpikes,batch_corr] = TrialByTrial(batchData([1,2,4])); % Function call
bin = 20;
[vectorized,sim_index] = cosine_similarity(batchSpikes,bin);
[z,mu,sigma] = zscore(sim_index);
figure('Name','Cosine-Similarity Index'); h = htmp(sim_index,100);
caxis([0 0.7]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
 figure('Name','Dice Correlation')
 for i = 1:size(batch_corr,3)
     subplot(2,3,i),h = htmp(batch_corr(:,:,i),20);caxis([0 0.4]);
 end
%% Plot all the Figures
figure('Name','Network Map'); NodeSize = 2;EdgeSize = 2;Cell_Map_Dice(AverageImage,Connected_ROI,ROIcentroid,NodeSize,EdgeSize)

%% Rotary Encoder
figure('Name','Pulse Data');plot(encoder_data.rotate_pulse);
figure('Name','Angular Distance');bar(encoder_data.ang_distance);
figure('Name','Angular Velocity');bar(encoder_data.ang_velocity,'FaceColor',[.16 .835 .384],'EdgeColor','none');
figure('Name','Avg. Angular Velocity');avgV = movmean(encoder_data.ang_velocity,2);bar(avgV,'FaceColor',[.16 .835 .384],'EdgeColor','none');

%%
A    = M2;
imwrite(A(:, :, 1), 'test.tiff');
for k = 2:size(A, 3)
  imwrite(A(:, :, k), 'test.tiff', 'WriteMode', 'append');
end
image_movie = mat2gray(M2);
implay(image_movie);

%%
figure,plot(strongConnections(:,8))
hold on,plot(weakConnections(:,8))
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
%%
plot_raster(1:120,Spikes(5,1:120))
% Have you tried using Multidimensional Scaling (MDS) to emebed the
% centroids in a 2 dimensional space for visualization?


