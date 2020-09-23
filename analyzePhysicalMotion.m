close all
mrQuit
clear all

crossValidation = 'trial';%'trial','run'
saveFolder = '~/noah/';
load([saveFolder 'decodePhysicalDirection_' crossValidation '.mat'],'crossValidation',...
    'r2thresh','ROIs','subNames','dataFolder',...
    'subRoiNans','subNumTrials','subNumScans',...
    'subAcc','subRandAcc','subRand95Acc','nperms','totalVox');
numSubs = length(subNames);
numRois = length(ROIs);
facealpha = 0.3;
randDist = squeeze(mean(subRandAcc,1));
meanAcc = squeeze(mean(subAcc,1));
maxRandAcc = max(randDist,[],2);
minRandAcc = min(randDist,[],2);
roiInd = 1:numRois;
for iroi=1:numRois
        rand95(iroi) = prctile(randDist(iroi,:),95);
        rand50(iroi) = prctile(randDist(iroi,:),50);
        rand05(iroi) = prctile(randDist(iroi,:),5);
end

rand50color = [0.2 0.2 0.2];
rand95color = [0.4 0 0];
meanColor = [0 0 0];
subColor = [0.6 0.6 0.6];

plotColors = {[0.2 0.2 1], [1 0.2 0.2], [0.2 1 0.2], [0.2 0.6 0.6]};
nanColor = [0 0 0];
figure(2)
roiInd = 1:numRois/2;
hemiROIs = {'lV1','lV2','lV3','lhV4','lLO1','lLO2','lTO1','lTO2','lV3a','lV3b','lV01','lV02'};

temp = plot([ roiInd; roiInd], [minRandAcc(roiInd)'; maxRandAcc(roiInd)'], 'linewidth',3,'color',subColor);
hold on
plot([ roiInd-0.3; roiInd+0.3], [rand95(roiInd); rand95(roiInd)], ':','linewidth',2,'color',subColor);
plot([ roiInd-0.3; roiInd+0.3], [rand50(roiInd); rand50(roiInd)], '-','linewidth',3,'color',subColor);
roiX = repmat(roiInd,numSubs,1);
plot(roiX+(rand(size(roiX))-0.5)/6,subAcc(:,roiInd),'o','markersize',5,'color',plotColors{1});
plot(roiInd, meanAcc(roiInd), '.','markersize',40,'color',plotColors{1});

xlim([0 length(roiInd)+1]);
% ylim([0.3 0.7]);
xlabel('ROI');
xticks([roiInd])
xticklabels(hemiROIs);
ylabel('decoding accuracy');
% L=legend([p(1) p(2) p(3)],'stim','eye','data');
% title([expName, ' ', directionOrIllusionStr, ' ', crossValidation ' cv, ' num2str(size(subRandAcc,1)) ' subjects, co>' num2str(locThresh)]);
set(gcf,'position',[100 1500 1000 400]);




figure(3)
roiInd = 1:numRois;

   temp = plot([ roiInd; roiInd], [minRandAcc(:)'; maxRandAcc(:)'], 'linewidth',3,'color',plotColors{1}.^0.5);
   hold on
   plot([ roiInd-0.3; roiInd+0.3], [rand95(:)'; rand95(:)'], 'linewidth',3,'color',plotColors{1}.^0.5);
   plot([ roiInd-0.3; roiInd+0.3], [rand05(:)'; rand05(:)'], 'linewidth',3,'color',plotColors{1}.^0.5);
   plot(roiInd, meanAcc(:), '.','markersize',15,'color',[0 0 0]);

xlim([0 numRois+1]);
% ylim([0.3 0.7]);
xlabel('ROI');
xticks([1:numRois])
xticklabels(ROIs);
ylabel('decoding accuracy');
% if expNum==1 | expNum==2
%     L=legend([p(1) p(2) p(3)],'stim','eye','data');
% elseif expNum==3
%     L=legend([p(1) p(2) p(3)],'local+global','local','data');
% end
% title([expName, ' ', directionOrIllusionStr, ' ', crossValidation ' cv, ' num2str(size(subRandAcc,1)) ' subjects, co>' num2str(locThresh)]);
set(gcf,'position',[100 1500 1000 400]);
