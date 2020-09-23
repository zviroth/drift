
close all
clear all
illusion = 0;

expNum=2;
crossValidation = 'trial';%'trial','run'
locThresh = 0.2; %only include voxels with coherence>locThresh, for each of the two localizers.

saveFolder = '~/noah/';
switch expNum
    case 1
        expName = '3conds';
    case 2
        expName = 'attn';
    case 3
        expName = '4conds';
    case 4
        expName = '3tdata';
end
thrsh = num2str(locThresh,'%.2f');
locStr = thrsh([1 3:4]);
directionOrIllusionStr = 'direction';
if illusion
    directionOrIllusionStr = 'illusion';
end
load([saveFolder 'decodeDrift_' directionOrIllusionStr '_' expName '_' crossValidation '_thresh' locStr '.mat'],'crossValidation',...
    'r2thresh','locThresh','ROIs','subNames','dataFolder',...
    'subRoiNans','subNumTrials','subNumScans',...
    'subAcc','subRandAcc','subRand95Acc','nperms','locSharedVox','totalVox');
if illusion
    subAcc = squeeze(mean(subAcc,4));
    subRandAcc = squeeze(mean(subRandAcc,4));
    
    
end
numRois = length(ROIs);
plotColors = {[0.6 0.6 1], [1 0.6 0.6]};
nanColor = [0 0 0];
facealpha = 0.3;
ilocTitle = {'stim','eye'};
randDist = squeeze(mean(subRandAcc,1));
meanAcc = squeeze(mean(subAcc,1));
for iroi=1:numRois
    for iloc=1:2
%         decodingAcc(isub,iroi,iloc) = res{iloc}.accuracy;
%         randAccDist(isub,iroi,iloc,:) = res{iloc}.randAccuracy;
        
        temp = subRandAcc(:,iroi,iloc,:);
%         randDist(iroi,iloc,:) = temp(:);
        
%         meanAcc(iroi,iloc) = mean(decodingAcc(:,iroi,iloc));
        rand95(iroi,iloc) = prctile(randDist(iroi,iloc,:),95);
    end
end

maxRandAcc = max(randDist,[],3);
minRandAcc = min(randDist,[],3);
roiInd = 1:numRois;

figure(1)
rows=2;
cols=1;
badRoi= zeros(2,numRois);

for iloc=1:2
    for iroi=1:numRois
       if sum(isnan( subRandAcc(:,iroi,iloc)))>0
          %mark this ROI
          badRoi(iloc,iroi) = 1;
          randDist(iroi,iloc,:) = randDist(1,iloc,:);
          facecolorVec(iloc,iroi,:) = nanColor;
       else
           facecolorVec(iloc,iroi,:) = plotColors{iloc};
       end
    end
    subplot(rows,cols,iloc)
%     violin(squeeze(randDist(:,iloc,:))','xlabel',ROIs,'mc',[]);
%     plot(1:numRois, rand95(:,iloc),'.','markersize',40,'color',[0 0 0]);

    

%     randDist(badRoi(iloc,:)==1,iloc,:) = ones(sum(badRoi(iloc,:)) randDist(1,iloc,:);
%     randDist(badRoi(iloc,:)==1,iloc,1) = 0.05;
%     bw = 
    myviolin(squeeze(randDist(:,iloc,:))','xlabel',ROIs,'facecolor',squeeze(facecolorVec(iloc,:,:)),'facealpha',facealpha);
%     myviolin(squeeze(randDist(badRoi(iloc,:)==0,iloc,:))','xlabel',ROIs,'facecolor',plotColors{iloc},'facealpha',0.1);

%     meanAcc(badRoi(iloc,:)==1,iloc) = 0.5;
    plot(1:numRois, meanAcc(:,iloc),'.','markersize',40,'color',plotColors{iloc});
%     plot(1:numRois, meanAcc(badRoi(iloc,:)==0,iloc),'.','markersize',40,'color',plotColors{iloc});

    ylim([0.3 0.7]);
    xticks(1:numRois);
    xticklabels(ROIs);
    title(ilocTitle{iloc});
    ylabel('decoding accuracy');
end
set(gcf,'position',[150 150 1100 600]);


figure(2)

for iloc=1:2
    shift = (iloc-1.5)*0.2;
   temp = plot([ roiInd+shift; roiInd+shift], [minRandAcc(:,iloc)'; maxRandAcc(:,iloc)'], 'linewidth',5,'color',plotColors{iloc});
   p(iloc) = temp(1);%for legend, keep plot properties of first ROI
   hold on
   plot([ roiInd+shift-0.3; roiInd+shift+0.3], [rand95(:,iloc)'; rand95(:,iloc)'], 'linewidth',3,'color',plotColors{iloc});
   p(3) = plot(roiInd+shift, meanAcc(:,iloc), '.','markersize',15,'color',[0 0 0]);
end
xlim([0 numRois+1]);
ylim([0.3 0.7]);
xlabel('ROI');
xticks([1:numRois])
xticklabels(ROIs);
ylabel('decoding accuracy');
L=legend([p(1) p(2) p(3)],'stim','eye','data');
title([expName ', ' crossValidation ' cv, ' num2str(size(subRandAcc,1)) ' subjects']);
set(gcf,'position',[100 1500 1000 400]);