close all
clear all

expNum=1;


r2thresh = 50;%r2 percentile
locThresh = 0.2;


saveFolder = '~/noah/';
switch expNum
    case 1
        expName = '3conds';
    case 2
        expName = 'attn';
end

load([saveFolder 'driftTS_' expName '.mat'],'ROIs','subNames','dataFolder',...
    'stimLoc','eyeLoc','driftBetas','subRoiNans','subNumTrials','subNumScans');
numSubs = length(subNames);
numRois = length(ROIs);

temp = trialReg.r2>prctile(trialReg.r2,50) & stimlocCorrData{1}.co(~roiNans)>locThresh & eyelocCorrData{1}.co(~roiNans)>locThresh;
locSharedVox(isub,iroi) = sum(temp);
totalVox(isub,iroi) = length(temp);

for isub=1:numSubs
    ntrials = subNumTrials(isub);
    for iroi=1:numRois
        
    end
end




% res{1}.totalVox = length(temp);
roiNans = subRoiNans{isub,iroi};
trialReg = driftBetas{isub,iroi};
for iloc=1:2
    switch iloc
        case 1
            locCorrData = stimLoc{isub,iroi};
        case 2
            locCorrData = eyeLoc{isub,iroi};
    end
    %     keyboard
    locCorrData.co = locCorrData.co(~roiNans);%data from concatenation excluded these voxels
    goodVox = trialReg.r2>prctile(trialReg.r2,50) & locCorrData.co>locThresh;
    
    locNumVox(isub,iroi,iloc) = sum(goodVox);
%     res{iloc}.numVox = sum(goodVox);
    [roiname ' loc ' num2str(iloc) ': using ' num2str(locNumVox(isub,iroi,iloc)) ' out of ' num2str(totalVox(isub,iroi)) ' voxels']
    if sum(goodVox)>0
        
%                 betas = trialReg.ehdr(goodVox,[1:nTrials (nTrials*2)+1:nTrials*3])'; %this is for 1 vs 3
        for iLR = 1:2
            switch iLR
                case 1
                    betas = trialReg.ehdr(:,[(nTrials)+1:nTrials*2 (nTrials*2)+1:nTrials*3])'; %this is for 2 vs 3
                case 2
                    betas = trialReg.ehdr(:,[1:nTrials (nTrials)+1:nTrials*2])'; %this is for 1 vs 2
            end
            %betas = trialReg.ehdr(:,[(nTrials)+1:nTrials*2 (nTrials*2)+1:nTrials*3])'; %this is for 2 vs 3
            %betas = trialReg.ehdr(:,[1:nTrials (nTrials)+1:nTrials*2])'; %this is for 1 vs 2
            betas = zscore(betas);
            
            % create a grouping variable
            group = cat(1, repmat('l', nTrials, 1), repmat('r', nTrials, 1));
            
            %     disp(sprintf('(plotROI): Computing accuracy for L vs R'));
            % create partitioning scheme
            CVO = cvpartition(group, 'KFold', nTrials*2);
            
            if strcmp(crossValidation, 'trial')
                err = zeros(CVO.NumTestSets,1);
                for i = 1:CVO.NumTestSets
                    trIdx = CVO.training(i);
                    teIdx = CVO.test(i);
                    ytest = classify(betas(teIdx,:),betas(trIdx,:),group(trIdx), 'diagLinear');
                    err(i) = sum(~strcmp(ytest,group(teIdx)));
                end
                accuracy = 1 - sum(err)/sum(CVO.TestSize);
                disp(sprintf('(plotROI): Classification accuracy: %f', accuracy));
                
            elseif strcmp(crossValidation, 'run')
                
                clear accuracy; clear ytest;
                nRuns = length(group)/2/4;
                runVec = repmat(1:nRuns, 4, 1);
                runVec = runVec(:);
                runVec = cat(1, runVec, runVec);
                for iRun = 1:nRuns
                    trIdx = runVec ~= iRun;
                    teIdx = runVec == iRun;
                    ytest(teIdx) = classify(betas(teIdx,:),betas(trIdx,:),group(trIdx), 'diagLinear');
                end
                accuracy = sum(ytest==group') / length(group);
                disp(sprintf('(plotROI): Classification accuracy: %f', accuracy));
            else
                keyboard
            end
            
            
            
            % New way of randomizing labels - only randomize ONCE per run, and have
            % equal amounts of L and R WITHIN run
            
            nScans = subNumScans(isub);
            
            % trialsPerCond will be the trials of BOTH CONDITIONS
            trialsPerCond = length(group) / nScans;
            groupByScan = cat(1, repmat('l', trialsPerCond/2, 1), repmat('r', trialsPerCond/2, 1));
            
            % We have to create 'length(nScans)' amount of randomized labels
            clear randLabels
            for iBoot = 1:nperms
                for iScan = 1:nScans
                    tempList = groupByScan(randperm(size(groupByScan,1)));
                    idx1 = (iScan*4)-3;
                    idx2 = iScan*4;
                    randLabels(idx1:idx2,iBoot) = tempList(1:4);
                    randLabels(idx1+length(group)/2 : idx2+length(group)/2,iBoot) = tempList(5:8);
                end
            end
            
            %         disp(sprintf('(plotROI): Computing bootstrap for L vs R'));
            %         disp(sprintf('(plotROI): (the new, Martin-version of bootstrapping!)'))
            disp(sprintf('(plotROI): Classification method = %swise', crossValidation));
            nBoot = 1000;
            bootAccuracy = zeros(nBoot, 1);
            for iBoot=1:nBoot;
                gNull = randLabels(:,iBoot);
                if strcmp(crossValidation, 'trial')
                    for i = 1:CVO.NumTestSets
                        trIdx = CVO.training(i);
                        teIdx = CVO.test(i);
                        ytest = classify(betas(teIdx,:),betas(trIdx,:),gNull(trIdx), 'diagLinear');
                        err(i) = sum(~strcmp(ytest,gNull(teIdx)));
                    end
                    bootAccuracy(iBoot) = 1 - sum(err)/sum(CVO.TestSize);
                    
                elseif strcmp(crossValidation, 'run')
                    clear ytest;
                    for iRun = 1:nRuns
                        trIdx = runVec ~= iRun;
                        teIdx = runVec == iRun;
                        ytest(teIdx) = classify(betas(teIdx,:),betas(trIdx,:),gNull(trIdx), 'diagLinear');
                    end
                    bootAccuracy(iBoot) = sum(ytest==group') / length(group);
                else
                    keyboard
                end
                
            end
            disp(sprintf('(plotROI): Bootstrapped threshold is: %f', prctile(bootAccuracy, 95)));
            
            % close all
            if toPlot
                h = figure(2);
                clf
                hist(bootAccuracy);
                vline(accuracy, 'k');
                vline(prctile(bootAccuracy, 95), 'r');
                title('black = accuracy, red = 95% of null distribution')
                %         filename = sprintf('Figures/nullDistribution_%s_%swise.pdf', roiname, crossValidation);
                %         savepdf(h, filename)
            end
            res{iloc}.accuracy(iLR) = accuracy;
            res{iloc}.randAccuracy(:,iLR) = bootAccuracy;
            res{iloc}.rand95accuracy(iLR) = prctile(bootAccuracy, 95);
        end
    else
        'ZERO GOOD VOXELS!!'
        res{iloc}.accuracy(1:2) = NaN;
        res{iloc}.randAccuracy(:,2) = NaN(nperms,2);
        res{iloc}.rand95accuracy(1:2) = NaN;%prctile(bootAccuracy, 95);
    end
end