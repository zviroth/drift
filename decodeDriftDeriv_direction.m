close all
clear all
nperms = 1000;
expNum=2;
crossValidation = 'trial';%'trial','run'
tic
r2thresh = 50;%r2 pbercentile
locThresh = 0.2;
thrsh = num2str(locThresh,'%.2f');
locStr = thrsh([1 3:4]);

saveFolder = '~/noah/';
switch expNum
    case 1
        expName = '3conds';
    case 2
        expName = 'attn';
end

load([saveFolder 'driftTS_deriv_' expName '.mat'],'ROIs','subNames','dataFolder',...
    'stimLoc','eyeLoc','driftBetas','driftBetasDeriv','subRoiNans','subNumTrials','subNumScans');
numSubs = length(subNames);
numRois = length(ROIs);



for isub=1:numSubs
    nTrials = subNumTrials(isub);
    for iroi=1:numRois
        roiname = ROIs{iroi};
        % res{1}.totalVox = length(temp);
        roiNans = subRoiNans{isub,iroi};
%         trialReg = driftBetas{isub,iroi};
        trialReg = driftBetasDeriv{isub,iroi};
        if ~isempty(trialReg)
            temp = trialReg.r2>prctile(trialReg.r2,50) & stimLoc{isub,iroi}.co(~roiNans)>locThresh & eyeLoc{isub,iroi}.co(~roiNans)>locThresh;
            locSharedVox(isub,iroi) = sum(temp);
            totalVox(isub,iroi) = length(temp);
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
                ['sub' num2str(isub) ' ' roiname ' loc ' num2str(iloc) ': using ' num2str(locNumVox(isub,iroi,iloc)) ' out of ' num2str(totalVox(isub,iroi)) ' voxels']
                if sum(goodVox)>0
                    
                    betas = trialReg.ehdr(goodVox,[1:nTrials (nTrials*2)+1:nTrials*3])'; %this is for 1 vs 3
                    betas = [betas trialReg.ehdr(goodVox,[nTrials*3+1:nTrials*4 (nTrials*5)+1:nTrials*6])']; %this is for 1 vs 3
%                             for iLR = 1:2
%                                 switch iLR
%                                     case 1
%                                         betas = trialReg.ehdr(:,[(nTrials)+1:nTrials*2 (nTrials*2)+1:nTrials*3])'; %this is for 2 vs 3
%                                     case 2
%                                         betas = trialReg.ehdr(:,[1:nTrials (nTrials)+1:nTrials*2])'; %this is for 1 vs 2
%                                 end
                    %betas = trialReg.ehdr(:,[(nTrials)+1:nTrials*2 (nTrials*2)+1:nTrials*3])'; %this is for 2 vs 3
                    %betas = trialReg.ehdr(:,[1:nTrials (nTrials)+1:nTrials*2])'; %this is for 1 vs 2
                    betas = zscore(betas);
                    
                    % create a grouping variable
                    group = cat(1, repmat('l', nTrials, 1), repmat('r', nTrials, 1));
                    
                    if strcmp(crossValidation, 'trial')
                        CVO = cvpartition(group, 'KFold', nTrials*2);
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
                    %                 nBoot = 1000;
                    bootAccuracy = zeros(nperms, 1);
                    
                    if strcmp(crossValidation, 'trial')
                        for iBoot=1:nperms
                            gNull = randLabels(:,iBoot);
                            for i = 1:CVO.NumTestSets
                                trIdx = CVO.training(i);
                                teIdx = CVO.test(i);
                                ytest = classify(betas(teIdx,:),betas(trIdx,:),gNull(trIdx), 'diagLinear');
                                err(i) = sum(~strcmp(ytest,gNull(teIdx)));
                            end
                            bootAccuracy(iBoot) = 1 - sum(err)/sum(CVO.TestSize);
                        end
                    elseif strcmp(crossValidation, 'run')
                        for iBoot=1:nperms
                            gNull = randLabels(:,iBoot);
                            clear ytest;
                            for iRun = 1:nRuns
                                trIdx = runVec ~= iRun;
                                teIdx = runVec == iRun;
                                ytest(teIdx) = classify(betas(teIdx,:),betas(trIdx,:),gNull(trIdx), 'diagLinear');
                            end
                            bootAccuracy(iBoot) = sum(ytest==group') / length(group);
                        end
                    else
                        keyboard
                    end

                    disp(sprintf('(plotROI): Bootstrapped threshold is: %f', prctile(bootAccuracy, 95)));
                    
                    subAcc(isub,iroi,iloc) = accuracy;
                    subRandAcc(isub,iroi,iloc,:) = bootAccuracy;
                    subRand95Acc(isub,iroi,iloc) = prctile(bootAccuracy, 95);
                    %         end
                else
                    'ZERO GOOD VOXELS!!'
                    subAcc(isub,iroi,iloc) = NaN;
                    subRandAcc(isub,iroi,iloc,:) = NaN(nperms,1);
                    subRand95Acc(isub,iroi,iloc) = NaN;%prctile(bootAccuracy, 95);
                end
            end
        else
            'ZERO GOOD VOXELS!!'
            subAcc(isub,iroi,1:2) = NaN;
            subRandAcc(isub,iroi,1:2,:) = NaN(2,nperms);
            subRand95Acc(isub,iroi,1:2) = NaN;%prctile(bootAccuracy, 95);
        end
    end
end

varinfo=whos('subRandAcc');
saveopt='';
if varinfo.bytes >= 2^31
    saveopt='-v7.3';
end

save([saveFolder 'decodeDriftDeriv_direction_' expName '_' crossValidation '_thresh' locStr '.mat'],'crossValidation',...
    'r2thresh','locThresh','ROIs','subNames','dataFolder',...
    'subRoiNans','subNumTrials','subNumScans',...
    'subAcc','subRandAcc','subRand95Acc','nperms','locSharedVox','totalVox');
toc