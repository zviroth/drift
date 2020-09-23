close all
clear all
nperms = 1000;

expNum=4;
crossValidation = 'run';%'trial','run'
tic
r2thresh = 50;%r2 pbercentile
% locThresh = 0.1;
% thrsh = num2str(locThresh,'%.2f');
% locStr = thrsh([1 3:4]);

saveFolder = '~/noah/';
switch expNum
    case 3
        nconds=4;
        expName = '4conds';
    case 4
        nconds=3;
        expName = '3tdata';
end

load([saveFolder 'driftTS_' expName '.mat'],'ROIs','subNames','dataFolder',...
    'driftBetas','subRoiNans','subNumTrials','subNumScans');

numSubs = length(subNames);
numRois = length(ROIs);



for isub=1:numSubs
    nTrials = subNumTrials(isub);
    for iroi=1:numRois
        roiname = ROIs{iroi};
        % res{1}.totalVox = length(temp);
        roiNans = subRoiNans{isub,iroi};
        trialReg = driftBetas{isub,iroi};
        
        if ~isempty(trialReg)
%             nTrials = size(trialReg.ehdr,2)/3;
            nTrials = (subNumTrials(isub,1) + subNumTrials(isub,3))/2;
            goodVox = trialReg.r2>prctile(trialReg.r2,50);
            totalVox(isub,iroi) = length(goodVox);
            numVox(isub,iroi) = sum(goodVox);
            ['sub' num2str(isub) ' ' roiname ': using ' num2str(numVox(isub,iroi)) ' voxels out of ' num2str(totalVox(isub,iroi)) ' voxels']
            if sum(goodVox)>0
                
                for idecode=1%:2%illusion or just local motion
%                     betas = trialReg.ehdr(goodVox,[1:nTrials (nTrials*2)+1:nTrials*3])'; %this is for 1 vs 3
                    betas = trialReg.ehdr(goodVox,[1:subNumTrials(isub,1) subNumTrials(isub,1)+subNumTrials(isub,2)+1:end])'; %this is for 1 vs 3

                    %                     if idecode==1
                    %                         betas = trialReg.ehdr(goodVox,[1:nTrials (nTrials)+1:nTrials*2])'; %this is for 1 vs 2
                    %                     else
                    %                         betas = trialReg.ehdr(goodVox,[nTrials*2+1:3*nTrials (3*nTrials)+1:nTrials*4])'; %this is for 3 vs 4
                    %                     end
                    betas = zscore(betas);
                    
                    % create a grouping variable
                    group = cat(1, repmat('l', subNumTrials(isub,1), 1), repmat('r', subNumTrials(isub,3), 1));
                    
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
                        
                        %                         clear accuracy; clear ytest;
                        %                         nRuns = length(group)/2/4;
                        %                         runVec = repmat(1:nRuns, 4, 1);
                        %                         runVec = runVec(:);
                        %                         runVec = cat(1, runVec, runVec);
                        %                         for iRun = 1:nRuns
                        %                             trIdx = runVec ~= iRun;
                        %                             teIdx = runVec == iRun;
                        %                             ytest(teIdx) = classify(betas(teIdx,:),betas(trIdx,:),group(trIdx), 'diagLinear');
                        %                         end
                        %                         accuracy = sum(ytest==group') / length(group);
                        %                         disp(sprintf('(plotROI): Classification accuracy: %f', accuracy));
                    else
                        keyboard
                    end
                    
                    nScans = subNumScans(isub);
                    
                    % trialsPerCond will be the trials of BOTH CONDITIONS
                    %                     trialsPerCond = length(group) / nScans;
                    trialsPerCond = 8;
                    groupByScan = cat(1, repmat('l', trialsPerCond/2, 1), repmat('r', trialsPerCond/2, 1));
                    
                    % We have to create 'length(nScans)' amount of randomized labels
                    clear randLabels
                    Cond1 = [];
                    Cond3 = [];
                    for iBoot = 1:1000
                        correctLength = 0;
                        while correctLength == 0
                            tempList = groupByScan(randperm(size(groupByScan,1)));
                            if length(Cond1)+4 <= subNumTrials(isub,1)
                                Cond1 = [Cond1; tempList(1:4)];
                            else
                                remaining = subNumTrials(isub,1)-length(Cond1);
                                if remaining == 1
                                    Cond1 = [Cond1; tempList(1)];
                                else
                                    Cond1 = [Cond1; tempList(1:subNumTrials(isub,1)-length(Cond1))];
                                end
                            end
                            if length(Cond3)+4 <= subNumTrials(isub,3)
                                Cond3 = [Cond3; tempList(5:end)];
                            else
                                remaining = subNumTrials(isub,3)-length(Cond3);
                                if remaining == 1
                                    Cond3 = [Cond3; tempList(end)];
                                else
                                    Cond3 = [Cond3; tempList(randperm(subNumTrials(isub,3)-length(Cond3)-1))];
                                end
                            end
                            
                            if length(Cond1) == subNumTrials(isub,1) && length(Cond3) == subNumTrials(isub,3)
                                randLabels(1:length(group),iBoot) = [Cond1; Cond3];
                                Cond1 = [];
                                Cond3 = [];
                                correctLength = 1;
                            end
                        end
                    end
%                     for iBoot = 1:nperms
%                         %                         for iScan = 1:nScans
%                         %                             tempList = groupByScan(randperm(size(groupByScan,1)));
%                         %                             idx1 = (iScan*4)-3;
%                         %                             idx2 = iScan*4;
%                         %                             randLabels(idx1:idx2,iBoot) = tempList(1:4);
%                         %                             randLabels(idx1+length(group)/2 : idx2+length(group)/2,iBoot) = tempList(5:8);
%                         %                         end
%                         for iScan = 1:nScans
%                             tempList = groupByScan(randperm(size(groupByScan,1)));
%                             idx1 = (iScan*trialsPerCond/2)-2;
%                             idx2 = iScan*trialsPerCond/2;
%                             randLabels(idx1:idx2,iBoot) = tempList(1:trialsPerCond/2);
%                             randLabels(idx1+length(group)/2 : idx2+length(group)/2,iBoot) = tempList((trialsPerCond/2)+1:end);
%                         end
%                     end
                    
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
                    
                    subAcc(isub,iroi,idecode) = accuracy;
                    subRandAcc(isub,iroi,idecode,:) = bootAccuracy;
                    subRand95Acc(isub,iroi,idecode) = prctile(bootAccuracy, 95);
                end
            else
                'ZERO GOOD VOXELS!!'
                subAcc(isub,iroi,1:2) = NaN;
                subRandAcc(isub,iroi,1:2,:) = NaN(nperms,1);
                subRand95Acc(isub,iroi,1:2) = NaN;%prctile(bootAccuracy, 95);
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
%%
save([saveFolder 'decodeDrift_direction_' expName '_' crossValidation '.mat'],'crossValidation',...
    'r2thresh','ROIs','subNames','dataFolder',...
    'subRoiNans','subNumTrials','subNumScans',...
    'subAcc','subRandAcc','subRand95Acc','nperms','numVox','totalVox');
toc