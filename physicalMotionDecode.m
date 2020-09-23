mrQuit
close all
clear all
nperms=1000;
r2thresh=50;
crossValidation = 'trial';%'trial','run'
dataFolder = '/Volumes/Drift/PhysicalMotion/';
subNames = {'s010020191107'};
saveFolder = '~/noah/';
tic
if exist('v','var')
    deleteView(v);
end

ROIs = {'lV1','lV2','lV3','lhV4','lLO1','lLO2','lTO1','lTO2','lV3a','lV3b','lV01','lV02',...
    'rV1','rV2','rV3','rhV4','rLO1','rLO2','rTO1','rTO2','rV3a','rV3b','rV01','rV02'};
% ROIs = {'lV1','rV1'};
numRois = length(ROIs);
numSubs = length(subNames);
saveFolder = '~/noah/';
cd(dataFolder);

nhdr = 3;
hdrlen = 1;
matchScan=1;

for isub=1:numSubs
    cd(subNames{isub});
    
    v=newView;
    
    concatGroupNum = viewGet(v,'groupNum','Concatenation');
    v = viewSet(v,'curGroup','Concatenation');
    erFilename =  dir('Concatenation/erAnal/erAnal.mat');
    
    clear d
    d = viewGet(v, 'd');
    if isempty(d)
        v=loadAnalysis(v,erFilename.name,fullfile(pwd,'/Concatenation/erAnal/'));
        d = viewGet(v, 'd');
    end
    dd=d;
%     [stimvol stimNames var] = getStimvol(v,{'dir=[0.7854 2.3562]'},'taskNum=2','phaseNum=1');
    [stimvol stimNames var] = getStimvol(v,'dir','taskNum=2','phaseNum=1');
%     keyboard
%     [stimvol stimNames var] = getStimvol(v,{'dir=[1 2]'},'taskNum=1','phaseNum=1');
    
    d.stimvol = stimvol;
    d = makescm(d, 16, 0, stimvol);%design matrix convolved with a 16 TR canonical HRF, with 2 conditions
    scm = makescm(v, 1, 0, dd.stimvol);%design matrix with delta functions, with 2 conditions
    % Find how many scans there are by the name of the Concatenation group
    concatInfo = viewGet(v, 'concatInfo');
    subNumScans(isub) = concatInfo.n;
    
    % construct an scm with individual trials
    trialscm = [];
    nconds=2;
    for iCond=1:nconds
        nTrials = length(dd.stimvol{iCond});
        subNumTrials(isub,iCond) = nTrials;
        for iTrial=1:nTrials
            temp = dd;
            temp.stimvol{iCond} = temp.stimvol{iCond}(iTrial);
            foo = makescm(v, 1, 0, temp.stimvol);
            trialscm = cat(2, trialscm, foo(:,iCond));
        end
    end
    
    for iroi=1:numRois
        roiname = ROIs{iroi};
        %save concatenated drift + 2 localizer timeseries for this ROI
        % get the stimvols from ev anal
        
        v = viewSet(v,'curGroup','Concatenation');
        rois = loadROITSeries(v, roiname,1,concatGroupNum,'keepNAN=true');
        roiNans = isnan(mean(rois.tSeries,2));
        rois.tSeries = rois.tSeries(~roiNans,:);
        rois.scanCoords = rois.scanCoords(:,~roiNans);
        rois.n = size(rois.tSeries,1);
        
        if rois.n>0
            %average across all ROI voxels?
            d = getr2timecourse(mean(rois.tSeries),d.nhdr,d.hdrlen,d.scm);%gives 2 deconvolved responses, one for L+R, one for no illusion
            
            % compute average response (from above), scale by norm
            hdr = mean(d.ehdr);
            %hdr = d.ehdr(1,:);
            hdr = hdr-hdr(1);
            hdr = hdr/norm(hdr);
            
            % trial-wise-scm, convolved with cannonical hdr
            trialscmhdr = conv2(trialscm,hdr');
            trialscmhdr = trialscmhdr(1:size(rois.tSeries,2),:);
            %get betas
            %         trialReg = getr2timecourse(rois.tSeries, size(trialscmhdr,2), hdrlen, trialscmhdr);
            driftBetas{isub,iroi} = getr2timecourse(rois.tSeries, size(trialscmhdr,2), hdrlen, trialscmhdr);
            
        else
            driftBetas{isub,iroi} = [];
        end
        subRoiNans{isub,iroi} = roiNans;
        
    end
    deleteView(v);
    cd ..

    
    nTrials = subNumTrials(isub,1);
    for iroi=1:numRois
        roiname = ROIs{iroi};
        % res{1}.totalVox = length(temp);
        roiNans = subRoiNans{isub,iroi};
        trialReg = driftBetas{isub,iroi};
        
        if ~isempty(trialReg)
%             nTrials = size(trialReg.ehdr,2)/3;
%             nTrials = (subNumTrials(isub,1) + subNumTrials(isub,2))/2;
            goodVox = trialReg.r2>prctile(trialReg.r2,r2thresh);
            totalVox(isub,iroi) = length(goodVox);
            numVox(isub,iroi) = sum(goodVox);
            ['sub' num2str(isub) ' ' roiname ': using ' num2str(numVox(isub,iroi)) ' voxels out of ' num2str(totalVox(isub,iroi)) ' voxels']
            if sum(goodVox)>0
                
                for idecode=1%:2%illusion or just local motion
%                     betas = trialReg.ehdr(goodVox,[1:nTrials (nTrials*2)+1:nTrials*3])'; %this is for 1 vs 3
%                     betas = trialReg.ehdr(goodVox,[1:subNumTrials(isub,1) subNumTrials(isub,1)+subNumTrials(isub,2)+1:end])'; %this is for 1 vs 3
                    betas = trialReg.ehdr(goodVox,:)'; %this is for 1 vs 2

                    betas = zscore(betas);
                    
                    % create a grouping variable
                    group = cat(1, repmat('l', subNumTrials(isub,1), 1), repmat('r', subNumTrials(isub,2), 1));
                    
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
                    %                     trialsPerCond = length(group) / nScans;
                    trialsPerCond = 8;
                    
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

toc
save([saveFolder 'decodePhysicalDirection_' crossValidation '.mat'],'crossValidation',...
    'r2thresh','ROIs','subNames','dataFolder',...
    'subRoiNans','subNumTrials','subNumScans',...
    'subAcc','subRandAcc','subRand95Acc','nperms','totalVox');

