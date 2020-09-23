function [ res ] = plotROI_localizer_illusion(v, roiname, crossValidation,locThresh)%, title )
toPlot=1;
nperms=1000;

% load the data
disp(sprintf('(plotROI): loading the data'));

% get the stimvols from ev anal
concatGroupNum = viewGet(v,'groupNum','Concatenation');
avgGroupNum = viewGet(v,'groupNum','Averages');
avgNscans = viewGet(v,'nscans',avgGroupNum);
for iscan=1:2%avgNscans
    s=viewGet(v,'stimfile',iscan,avgGroupNum);
    stimStr = strfind(s{1}.task{1}.taskFilename,'Still');
    eyeStr = strfind(s{1}.task{1}.taskFilename,'OccM');
    if ~isempty(stimStr)%strcmp('driftIllusion_Still_Control.m',s{1}.task{1}.taskFilename)
        stimlocScan = iscan;
    elseif ~isempty(eyeStr)%strcmp('driftIllusion_OccM_Control.m',s{1}.task{1}.taskFilename)
        eyelocScan = iscan;
    else
        keyboard
    end
end
%     matchScan = viewGet(v,'curScan');
matchScan=1;
stimlocCorrData = loadROIcoranalMatching(v, roiname, stimlocScan, avgGroupNum, matchScan, concatGroupNum);
eyelocCorrData = loadROIcoranalMatching(v, roiname, eyelocScan, avgGroupNum, matchScan, concatGroupNum);

v = viewSet(v,'curGroup','Concatenation');
rois = loadROITSeries(v, roiname,1,concatGroupNum,'keepNAN=true');
roiNans = isnan(mean(rois.tSeries,2));
rois.tSeries = rois.tSeries(~roiNans,:);
rois.scanCoords = rois.scanCoords(:,~roiNans);
rois.n = size(rois.tSeries,1);

d = viewGet(v, 'd');
if isempty(d)
    v=loadAnalysis(v,'erAnal/erAnal.mat');
    d = viewGet(v, 'd');
end
if isempty(d)
    v=loadAnalysis(v,'erAnal/erAnal_24sec.mat');
    d = viewGet(v, 'd');
end
if isempty(d)
    v=loadAnalysis(v,'erAnal/erAnal_24.mat');
    d = viewGet(v, 'd');
end
if isempty(d)
    keyboard
end
[stimvol stimNames var] = getStimvol(v,{'groupNum=[1 3]','groupNum=2'},'taskNum=1','phaseNum=1');
d.stimvol = stimvol;
d = makescm(d, 16, 0, stimvol);

d = getr2timecourse(mean(rois.tSeries),d.nhdr,d.hdrlen,d.scm);

t = 0:1.5:24-1.5;
Colors = [[0/255 171/255 75/255]; [0/255 126/255 200/255]];


if toPlot
    h = figure(1);
    clf
    hold on
    %plot it so that the no illusion (blue) is under the illusion (green)
    i=2;
    myerrorbar(t, d.ehdr(i,:)', 'yError', d.ehdrste(i,:)', 'Color', Colors(i,:), 'MarkerSize', 20, 'LineWidth', 3.5);
    i=1;
    myerrorbar(t, d.ehdr(i,:)', 'yError', d.ehdrste(i,:)', 'Color', Colors(i,:), 'MarkerSize', 20, 'LineWidth', 3.5);
    
    % for i=1:size(d.ehdr,1)
    %     myerrorbar(t, d.ehdr(i,:)', 'yError', d.ehdrste(i,:)', 'Color', Colors(i,:), 'MarkerSize', 20, 'LineWidth', 3.5)
    %     hold on
    % end
    h1=plot(0,d.ehdr(1,1), 'Color', Colors(1,:), 'MarkerSize', 1);
    h2=plot(0,d.ehdr(2,1), 'Color', Colors(2,:), 'MarkerSize', 1);
    % legend([h1 h2],{'Illusion','No Illusion'});
    % leg=legend;
    % leg.FontSize = 19;
    % hl = findobj(gca,'type','line');
    % set(hl,'LineWidth',2.5);
    axis square; axis tight
    YLabel = sprintf('fMRI response \n(%%change image intensity)');
    drawPublishAxis('xTick', [0 24], 'yTick', [-2.5 0 2.5], 'xLabel', 'Time (sec)', ...
        'yLabel',  YLabel,  ...
        'figSize', [17.6], 'labelFontSize', 29);
    legend off
    h.Position=[5   5   22.0839   17.6389];
    % filename = sprintf('Figures/tsPlot_13v2_%s.pdf', roiname);
    % savepdf(h, filename);
end

%% compute betas for three conditions using averaged hdr (computed from er-anal, above)
% grab the event-related d structure from a running mrLoadRet
dd = viewGet(v, 'd');

% construct a stimulus-convolution matrix with 1 time point (indicated
% trial onsets
scm = makescm(v, 1, 0, dd.stimvol);

% compute average response (from above), scale by norm
hdr = mean(d.ehdr);
%hdr = d.ehdr(1,:);
hdr = hdr-hdr(1);
hdr = hdr/norm(hdr);


% convolve scm with hemodynamic response
scm = conv2(scm, hdr');
scm = scm(1:size(dd.scm,1),:);

% compute betas
nhdr = 3;
hdrlen = 1;
dd2 = getr2timecourse(mean(rois.tSeries),nhdr,hdrlen,scm(1:size(rois.tSeries,2),:));

%% compute betas for individual trials
% construct an scm with individual trials

trialscm = [];
nTrials = length(dd.stimvol{1});
for iCond=1:3
    for iTrial=1:nTrials
        temp = dd;
        temp.stimvol{iCond} = temp.stimvol{iCond}(iTrial);
        foo = makescm(v, 1, 0, temp.stimvol);
        trialscm = cat(2, trialscm, foo(:,iCond));
    end
end


% trial-wise-scm, convolved with cannonical hdr
trialscmhdr = conv2(trialscm,hdr');
trialscmhdr = trialscmhdr(1:size(rois.tSeries,2),:);

trialReg = getr2timecourse(rois.tSeries, size(trialscmhdr,2), hdrlen, trialscmhdr);

% now do classification!
% extract the betas for condition 'l' and condition 'r'
% keyboard

temp = trialReg.r2>prctile(trialReg.r2,50) & stimlocCorrData{1}.co(~roiNans)>locThresh & eyelocCorrData{1}.co(~roiNans)>locThresh;
res{1}.sharedVox = sum(temp);
res{1}.totalVox = length(temp);

for iloc=1:2
    switch iloc
        case 1
            locCorrData = stimlocCorrData{1};
        case 2
            locCorrData = eyelocCorrData{1};
    end
    %     keyboard
    locCorrData.co = locCorrData.co(~roiNans);%data from concatenation excluded these voxels
    
    goodVox = trialReg.r2>prctile(trialReg.r2,50) & locCorrData.co>locThresh;

    res{iloc}.numVox = sum(goodVox);
    [roiname ' loc ' num2str(iloc) ': using ' num2str(res{iloc}.numVox) ' out of ' num2str(res{1}.totalVox) ' voxels']
    if sum(goodVox)>0
        
%         betas = trialReg.ehdr(goodVox,[1:nTrials (nTrials*2)+1:nTrials*3])'; %this is for 1 vs 3
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
        
        % Find how many scans there are by the name of the Concatenation group
        concatInfo = viewGet(v, 'concatInfo');
        nScans = concatInfo.n;
        
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

