function [ res ] = plotROI_localizer_roiList(v, ROIs, crossValidation,locThresh)%, title )
%UNTITLED After specifying the name of an ROI, this program creates the
%(illusion vs no illusion) bar graph that was used for NIH Post Bac Poster
%Day 2018, and can be used for similar analysis between conditions of fMRI
%data
% locThresh=0.3;

%   Detailed explanation goes here
%roiname = 'lV1stimloc';
%roiname = 'lMTplus';
nperms=1000;
if ~exist('Figures', 'dir')
    mkdir('Figures')
end
filename = sprintf('Anal/roiAnal_%s.mat', roiname);
filename = 'randomfoo123123123.mat';
if exist(filename, 'file') == 2
    keyboard
    disp(sprintf('(plotROI): Data already created for %s. Loading...', roiname))
    disp('  ');
    load(filename)
    disp('Data loaded. Analyzing...');
    disp('  ');
    trialReg = roiAnal.trialReg;
    nTrials = size(roiAnal.trialReg.ehdr,2)/3;
    %clear roiAnal;
else
    
    % load the data
    disp(sprintf('(plotROI): loading the data'));
    
    
    
    % get the stimvols from ev anal
    concatGroupNum = viewGet(v,'groupNum','Concatenation');
    avgGroupNum = viewGet(v,'groupNum','Averages');
    avgNscans = viewGet(v,'nscans',avgGroupNum);
    for iscan=1:avgNscans
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
        keyboard
    end
    [stimvol stimNames var] = getStimvol(v,{'groupNum=[1 3]','groupNum=2'},'taskNum=1','phaseNum=1');
    d.stimvol = stimvol;
    d = makescm(d, 16, 0, stimvol);
    
    d = getr2timecourse(mean(rois.tSeries),d.nhdr,d.hdrlen,d.scm);
    
    t = 0:1.5:24-1.5;
    Colors = [[0/255 171/255 75/255]; [0/255 126/255 200/255]];
    
    
    
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
    filename = sprintf('Figures/tsPlot_13v2_%s.pdf', roiname);
    savepdf(h, filename);
    
    
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
    
    % % alternative method is to grab it from the raw data
    % % load the time series from the roi
    % ts = loadROITSeries(getMLRView, 'roisDeep');
    % % average across voxels
    % hrf = mean(ts.tSeries);
    % % average across all trials
    % hrf = mean(reshape(hrf, 16, length(hrf)/16), 2);
    % % remove mean
    % hrf = hrf - hrf(1); %mean(hrf);
    % % normalize
    % hrf = hrf / norm(hrf');
    % hdr = hrf';
    
    
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
end
% now do classification!
% extract the betas for condition 'l' and condition 'r'
% keyboard

temp = trialReg.r2>prctile(trialReg.r2,50) & stimlocCorrData{1}.co>locThresh & eyelocCorrData{1}.co>locThresh;
res{1}.sharedVox = sum(temp);

for iloc=1:2
    switch iloc
        case 1
            locCorrData = stimlocCorrData{1};
        case 2
            locCorrData = eyelocCorrData{1};
    end
    %     keyboard
    locCorrData.co = locCorrData.co(~roiNans);%data from concatenation excluded these voxels
    %     locCorrData.co(roiNans) = [];%data from concatenation excluded these voxels
    
    %     stimlocCorrData{1}.co
    %     goodVox = trialReg.r2>prctile(trialReg.r2,50);
    goodVox = trialReg.r2>prctile(trialReg.r2,50) & locCorrData.co>locThresh;
    keyboard
    res{iloc}.numVox = sum(goodVox);
    if sum(goodVox)>0
        
        betas = trialReg.ehdr(goodVox,[1:nTrials (nTrials*2)+1:nTrials*3])'; %this is for 1 vs 3
        
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
        
        disp(sprintf('(plotROI): Computing bootstrap for L vs R'));
        disp(sprintf('(plotROI): (the new, Martin-version of bootstrapping!)'))
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
        h = figure(2);
        clf
        hist(bootAccuracy);
        vline(accuracy, 'k');
        vline(prctile(bootAccuracy, 95), 'r');
        title('black = accuracy, red = 95% of null distribution')
        filename = sprintf('Figures/nullDistribution_%s_%swise.pdf', roiname, crossValidation);
        savepdf(h, filename)
        
        res{iloc}.accuracy = accuracy;
        res{iloc}.randAccuracy = bootAccuracy;
        res{iloc}.rand95accuracy = prctile(bootAccuracy, 95);
    else
        'ZERO GOOD VOXELS!!'
        res{iloc}.accuracy = NaN;
        res{iloc}.randAccuracy = NaN(nperms,1);
        res{iloc}.rand95accuracy = NaN;%prctile(bootAccuracy, 95);
    end
end

% if ~exist('Anal', 'dir')
%     mkdir('Anal')
% end
% filename = sprintf('Anal/roiAnal_%s_%swise.mat', roiname, crossValidation);
% filename = 'tempfoo.mat'; %needs to be saved over while i'm debugging
% if exist(filename, 'file') == 2
%     roiAnal.bootAccuracy = bootAccuracy; %save over the incorrect way
%     roiAnal.accuracy = accuracy;
%     disp(' ');
%     disp('Overwriting old bootAccuracy');
%     disp('Overwriting old bootAccuracy.');
%     disp('Overwriting old bootAccuracy..');
%     disp('Overwriting old bootAccuracy... complete');
% elseif exist(filename, 'file') ~= 2 && strcmp(crossValidation, 'run')
%     roiAnal.accuracy = accuracy;
%     roiAnal.bootAccuracy = bootAccuracy;
%     roiAnal.trialReg = trialReg;
%     %not sure if we can always grab these next 2 lines ...?
%         roiAnal.d = d;
%     roiAnal.dd2 = dd2;
% else
%     roiAnal.accuracy = accuracy;
%     roiAnal.bootAccuracy = bootAccuracy;
%     roiAnal.d = d;
%     roiAnal.dd2 = dd2;
%     roiAnal.trialReg = trialReg;
%     filename = sprintf('Anal/roiAnal_%s', roiname);
% end
% varinfo=whos('roiAnal');
% saveopt='';
% if varinfo.bytes >= 2^31
%     saveopt='-v7.3';
% end
% filename = sprintf('Anal/roiAnal_%s_%swise.mat', roiname, crossValidation);
%
% save(filename,'roiAnal',saveopt);
