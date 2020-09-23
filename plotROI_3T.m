function [ h ] = plotROI_3T(v, roiname, crossValidation)%, title )
%UNTITLED After specifying the name of an ROI, this program creates the
%(illusion vs no illusion) bar graph that was used for NIH Post Bac Poster
%Day 2018, and can be used for similar analysis between conditions of fMRI
%data


%   Detailed explanation goes here
%roiname = 'lV1stimloc';
%roiname = 'lMTplus';
if ~exist('Figures', 'dir')
    mkdir('Figures')
end
filename = sprintf('Anal/roiAnal_%s.mat', roiname);
filename = 'randomfoo123123123.mat';
if exist(filename, 'file') == 2
    keyboard
    disp(sprintf('(plotROI_3T): Data already created for %s. Loading...', roiname))
    disp('  ');
    load(filename)
    disp('Data loaded. Analyzing...');
    disp('  ');
    trialReg = roiAnal.trialReg;
    nTrials = size(roiAnal.trialReg.ehdr,2)/3;
    %clear roiAnal;
else
    
    % load the data
    disp(sprintf('(plotROI_3T): loading the data'));
    rois = loadROITSeries(v, roiname);
    
    % get the stimvols from ev anal
    d = viewGet(v, 'd');
    if isempty(d)
        keyboard
    end
    [stimvol stimNames var] = getStimvol(v,{'groupNum=[1 3]','groupNum=2'},'taskNum=1','phaseNum=1');
    d.stimvol = stimvol;
    d = makescm(d, 12, 0, stimvol);
    
    d = getr2timecourse(mean(rois.tSeries),d.nhdr,d.hdrlen,d.scm);
    
    t = 0:2:24-2;
    Colors = [[0/255 171/255 75/255]; [0/255 126/255 200/255]];
    
    h = figure(1);
    
    clf
    hold on
    %plot it so that the no illusion (blue) is under the illusion (green)
    i=2;
    myerrorbar(t, d.ehdr(i,:)', 'yError', d.ehdrste(i,:)', 'Color', Colors(i,:), 'MarkerSize', 20, 'LineWidth', 3.5)
    i=1;
    myerrorbar(t, d.ehdr(i,:)', 'yError', d.ehdrste(i,:)', 'Color', Colors(i,:), 'MarkerSize', 20, 'LineWidth', 3.5)
    
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
    drawPublishAxis('xTick', [0 24], 'yTick', [-1.5 0 1.5], 'xLabel', 'Time (sec)', ...
        'yLabel',  YLabel,  ...
        'figSize', [17.6], 'labelFontSize', 29);
    h.Position=[63.3236   25.9997   22.0839   17.6389];
    filename = sprintf('Figures/tsPlot_13v2_%s.pdf', roiname);
    savepdf(h, filename)
    
    
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
        for iTrial=1:length(dd.stimvol{iCond})
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
goodVox = trialReg.r2>prctile(trialReg.r2,50);
betas = trialReg.ehdr(goodVox,[1:length(dd.stimvol{1}) length(dd.stimvol{1})+length(dd.stimvol{2})+1:end])'; %this is for 1 vs 3

%betas = trialReg.ehdr(:,[(nTrials)+1:nTrials*2 (nTrials*2)+1:nTrials*3])'; %this is for 2 vs 3
%betas = trialReg.ehdr(:,[1:nTrials (nTrials)+1:nTrials*2])'; %this is for 1 vs 2
betas = zscore(betas);

% create a grouping variable
%group = cat(1, repmat('l', nTrials, 1), repmat('r', nTrials, 1));
group = cat(1, repmat('l', length(dd.stimvol{1}),1), repmat('r', length(dd.stimvol{3}), 1));

disp(sprintf('(plotROI_3T): Computing accuracy for L vs R'));
% create partitioning scheme
CVO = cvpartition(group, 'KFold', length(group));

if strcmp(crossValidation, 'trial')
    
    
    for i = 1:CVO.NumTestSets
        trIdx = CVO.training(i);
        teIdx = CVO.test(i);
        ytest(teIdx) = classify(betas(teIdx,:),betas(trIdx,:),group(trIdx), 'diagLinear');
    end
    accuracy = sum(ytest==group') / length(group);
    disp(sprintf('(plotROI_3T): Classification accuracy: %f', accuracy));
    
    
    
    
    
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
    disp(sprintf('(plotROI_3T): Classification accuracy: %f', accuracy));
else
    keyboard
end


% Old method of computing permutation test that randomizes once per BOOT
% disp(sprintf('(plotROI_3T): Computing bootstrap for L vs R'));
% nBoot = 1000;
% bootAccuracy = zeros(nBoot, 1);
% for iBoot=1:nBoot;
%     gNull = group(randperm(size(group,1)));
%     for i = 1:CVO.NumTestSets
%         trIdx = CVO.training(i);
%         teIdx = CVO.test(i);
%         ytest = classify(betas(teIdx,:),betas(trIdx,:),gNull(trIdx), 'diagLinear');
%         err(i) = sum(~strcmp(ytest,gNull(teIdx)));
%     end
%     bootAccuracy(iBoot) = 1 - sum(err)/sum(CVO.TestSize);
%     disppercent(iBoot/1000);
% end
% disppercent(inf);
% disp(sprintf('Bootstrapped threshold is: %f', prctile(bootAccuracy, 95)));

% G = cvpartition(group, 'HoldOut', 0.5);
% accuracy = zeros(1000,1); bootAccuracy = zeros(1000,1);
% for iPerm=1:1000
%     CVO = repartition(G);
%     gNull = group(randperm(size(group,1)));
%     trIdx = CVO.training(i);
%     teIdx = CVO.test(i);
%     ytest = classify(betas(teIdx,:),betas(trIdx,:),group(trIdx), 'diagLinear');
%     accuracy(iPerm) = sum(ytest == group(teIdx)) / sum(teIdx);
%     ytest = classify(betas(teIdx,:),betas(trIdx,:),gNull(trIdx), 'diagLinear');
%     bootAccuracy(iPerm) = sum(ytest == group(teIdx)) / sum(teIdx);
% end
%
% disp(sprintf('Classification accuracy: %f', accuracy));

% New way of randomizing labels - only randomize ONCE per run, and have
% equal amounts of L and R WITHIN run

% Find how many scans there are by the name of the Concatenation group
concatInfo = viewGet(v, 'concatInfo');
nScans = concatInfo.n;

% trialsPerCond will be the trials of BOTH CONDITIONS
trialsPerCond = 8;
groupByScan = cat(1, repmat('l', trialsPerCond/2, 1), repmat('r', trialsPerCond/2, 1));

if length(dd.stimvol{1}) - length(dd.stimvol{3}) > 4
    keyboard
end

% We have to create 'length(nScans)' amount of randomized labels
clear randLabels
%randLabels = zeros(1,1000);
Cond1 = [];
Cond3 = [];
for iBoot = 1:1000
    correctLength = 0;
    while correctLength == 0
        tempList = groupByScan(randperm(size(groupByScan,1)));
        if length(Cond1)+4 <= length(dd.stimvol{1})
            Cond1 = [Cond1; tempList(1:4)];
        else
            remaining = length(dd.stimvol{1})-length(Cond1);
            if remaining == 1
                Cond1 = [Cond1; tempList(1)];
            else
                Cond1 = [Cond1; tempList(1:length(dd.stimvol{1})-length(Cond1))];
            end
        end
        if length(Cond3)+4 <= length(dd.stimvol{3})
            Cond3 = [Cond3; tempList(5:end)];
        else
            remaining = length(dd.stimvol{3})-length(Cond3);
            if remaining == 1
                Cond3 = [Cond3; tempList(end)];
            else
                Cond3 = [Cond3; tempList(randperm(length(dd.stimvol{3})-length(Cond3)-1))];
            end
        end
        
        if length(Cond1) == length(dd.stimvol{1}) && length(Cond3) == length(dd.stimvol{3})
            randLabels(1:length(group),iBoot) = [Cond1; Cond3];
            Cond1 = [];
            Cond3 = [];
            correctLength = 1;
        end
    end
    %     correctLength = 0;
    %     Counter = 0;
    %     while correctLength == 0
    %         thisLoopSize = size(randLabels,1);
    %         Counter = Counter + 1;
    %         tempList = groupByScan(randperm(size(groupByScan,1)));
    %         idx1 = (Counter*4)-3;
    %         idx2 = Counter*4;
    %
    %
    %
    %
    %         for iScan = 1:nScans
    %             tempList = groupByScan(randperm(size(groupByScan,1)));
    %             %         idx1 = (iScan*4)-3;
    %             %         idx2 = iScan*4;
    %             %         randLabels(idx1:idx2,iBoot) = tempList(1:4);
    %             %         randLabels(idx1+length(group)/2 : idx2+length(group)/2,iBoot) = tempList(5:8);
    %             randLabels(:,iBoot) = group(randperm(size(group,1)));
    %         end
end

disp(sprintf('(plotROI_3T): Computing bootstrap for L vs R'));
disp(sprintf('(plotROI_3T): (the new, Martin-version of bootstrapping!)'))
disp(sprintf('(plotROI_3T): Classification method = %swise', crossValidation));
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
disp(sprintf('(plotROI_3T): Bootstrapped threshold is: %f', prctile(bootAccuracy, 95)));

close all
h = figure(2);
hist(bootAccuracy)
vline(accuracy, 'k');
vline(prctile(bootAccuracy, 95), 'r');
title('black = accuracy, red = 95% of null distribution')
filename = sprintf('Figures/nullDistribution_%s_%swise.pdf', roiname, crossValidation);
savepdf(h, filename)


if ~exist('Anal', 'dir')
    mkdir('Anal')
end
filename = sprintf('Anal/roiAnal_%s_%swise.mat', roiname, crossValidation);
filename = 'tempfoo.mat'; %needs to be saved over while i'm debugging
if exist(filename, 'file') == 2
    roiAnal.bootAccuracy = bootAccuracy; %save over the incorrect way
    roiAnal.accuracy = accuracy;
    disp(' ');
    disp('Overwriting old bootAccuracy');
    disp('Overwriting old bootAccuracy.');
    disp('Overwriting old bootAccuracy..');
    disp('Overwriting old bootAccuracy... complete');
elseif exist(filename, 'file') ~= 2 && strcmp(crossValidation, 'run')
    roiAnal.accuracy = accuracy;
    roiAnal.bootAccuracy = bootAccuracy;
    roiAnal.trialReg = trialReg;
    %not sure if we can always grab these next 2 lines ...?
        roiAnal.d = d;
    roiAnal.dd2 = dd2;
else
    roiAnal.accuracy = accuracy;
    roiAnal.bootAccuracy = bootAccuracy;
    roiAnal.d = d;
    roiAnal.dd2 = dd2;
    roiAnal.trialReg = trialReg;
    filename = sprintf('Anal/roiAnal_%s', roiname);
end
varinfo=whos('roiAnal');
saveopt='';
if varinfo.bytes >= 2^31
    saveopt='-v7.3';
end
filename = sprintf('Anal/roiAnal_%s_%swise.mat', roiname, crossValidation);

save(filename,'roiAnal',saveopt);

%save(filename, 'roiAnal')
