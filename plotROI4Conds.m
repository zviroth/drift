function [ h ] = plotROI4Conds(v, roiname, Conds, crossValidation)%, title )
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
filename = sprintf('Anal/roiAnal_%s_%d_%d.mat', roiname, Conds(1), Conds(2));
%filename = 'foo23423324.mat';
if exist(filename, 'file') == 2
    disp(sprintf('(plotROI4Conds): Data already created for %s. Loading...', roiname))
    disp('  ');
    load(filename)
    disp('Data loaded. Analyzing...');
    disp('  ');
    trialReg = roiAnal.trialReg;
    nTrials = size(roiAnal.trialReg.ehdr,2)/4; %this is divided by 4 here, but divided by 3 in plotROI. Confirmed! Don't get confused! 
    %clear roiAnal;
else
% %     % For the data collected in the 'Attn_4_Conds', we need to figure out
% %     % which scans were the Attention 4 condition task, using:
% %     v = viewSet(v, 'curGroup', 'Raw');
% %     
% %     % figure out which scans to run on
% %     % assumed the 'description' field has been set with Run1 to RunN
% %     nScans = viewGet(v, 'nscans');
% %     whichScans = [];
% %     for iScan=1:nScans
% %         if ~isempty(strfind(viewGet(v, 'description', iScan), 'Attention'))
% %             whichScans = [whichScans iScan];
% %         end
% %     end
% %     
% %     % Go back to Concatenation group
% %     v = viewSet(v, 'curGroup', 'Concatenation');
% %     disp(sprintf('(plotROI4Conds): Found %d scans', length(whichScans)))
    
    % load the data
    disp(sprintf('(plotROI4Conds): loading the data'));
    v = viewSet(v,'curGroup','Concatenation');
    rois = loadROITSeries(v, roiname);
    
% for iRun = 1:length(whichScans)
%   v = viewSet(v,'curScan',whichScans(iRun));
%   [stimvol stimNames var] = getStimvol(v,{groupNum1, groupNum2},'taskNum=1','phaseNum=1');
%   design{iRun} = zeros(140,3);
%   design{iRun}(stimvol{1},1) = 1;
%   design{iRun}(stimvol{2},2) = 1;
%   design{iRun}(stimvol{3},3) = 1;
%   design{iRun} = design{iRun}(junkFrames+1:end,:);
% end
    
    % get the stimvols from ev anal
    d = viewGet(v, 'd');
    if isempty(d)
%         v=loadAnalysis(v);
        v=loadAnalysis(v,'erAnal/erAnal.mat');
        d = viewGet(v, 'd');
    end
    if isempty(d)
        keyboard
    end
    groupNum1 = sprintf('groupNum = %d', Conds(1));
    groupNum2 = sprintf('groupNum = %d', Conds(2));
    [stimvol stimNames var] = getStimvol(v,{groupNum1, groupNum2},'taskNum=1','phaseNum=1');
    d.stimvol = stimvol;
    d = makescm(d, 16, 0, stimvol);
    
    d = getr2timecourse(mean(rois.tSeries),d.nhdr,d.hdrlen,d.scm);
    
    t = 0:1.5:24-1.5;
    Colors = [[0/255 171/255 75/255]; [0/255 126/255 200/255]];
    

    %% compute betas for 4 conditions using averaged hdr (computed from er-anal, above)
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
    nhdr = 4;
    hdrlen = 1;
    dd2 = getr2timecourse(mean(rois.tSeries),nhdr,hdrlen,scm(1:size(rois.tSeries,2),:));
    
    %% compute betas for individual trials
    % construct an scm with individual trials
    trialscm = [];
    nTrials = length(dd.stimvol{1});
    for iCond=1:4
        for iTrial=1:nTrials
            thisstimvol = zeros(length(rois.tSeries),1);
            thisstimvol(dd.stimvol{iCond}(iTrial)) = 1;
            trialscm = cat(2, trialscm, thisstimvol);
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

if isequal(Conds, [1 2])
    betas = trialReg.ehdr(:,[1:nTrials (nTrials)+1:nTrials*2])'; %this is for 1 vs 2
elseif isequal(Conds, [3 4])
    betas = trialReg.ehdr(:,[(nTrials*2)+1:nTrials*3 (nTrials*3)+1:nTrials*4])'; %this is for 3 vs 4
else
    keyboard
end
%betas = trialReg.ehdr(goodVox,[1:nTrials (nTrials*2)+1:nTrials*3])'; %this is for 1 vs 3
%betas = trialReg.ehdr(:,[(nTrials)+1:nTrials*2 (nTrials*2)+1:nTrials*3])'; %this is for 2 vs 3
%betas = trialReg.ehdr(:,[1:nTrials (nTrials)+1:nTrials*2])'; %this is for 1 vs 2
betas = zscore(betas);

% create a grouping variable
group = cat(1, repmat('l', nTrials, 1), repmat('r', nTrials, 1));
if isequal(Conds, [1 2])
    disp(sprintf('(plotROI4Conds): Computing accuracy for L vs R'));
elseif isequal(Conds, [3 4])
    disp(sprintf('(plotROI4Conds): Computing accuracy for LOCAL L vs R'));
end
% create partitioning scheme
CVO = cvpartition(group, 'KFold', nTrials*2);
%CVO = cvpartition(group, 'LeaveOut', 1);
%CVO = cvpartition(group, 'HoldOut', 0.5);

disp(sprintf('crossValidation method = %s', crossValidation));
if strcmp(crossValidation, 'trial')
err = zeros(CVO.NumTestSets,1);
for i = 1:CVO.NumTestSets
    trIdx = CVO.training(i);
    teIdx = CVO.test(i);
    ytest = classify(betas(teIdx,:),betas(trIdx,:),group(trIdx), 'diagLinear');
    err(i) = sum(~strcmp(ytest,group(teIdx)));
end
accuracy = 1 - sum(err)/sum(CVO.TestSize);
disp(sprintf('(plotROI4Conds): Classification accuracy: %f', accuracy));
elseif strcmp(crossValidation, 'run')
    clear accuracy; clear ytest;
    nRuns = length(group)/2/3;
    runVec = repmat(1:nRuns, 3, 1);
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
    disp('you haven''t given a proper crossValidation arguement!')
    keyboard
end
    


%
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
trialsPerCond = length(group) / nScans;
groupByScan = cat(1, repmat('l', trialsPerCond/2, 1), repmat('r', trialsPerCond/2, 1));

% We have to create 'length(nScans)' amount of randomized labels
clear randLabels
for iBoot = 1:1000
    for iScan = 1:nScans
        tempList = groupByScan(randperm(size(groupByScan,1)));
        idx1 = (iScan*trialsPerCond/2)-2;
        idx2 = iScan*trialsPerCond/2;
        randLabels(idx1:idx2,iBoot) = tempList(1:trialsPerCond/2);
        randLabels(idx1+length(group)/2 : idx2+length(group)/2,iBoot) = tempList((trialsPerCond/2)+1:end);
    end
end

disp(sprintf('(plotROI4Conds): Computing bootstrap for L vs R'));
disp(sprintf('(plotROI4Conds): (the new, Martin-version of bootstrapping!)'))
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


close all
h = figure(2);
hist(bootAccuracy)
vline(accuracy, 'k')
vline(prctile(bootAccuracy, 95), 'r')
title('black = accuracy, red = 95% of null distribution')
filename = sprintf('Figures/nullDistribution_%s_%d_%d_%swise.pdf', roiname, Conds(1), Conds(2), crossValidation);
savepdf(h, filename)

if ~exist('Anal', 'dir')
    mkdir('Anal')
end


filename = sprintf('Anal/roiAnal_%s_%d_%d_%swise.mat', roiname, Conds(1), Conds(2), crossValidation);
%filename = 'foo234234.mat'; %make sure to save it no matter what :)
if exist(filename, 'file') == 2
    roiAnal.bootAccuracy = bootAccuracy; %save over the incorrect way
    disp(' ');
    disp('Overwriting old bootAccuracy');
    disp('Overwriting old bootAccuracy.');
    disp('Overwriting old bootAccuracy..');
    disp('Overwriting old bootAccuracy... complete');
elseif exist(filename, 'file') ~= 2 && strcmp(crossValidation, 'run')
    roiAnal.accuracy = accuracy;
    roiAnal.bootAccuracy = bootAccuracy;
    roiAnal.trialReg = trialReg;
else
    roiAnal.accuracy = accuracy;
    roiAnal.bootAccuracy = bootAccuracy;
    roiAnal.d = d;
    roiAnal.dd2 = dd2;
    roiAnal.trialReg = trialReg;
    %filename = sprintf('Anal/roiAnal_%s_%d_%d.mat', roiname, Conds(1), Conds(2));
end
varinfo=whos('roiAnal');
saveopt='';
if varinfo.bytes >= 2^31
    saveopt='-v7.3';
end
save(filename,'roiAnal',saveopt);

%save(filename, 'roiAnal')