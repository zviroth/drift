analyzeData = 0;
IllusionVsNoIllusion = 0;
plotLvsR = 0;
plotScatter = 0;
plotModulationIndex = 1;
decodeData = 0;
Analysis = 'indiv';
%Analysis = 'group';
%Analysis = 'sanity';
%Analysis = 'indiv+sanity';
%Analysis = 'Attn_4_Conds';
%Analysis = 'wang+benson';
%Analysis = '3T';

NumberofAnalyses = 4; % not using the 95th percentile anymore, but now it's 100 - (.05/NumberofAnalyses);
crossValidation = 'run';

[SubjDir] = driftSubjDir( 'ROIs');
Conds = [];
disp('  '); disp(sprintf('This analyses group = %s', Analysis)); disp('  ');
if strcmp(Analysis, 'indiv')
    ROInames = {'stimloc_lh_individual_decoding_mt', 'stimloc_lh_individual_decoding_v1', 'stimloc_lh_individual_decoding_ips'};
    %ROInames = {'stimloc_lh_individual_decoding_ips'};
%     illvsNoIll.mtMat = []; illvsNoIll.v1Mat = []; illvsNoIll.ipsMat = [];
%     DecodingLvsR.mt = []; DecodingLvsR.v1 = []; DecodingLvsR.ips = [];
elseif strcmp(Analysis, 'group')
    ROInames = {'stimloc_lh_group_decoding_mt', 'stimloc_lh_group_decoding_ips', 'stimloc_lh_group_decoding_sts', 'stimloc_lh_group_decoding_po', 'stimloc_lh_group_decoding_spl'};
    illvsNoIll.mtMat = []; illvsNoIll.stsMat = []; illvsNoIll.ipsMat = []; illvsNoIll.poMat = []; illvsNoIll.splMat = [];
    DecodingLvsR.mt = []; DecodingLvsR.sts = []; DecodingLvsR.ips = []; DecodingLvsR.po = []; DecodingLvsR.spl = [];
elseif strcmp(Analysis, 'sanity')
    ROInames = {'control_v1_low_eccen_indiv'};
    illvsNoIll.v1CtrlMat = [];
    DecodingLvsR.v1Ctrl = [];
elseif strcmp(Analysis, 'indiv+sanity')
    ROInames = {'stimloc_lh_individual_decoding_mt', 'stimloc_lh_individual_decoding_v1', ...
        'stimloc_lh_individual_decoding_ips','control_v1_low_eccen_indiv'};
    
    SubjDir = [SubjDir, driftSubjDir( 'Attn_4_Conds')]; % Add the 4 Condition people to the attention version of the experiment
elseif strcmp(Analysis, '3T')
    
    ROInames = {'x3T_mt', 'x3T_v1', 'x3T_v3a'};
    %ROInames = {'x3T_v1'};
    NumberofAnalyses = 3;
    [SubjDir] = driftSubjDir( '3T');
elseif strcmp(Analysis, 'wang+benson')
    ROInames = {'benson_rh_stimlocfromlh_v1', 'lh_wang2015atlas_IPS0', ...
        'lh_wang2015atlas_IPS1','lh_wang2015atlas_IPS2', ...
        'lh_wang2015atlas_IPS3', 'lh_wang2015atlas_IPS4', ...
        'lh_wang2015atlas_IPS5'};
elseif strcmp(Analysis, 'Attn_4_Conds')
    [SubjDir] = driftSubjDir( 'Attn_4_Conds');
    ROInames = {'avg_wangAtlas_4Conds_v1', 'avg_wangAtlas_4Conds_mt', 'avg_wangAtlas_4Conds_v3a'};
    %SubjDir{4} = '';
    %SubjDir = SubjDir(5);
    Conds = [1 2];
    %Conds = [3 4];
else
    keyboard
end


% illvsNoIll.mtMat = []; illvsNoIll.v1Mat = []; illvsNoIll.ipsMat = [];
% DecodingLvsR.mt.accuracy = []; DecodingLvsR.v1.accuracy = []; DecodingLvsR.ips.accuracy = [];
% DecodingLvsR.mt.allEhdr = []; DecodingLvsR.v1.allEhdr = [];
% DecodingLvsR.ips.allEhdr = [];
% DecodingLvsR.mt.r2 = []; DecodingLvsR.v1.r2 = [];
% DecodingLvsR.ips.r2 = [];
Colors = [[0/255 171/255 75/255]; [0/255 126/255 200/255]; ...
    [127/255 127/255 127/255] ; [161/255 106/255 173/255]; ...
    [54/255 194/255 214/255]; [127/255 127/255 127/255];];

SubjCounter = 0;

%clean up the ROInames
ROInamesShort = {};
for i = 1:length(ROInames)
    temp = strsplit(ROInames{i}, '_');
    ROInamesShort{i,1} = temp{end};
end

clear roiData
for iSubj = 1:length(SubjDir)
    if ~isempty(SubjDir{iSubj})
        cd(SubjDir{iSubj})
        disp(sprintf('Analyzing Subj %d: %s', iSubj, getLastDir(pwd)));
        
        if analyzeData == 1
            v = newView;
            v = viewSet(v, 'curGroup', 'Concatenation');
            erAnalOptions = {'erAnal_24sec', 'erAnal24', 'erAnal_24', 'erAnal'};
            
            for iAnal = 1:length(erAnalOptions)
                if exist(sprintf('%s/Concatenation/erAnal/%s.mat', SubjDir{iSubj}, erAnalOptions{iAnal}))
                    disp('-----');
                    disp(sprintf('Loading analyses %s for %s', erAnalOptions{iAnal}, getLastDir(pwd)));
                    disp('-----');
                    disp('  '); %line break
                    thisfileloc = sprintf('erAnal/%s', erAnalOptions{iAnal});
                    v = loadAnalysis(v, thisfileloc);
                    break
                elseif iAnal == length(erAnalOptions)
                    keyboard
                    %couldn't find any analyses
                end
            end
            %v = loadAnalysis(v, 'erAnal');
            for iROI = 1:length(ROInames)
                if exist(sprintf('ROIs/%s.mat', ROInames{iROI}))
                    if strcmp(Analysis, '3T')
                        %we need to use a different function as the formation of 3T data sucks! (and is different than 7T)
                        plotROI_3T(v, ROInames{iROI}, crossValidation)
                    elseif ~strcmp(Analysis, 'Attn_4_Conds')
                        plotROI(v, ROInames{iROI}, crossValidation)
                    elseif strcmp(Analysis, 'Attn_4_Conds')
                        plotROI4Conds(v, ROInames{iROI}, [1 2], crossValidation) %Conditions 1 vs 2 (illusion left vs illusion right)
                        plotROI4Conds(v, ROInames{iROI}, [3 4], crossValidation) %Conditions 3 vs 4 (local left vs local right)
                    end
                else
                    disp(sprintf('%s does not exist!', ROInames{iROI}))
                end
            end
            deleteView(v)
            close all
        elseif IllusionVsNoIllusion == 1 || plotLvsR == 1 || plotScatter == 1 || decodeData == 1 || plotModulationIndex == 1
            cd('Anal');
            SubjCounter = SubjCounter + 1;
            
            for iROI = 1:length(ROInames)
                if ~strcmp(Analysis, 'Attn_4_Conds')
                    if isempty(crossValidation)
                        filename = sprintf('roiAnal_%s.mat', ROInames{iROI});
                    else
                        filename = sprintf('roiAnal_%s_%swise.mat', ROInames{iROI}, crossValidation);
                    end
                else
                    filename = sprintf('roiAnal_%s_%d_%d_%swise.mat', ROInames{iROI}, 1, 2, crossValidation);
                    filename_34 = sprintf('roiAnal_%s_%d_%d_%swise.mat', ROInames{iROI}, 3, 4, crossValidation);
                end
                
                if exist(filename) == 2
                    load(filename)
                    roiData{SubjCounter,iROI} = roiAnal;
                    if strcmp(Analysis, 'Attn_4_Conds')
                        load(filename_34)
                        roiData{SubjCounter,iROI+length(ROInames)} = roiAnal;
                    end
                else
                    if strcmp(ROInamesShort{iROI}, 'ips')
                        ROInamesShort{iROI} = 'v3a';
                    end
                    
                    % Try to find an 'equivalent' ROI (i.e. MT+ defined by
                    % stim localizer, vs MT+ defined by wangAtlas)
                    thisfile = dir(sprintf('%sAnal/roiAnal_avg_wangAtlas*%s_1_2.mat', SubjDir{iSubj}, ROInamesShort{iROI}));
                    if size(thisfile,1) > 1 || isempty(thisfile)
                        disp(sprintf('none loaded for %s + %s', Analysis, ROInames{iROI}))
                    else
                        temp = load(fullfile(SubjDir{iSubj}, 'Anal', thisfile.name));
                        roiData{SubjCounter,iROI} = temp.roiAnal;
                        disp(sprintf('For subject %s, used %s_%s instead of %s', getLastDir(SubjDir{iSubj}), 'roiAnal_avg_wangAtlas',ROInamesShort{iROI},  ROInames{iROI}))
                    end
                    
                    %clean up the ROInames
                    ROInamesShort = {};
                    for i = 1:length(ROInames)
                        temp = strsplit(ROInames{i}, '_');
                        ROInamesShort{i,1} = temp{end};
                    end
                    
                end
                
            end
        else
            keyboard
        end
    end
end
if decodeData ~= 1
    keyboard
end
if ~strcmp(Analysis, 'Attn_4_Conds')
    Decision = strfind(SubjDir, 'Decision')';
    Decision = find(~cellfun(@isempty,Decision));
    Attention = strfind(SubjDir, 'Att')';
    Attention = find(~cellfun(@isempty,Attention));
end

if decodeData == 1
    % For each ROI, find the following information:
    % 1) Get 50% and 95% prctile of combined null distribution (Martin's method)
    % clear variables
    bootMat = []; accMat = [];
    disp('  '); disp('Compiling average accuracy matrix (accMat) and bootstrap matrix (bootMat)');
    if strcmp(Analysis, 'Attn_4_Conds') && length(ROInames)*2~=size(roiData,2)
        keyboard % the for loop below needs to be changed back to uncommented
    end
    
    %for iROI = 1:length(ROInames)
    for iROI = 1:size(roiData,2)
        clear randMat;
        clear randMatAttn;
        clear randMatDec;
        for iRand = 1:10000
            for iSubj = 1:size(roiData,1)
                if ~isempty(roiData{iSubj,iROI})
                    randMat(iRand, iSubj) = roiData{iSubj,iROI}.bootAccuracy(randperm(1000,1));
                end
            end
            if ~strcmp(Analysis, 'Attn_4_Conds')
                % Get specific bootMat for Decision
                for iSubj = 1:size(Decision)
                    if ~isempty(roiData{Decision(iSubj),iROI})
                        randMatDec(iRand, iSubj) = roiData{Decision(iSubj),iROI}.bootAccuracy(randperm(1000,1));
                    end
                end
                
                % Get specific bootmat for Attention
                for iSubj = 1:size(Attention)
                    if ~isempty(roiData{Attention(iSubj),iROI})
                        randMatAttn(iRand, iSubj) = roiData{Attention(iSubj),iROI}.bootAccuracy(randperm(1000,1));
                    end
                end
            end
        end
        newPrctile = 100 - (5/NumberofAnalyses);
        bootMat(iROI,1) = prctile(mean(randMat,2),50);
        bootMat(iROI,2) = prctile(mean(randMat,2),newPrctile);
        
        if ~strcmp(Analysis, 'Attn_4_Conds')
            if ~isempty(Decision)
                bootMatDec(iROI,1) = prctile(mean(randMatDec,2),50);
                bootMatDec(iROI,2) = prctile(mean(randMatDec,2),newPrctile);
            end
            if ~isempty(Attention)
                bootMatAttn(iROI,1) = prctile(mean(randMatAttn,2),50);
                bootMatAttn(iROI,2) = prctile(mean(randMatAttn,2),newPrctile);
            end
        end
        
        for iSubj = 1:size(roiData,1)
            if ~isempty(roiData{iSubj,iROI})
                accMat(iSubj, iROI) = roiData{iSubj, iROI}.accuracy;
            else
                accMat(iSubj, iROI) = NaN;
            end
        end
    end
    
        accMat = accMat * 100;
    bootMat = bootMat * 100;
    
    if ~strcmp(Analysis, 'Attn_4_Conds')
        if ~isempty(Decision)
            bootMatAttn = bootMatAttn * 100;
            bootMatDec = bootMatDec * 100;
            accMatDec = accMat(Decision,:);
            accMatAttn = accMat(Attention,:);
        end
        
    else %compare decoding of illusion trajectory with decoding of local motion
        % Step 1, compute difference of the 2 for each subject
        for iROI = 1:size(roiData,2)/2
            accMatDiff(:,iROI) = accMat(:,iROI) - accMat(:,(size(roiData,2)/2)+iROI);
        end
        % Step 2, average it together by simply doing mean(accMatDiff,1)
        
        % Step 3, compute a null distribution of differences. For each
        % subject, take random samples by doing (trajectory) - (local
        % motion), average across subjects, repeat 10,000
        
        for iROI = 1:size(roiData,2)/2
            for iRand = 1:10000
                for iSubj = 1:size(roiData,1)
                    if ~isempty(roiData{iSubj,iROI})
                        randomIllusion(iRand, iSubj) = roiData{iSubj,iROI}.bootAccuracy(randperm(1000,1));
                        randomLocalMotion(iRand, iSubj) = roiData{iSubj, (size(roiData,2)/2)+iROI}.bootAccuracy(randperm(1000,1)); 
                    end
                end
            end
            
            temp = randomIllusion - randomLocalMotion;
            temp2(:,iROI) = mean(temp,2);
            bootMatDiff(iROI,1) = prctile(mean(temp2(:,iROI),2),50);
            bootMatDiff(iROI,2) = prctile(mean(temp2(:,iROI),2),newPrctile);
        end
        
    end
    
    
    
    keyboard
    if strcmp(Analysis, 'indiv')
        
        % Eli wants the order to be V1 / MT+ / V3A, so change things up
        ROInamesShort = {'V1', 'MT+', 'V3A'}';
        accMat = [accMat(:,2), accMat(:,1), accMat(:,3)];
        keyboard
        figure
        h = mybar(nanmean(accMat,1)', ...
            'yError', nanstd(accMat,1)'/length(accMat), ...
            'groupLabels', ROInamesShort, ...
            'xLabelText', 'visual area', ...
            'yLabelText', 'decoding accuracy (%)', ...
            'withinGroupColors', {Colors(1,:)}, ...
            'yAxisMin', 0, 'yAxisMax', 60);
        
        hline(mean(bootMat(:,1)));
        hline(mean(bootMat(:,2)));
        
        %leg = legend; leg.FontSize = 16;
        axis square
        h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
        hFig1 = figure(n); %handle of the most recent figure
        set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
        
        % then, break the data up into attn or decision data
        keyboard
        figure
        h = mybar([mean(accMatDec,1)', nanmean(accMatAttn,1)'], ...
            'yError', [std(accMatDec,1)'/length(Decision), ...
            nanstd(accMatAttn,1)'/length(Attention)], ...
            'groupLabels', ROInamesShort, ...
            'xLabelText', 'visual area', ...
            'yLabelText', 'decoding accuracy (%)', ...
            'withinGroupColors', {Colors(4,:), Colors(5,:)}, ...
            'yAxisMin', 0, 'yAxisMax', 60, ...
            'withinGroupLabels', {'decision task', 'attention task'});
        
        hline(mean(bootMat(:,1)));
        hline(mean(bootMat(:,2)));
        
        leg = legend; leg.FontSize = 16;
        axis square
        h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
        hFig1 = figure(n); %handle of the most recent figure
        set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
    elseif strcmp(Analysis, 'indiv+sanity')
        % Eli wants the order to be V1 / MT+ / V3A, so change things up
        ROInamesShort = {'V1', 'MT+', 'V3A', 'DUMMY'}';
        accMat = [accMat(:,2), accMat(:,1), accMat(:,3), accMat(:,4)];
        
        figure
        h = mybar(nanmean(accMat,1)', ...
            'yError', nanstd(accMat,1)'/length(accMat), ...
            'groupLabels', ROInamesShort, ...
            'yLabelText', 'decoding accuracy (%)', ...
            'withinGroupColors', {Colors(1,:)}, ...
            'yAxisMin', 0, 'yAxisMax', 60);
        
        hline(mean(bootMat(:,1)));
        hline(mean(bootMat(:,2)));
        
        axis square
        h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
        hFig1 = figure(n); %handle of the most recent figure
        set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
        savepdf(h, 'roiDecoding_stimlocs+V1.pdf')
        
        % then, break the data up into attn or decision data
        
        accMatDec = accMat(Decision,:);
        accMatAttn = accMat(Attention,:);
        figure
        h = mybar([mean(accMatDec,1)', nanmean(accMatAttn,1)'], ...
            'yError', [std(accMatDec,1)'/length(Decision), ...
            nanstd(accMatAttn,1)'/length(Attention)], ...
            'groupLabels', ROInamesShort, ...
            'xLabelText', 'visual area', ...
            'yLabelText', 'decoding accuracy (%)', ...
            'withinGroupColors', {Colors(4,:), Colors(5,:)}, ...
            'yAxisMin', 0, 'yAxisMax', 60, ...
            'withinGroupLabels', {'decision task', 'attention task'});
        
        hline(mean(bootMat(:,1)));
        hline(mean(bootMat(:,2)));
        
        leg = legend; leg.FontSize = 16;
        axis square
        h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
        hFig1 = figure(n); %handle of the most recent figure
        set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
        savepdf(hFig1, 'roiDecoding_stimlocs_AttnDec_Separate.pdf')
        
        keyboard
    else
        figure
        h = mybar(mean(accMat,1)', ...
            'yError', std(accMat)'/length(accMat), ...
            'groupLabels', ROInamesShort, ...
            'xLabelText', 'visual area', ...
            'yLabelText', 'decoding accuracy (%)', ...
            'withinGroupColors', {Colors(1,:)}, ...
            'yAxisMin', 0, 'yAxisMax', 60);
        
        hline(mean(bootMat(:,1)));
        hline(mean(bootMat(:,2)));
        
        %leg = legend; leg.FontSize = 16;
        axis square
        h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
        hFig1 = figure(n); %handle of the most recent figure
        set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
        if isequal(Conds, [3 4])
            savepdf(hFig1, 'roiDecoding_4Conds_LocalLeft_LocalRight.pdf')
        elseif isequal(Conds, [1 2])
            savepdf(hFig1, 'roiDecoding_4Conds_IllusionLeft_IllusionRIght.pdf')
        end
        
        keyboard
        if ~isempty(Decision) || ~isempty(Attention)
            % Then, let's do a graph that breaks that up by attention or
            % decision data
            
            accMatDec = accMat(Decision,:);
            accMatAttn = accMat(Attention,:);
            
            figure
            h = mybar([mean(accMatDec,1)', nanmean(accMatAttn,1)'], ...
                'yError', [std(accMatDec,1)'/length(Decision), ...
                std(accMatAttn,1)'/length(Attention)], ...
                'groupLabels', ROInamesShort, ...
                'xLabelText', 'visual area', ...
                'yLabelText', 'decoding accuracy (%)', ...
                'withinGroupColors', {Colors(4,:), Colors(5,:)}, ...
                'yAxisMin', 0, 'yAxisMax', 60, ...
                'withinGroupLabels', {'decision task', 'attention task'});
            
            hline(mean(bootMat(:,1)));
            hline(mean(bootMat(:,2)));
            
            %leg = legend; leg.FontSize = 16;
            axis square
            h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
            hFig1 = figure(n); %handle of the most recent figure
            set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
            keyboard
        end
        
        
    end
elseif IllusionVsNoIllusion == 1
    if ~strcmp(Analysis, 'Attn_4_Conds')
        for iROI = 1:length(ROInames)
            for iSubj = 1:size(roiData,1)
                if ~isempty(roiData{iSubj,iROI})
                    illusionMat(iSubj, iROI, 1) = mean([roiData{iSubj,iROI}.dd2.ehdr(1), roiData{iSubj,iROI}.dd2.ehdr(3)]);
                    illusionMat(iSubj, iROI, 2) = roiData{iSubj, iROI}.dd2.ehdr(2);
                    
                    illusionLvR(iSubj, iROI, 1) = roiData{iSubj, iROI}.dd2.ehdr(1);
                    illusionLvR(iSubj, iROI, 2) = roiData{iSubj, iROI}.dd2.ehdr(3);
                else
                    illusionMat(iSubj, iROI, 1) = NaN;
                    illusionMat(iSubj, iROI, 2) = NaN;
                    
                    illusionLvR(iSubj, iROI, 1) = NaN;
                    illusionLvR(iSubj, iROI, 2) = NaN;
                end
            end
        end
    else
        for iROI = 1:size(roiData,2)
            for iSubj = 1:size(roiData,1)
                illusionMat(iSubj, iROI ,1) = mean([roiData{iSubj,iROI}.dd2.ehdr(1), roiData{iSubj,iROI}.dd2.ehdr(3)]);
                illusionMat(iSubj, iROI, 2) = roiData{iSubj, iROI}.dd2.ehdr(2);
            end
        end
    end
    keyboard

    
    %running Anovas
    if strcmp(Analysis, 'Attn_4_Conds')
        keyboard
        disp('------ Anova Illusion vs. Local Motion Only ------')
        
        illVsLocal = zeros(size(illusionMat));
        illVsLocal(:,1:3,:) = 1;
        Conditions = zeros(size(illusionMat));
        Conditions(:,:,2) = 1;
        

    else
        disp('------ T-Tests Illusion vs. No Illusion ------')
        for iROI=1:3
            [h,p, ci, stats] = ttest(illusionMat(Decision,iROI,1),illusionMat(Decision,iROI,2));
            disp(sprintf('ttest for ROI %s & Attend Stimulus (Decision): t(%d)=%.2f, p=%.4f', ROInamesShort{iROI}, length(Decision)-1, stats.tstat, p))
            disp('   ');
            [h,p, ci, stats] = ttest(illusionMat(Attention,iROI,1),illusionMat(Attention,iROI,2));
            disp(sprintf('ttest for ROI %s & Attend Target (Attention): t(%d)=%.2f, p=%.4f', ROInamesShort{iROI}, length(Attention)-1, stats.tstat, p))
            disp('   ');
            [h,p, ci, stats] = ttest(illusionMat(:,iROI,1),illusionMat(:,iROI,2));
            disp(sprintf('ttest for ROI %s & All Data: t(%d)=%.2f, p=%.4f', ROInamesShort{iROI}, size(illusionMat,1)-1, stats.tstat, p))
            disp('   ');
            disp('----------')
        end
    end
    
    
    Task=zeros(size(illusionMat));
    Task(Attention,:,:)=1;
    Conditions = zeros(size(illusionMat));
    Conditions(:,:,2) = 1;
    ROIs = zeros(size(illusionMat));
    ROIs(:,1,:) = 1; ROIs(:,2,:) = 2; ROIs(:,3,:) = 3;
    ROI = [illusionMat(:,1,1); illusionMat(:,1,2)];
   
    P = anovan(ROI, {[Task(:,1,1); Task(:,1,2)], [Conditions(:,1,1); Conditions(:,1,2)]},'model','interaction','varnames',{'Task','Conditions'});
    
    % T test between L and R
    disp('------ T-Tests Left vs. Right ------')
    for iROI = 1:3
        [h,p, ci, stats] = ttest(illusionLvR(:,iROI,1),illusionLvR(:,iROI,2));
        disp(sprintf('ttest for ROI %s & All Data: t(%d)=%.2f, p=%.4f', ROInamesShort{iROI}, size(illusionMat,1)-1, stats.tstat, p))
        disp('   ')
    end
        
    % Anova between condition and ROI to test for interaction
    P = anovan(illusionMat(:), {ROIs(:), Conditions(:)}, 'model','interaction','varnames',{'ROIs','Conditions'});
    
    
    keyboard
    
    keyboard
    % Eli wants the order to be V1 / MT+ / V3A, so change things up
    ROInamesShort = {'V1', 'MT+', 'V3A'}';
    illusionMat(:,:,1) = [illusionMat(:,2,1), illusionMat(:,1,1), illusionMat(:,3,1)];
    illusionMat(:,:,2) = [illusionMat(:,2,2), illusionMat(:,1,2), illusionMat(:,3,2)];
    keyboard
    figure
    h = mybar([nanmean(illusionMat(:,:,1),1)', nanmean(illusionMat(:,:,2),1)'], ...
        'yError', [nanstd(illusionMat(:,:,1),1)'/size(illusionMat,1), ...
        nanstd(illusionMat(:,:,2),1)'/size(illusionMat,1)], ...
        'groupLabels', ROInamesShort, ...
        'yLabelText', sprintf('fMRI response \n(%%change image intensity)'), ...
        'withinGroupColors', {Colors(1,:), Colors(2,:)}, ...
        'yAxisMin', 0, 'yAxisMax', 4, ...
        'withinGroupLabels', {'illusion', 'no illusion'});
    leg = legend; leg.FontSize = 16;
    axis square
    h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
    hFig1 = figure(n); %handle of the most recent figure
    set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
    
    savepdf(h, 'boldResponse_IllvsNoIll.pdf')
    
    keyboard
    elseif plotModulationIndex == 1
        
        if ~strcmp(Analysis, 'Attn_4_Conds')
            for iROI = 1:length(ROInames)
                for iSubj = 1:size(roiData,1)
                    if ~isempty(roiData{iSubj,iROI})
                        illusionMat(iSubj, iROI, 1) = mean([roiData{iSubj,iROI}.dd2.ehdr(1), roiData{iSubj,iROI}.dd2.ehdr(3)]);
                        illusionMat(iSubj, iROI, 2) = roiData{iSubj, iROI}.dd2.ehdr(2);
                        
                        illusionLvR(iSubj, iROI, 1) = roiData{iSubj, iROI}.dd2.ehdr(1);
                        illusionLvR(iSubj, iROI, 2) = roiData{iSubj, iROI}.dd2.ehdr(3);
                    else
                        illusionMat(iSubj, iROI, 1) = NaN;
                        illusionMat(iSubj, iROI, 2) = NaN;
                        
                        illusionLvR(iSubj, iROI, 1) = NaN;
                        illusionLvR(iSubj, iROI, 2) = NaN;
                    end
                end
            end
        end
        
            % Eli wants the order to be V1 / MT+ / V3A, so change things up
    ROInamesShort = {'V1', 'MT+', 'V3A'}';
    illusionMat(:,:,1) = [illusionMat(:,2,1), illusionMat(:,1,1), illusionMat(:,3,1)];
    illusionMat(:,:,2) = [illusionMat(:,2,2), illusionMat(:,1,2), illusionMat(:,3,2)];
    
        % Also, run 2X2 anova for EACH ROI (task and condition as factors)
        for iROI = 1:length(ROInames)
            ROIdata = [illusionMat(:,iROI,1),  illusionMat(:,iROI,2)];
            factor_task = zeros(size(ROIdata));
            factor_task(Decision,:) = 1; factor_task(Attention,:) = 2;
            factor_cond = zeros(size(ROIdata));
            factor_cond(:,1) = 1; factor_cond(:,2) = 2;
            P = anovan(ROIdata(:), {factor_task(:), factor_cond(:)}, 'model','interaction','varnames',{'Task','Conditions'});
            
            %also do post-hoc t-tests to show an effect of condition
            [h,p, ci, stats] = ttest(illusionMat(Decision,iROI,1),illusionMat(Decision,iROI,2));
            disp(sprintf('ttest for ROI %s & Attend Stimulus (Decision): t(%d)=%.2f, p=%.4f', ROInamesShort{iROI}, length(Decision)-1, stats.tstat, p))
            [h,p, ci, stats] = ttest(illusionMat(Attention,iROI,1),illusionMat(Attention,iROI,2));
            disp(sprintf('ttest for ROI %s & Attend Target (Attention): t(%d)=%.2f, p=%.4f', ROInamesShort{iROI}, length(Attention)-1, stats.tstat, p))
            [h,p, ci, stats] = ttest(illusionMat(:,iROI,1),illusionMat(:,iROI,2));
            disp(sprintf('ttest for ROI %s & ALL DATA: t(%d)=%.2f, p=%.4f', ROInamesShort{iROI}, length(Attention)-1, stats.tstat, p))
            disp('  ');
            disp('----------');
        end
        
        % (Figure 2D) 1) replot data with 4 bars per ROI (illusion, no
        % illusion, for decision and attention)
        
        [mean(illusionMat(Decision,1,1)), mean(illusionMat(Decision,1,2)), ...
         mean(illusionMat(Attention,1,1)), mean(illusionMat(Attention,1,2)); ...
         mean(illusionMat(Decision,2,1)), mean(illusionMat(Decision,2,2)), ...
         mean(illusionMat(Attention,2,1)), mean(illusionMat(Attention,2,2)); ...
         mean(illusionMat(Decision,3,1)), mean(illusionMat(Decision,3,2)), ...
         nanmean(illusionMat(Attention,3,1)), nanmean(illusionMat(Attention,3,2))];
     
     stdForGraph = [std(illusionMat(Decision,1,1)), std(illusionMat(Decision,1,2)), ...
         std(illusionMat(Attention,1,1)), std(illusionMat(Attention,1,2)); ...
         std(illusionMat(Decision,2,1)), std(illusionMat(Decision,2,2)), ...
         std(illusionMat(Attention,2,1)), std(illusionMat(Attention,2,2)); ...
         std(illusionMat(Decision,3,1)), std(illusionMat(Decision,3,2)), ...
         nanstd(illusionMat(Attention,3,1)), nanstd(illusionMat(Attention,3,2))];
     stdForGraph(:,1) = stdForGraph(:,1) / sqrt(length(Decision));
     stdForGraph(:,2) = stdForGraph(:,2) / sqrt(length(Attention));
     stdForGraph(:,3) = stdForGraph(:,3) / sqrt(length(Decision));
     stdForGraph(:,4) = stdForGraph(:,4) / sqrt(length(Attention));
           keyboard
            
        figure
        h = mybar([mean(illusionMat(Decision,1,1)), mean(illusionMat(Decision,1,2)), ...
         mean(illusionMat(Attention,1,1)), mean(illusionMat(Attention,1,2)); ...
         mean(illusionMat(Decision,2,1)), mean(illusionMat(Decision,2,2)), ...
         mean(illusionMat(Attention,2,1)), mean(illusionMat(Attention,2,2)); ...
         mean(illusionMat(Decision,3,1)), mean(illusionMat(Decision,3,2)), ...
         nanmean(illusionMat(Attention,3,1)), nanmean(illusionMat(Attention,3,2))], ...
            'yError', stdForGraph, ...
            'groupLabels', {ROInamesShort{1}, ROInamesShort{2}, ROInamesShort{3}}, ...
            'xLabelText', 'ROI', ...
            'yLabelText', 'FILL IN', ...
            'withinGroupColors', {Colors(4,:), Colors(4,:), Colors(5,:), Colors(5,:)}, ...
            'yAxisMin', 0, 'yAxisMax', 4, ...
            'withinGroupLabels', {'Attend stimulus', 'Attend stimulus', 'Attend target', 'Attend target'});
        %'withinGroupColors', {ExpData.Inds.Spat.Color.Scram, ExpData.Inds.Ret.Color.Image}, ...
        leg = legend; leg.FontSize = 14;
        axis square
        h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
        hFig1 = figure(n); %handle of the most recent figure
        set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
                savepdf(hFig1, 'Fig2D.pdf')

        
        keyboard
        if ~strcmp(Analysis, 'Attn_4_Conds')
            for iROI = 1:length(ROInames)
                for iSubj = 1:size(roiData,1)
                    if ~isempty(roiData{iSubj,iROI})
                        illVsNoIll(iSubj,iROI) = roiData{iSubj,iROI}.dd2.ehdr
                    end
                end
            end
        end
            
            
        keyboard
        
        
    modIndex.mtDec = (illvsNoIll.mtMat(Decision,1)-illvsNoIll.mtMat(Decision,2)) ./ (illvsNoIll.mtMat(Decision,1)+illvsNoIll.mtMat(Decision,2));
    modIndex.mtAttn = (illvsNoIll.mtMat(Attention,1)-illvsNoIll.mtMat(Attention,2)) ./ (illvsNoIll.mtMat(Attention,1)+illvsNoIll.mtMat(Attention,2));
    modIndex.v1Dec = (illvsNoIll.v1Mat(Decision,1)-illvsNoIll.v1Mat(Decision,2)) ./ (illvsNoIll.v1Mat(Decision,1)+illvsNoIll.v1Mat(Decision,2));
    modIndex.v1Attn = (illvsNoIll.v1Mat(Attention,1)-illvsNoIll.v1Mat(Attention,2)) ./ (illvsNoIll.v1Mat(Attention,1)+illvsNoIll.v1Mat(Attention,2));
    modIndex.ipsDec = (illvsNoIll.ipsMat(Decision,1)-illvsNoIll.ipsMat(Decision,2)) ./ (illvsNoIll.ipsMat(Decision,1)+illvsNoIll.ipsMat(Decision,2));
    modIndex.ipsAttn = (illvsNoIll.ipsMat(Attention,1)-illvsNoIll.ipsMat(Attention,2)) ./ (illvsNoIll.ipsMat(Attention,1)+illvsNoIll.ipsMat(Attention,2));

    figure
    h = mybar([mean(modIndex.v1Dec), mean(modIndex.v1Attn); ...
        mean(modIndex.mtDec), mean(modIndex.mtAttn); ...
        mean(modIndex.ipsDec), nanmean(modIndex.ipsAttn)], ...
        'yError', [std(modIndex.v1Dec)/length(modIndex.v1Dec), std(modIndex.v1Attn)/length(modIndex.v1Attn); ...
        std(modIndex.mtDec)/length(modIndex.mtDec), std(modIndex.mtAttn)/length(modIndex.mtAttn); ...
        std(modIndex.ipsDec)/length(modIndex.ipsDec), nanstd(modIndex.ipsAttn)/length(modIndex.ipsAttn)], ...
        'groupLabels', {'v1', 'mt', 'ips'}, ...
        'yLabelText', 'Modulation index', ...
        'withinGroupColors', {Colors(4,:), Colors(5,:)}, ...
        'yAxisMin', -.1, 'yAxisMax', .3);
else
    keyboard
end
keyboard




%
%
%
%
% keyboard
% ROInamesShort = {};
% for i = 1:length(ROInames)
%     temp = strsplit(ROInames{i}, '_');
%     ROInamesShort{i,1} = temp{end};
% end
% if plotIllvsNoIll == 1
%     if strcmp(Analysis, 'indiv')
%         keyboard
%         %plot the data
%         figure
%         h = mybar([mean(illvsNoIll.v1Mat(:,1)), mean(illvsNoIll.v1Mat(:,2)); mean(illvsNoIll.mtMat(:,1)), mean(illvsNoIll.mtMat(:,2)); nanmean(illvsNoIll.ipsMat(:,1)), nanmean(illvsNoIll.ipsMat(:,2))], ...
%             'groupLabels', {'v1', 'mt', 'ips'}, ...
%             'xLabelText', 'visual area', ...
%             'yLabelText', 'fMRI response (% change img intensity)', ...
%             'withinGroupColors', {Colors(1,:), Colors(2,:)}, ...
%             'yAxisMin', 0, 'yAxisMax', 4, ...
%             'yTick', [0 2 4], ...
%             'withinGroupLabels', {'illusion', 'no illusion'});
%         leg = legend; leg.FontSize = 16;
%         axis square
%         h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
%         hFig1 = figure(n); %handle of the most recent figure
%         set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
%
%     elseif strcmp(Analysis, 'group')
%         figure
%         h = mybar([mean(illvsNoIll.mtMat(:,1)), mean(illvsNoIll.mtMat(:,2)); mean(illvsNoIll.ipsMat(:,1)), mean(illvsNoIll.ipsMat(:,2)); ...
%             nanmean(illvsNoIll.stsMat(:,1)), nanmean(illvsNoIll.stsMat(:,2)) ; ...
%             nanmean(illvsNoIll.poMat(:,1)), nanmean(illvsNoIll.poMat(:,2)) ; ...
%             nanmean(illvsNoIll.splMat(:,1)), nanmean(illvsNoIll.splMat(:,2))], ...
%             'groupLabels', ROInamesShort', ...
%             'xLabelText', 'visual area', ...
%             'yLabelText', 'fMRI response (% change img intensity)', ...
%             'withinGroupColors', {Colors(1,:), Colors(2,:)}, ...
%             'yAxisMin', 0, 'yAxisMax', 4, ...
%             'yTick', [0 2 4], ...
%             'withinGroupLabels', {'illusion', 'no illusion'});
%         leg = legend; leg.FontSize = 16;
%         axis square
%         h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
%         hFig1 = figure(n); %handle of the most recent figure
%         set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
%
%     end
% elseif plotModulationIndex == 1
%     modIndex.mtDec = (illvsNoIll.mtMat(Decision,1)-illvsNoIll.mtMat(Decision,2)) ./ (illvsNoIll.mtMat(Decision,1)+illvsNoIll.mtMat(Decision,2));
%     modIndex.mtAttn = (illvsNoIll.mtMat(Attention,1)-illvsNoIll.mtMat(Attention,2)) ./ (illvsNoIll.mtMat(Attention,1)+illvsNoIll.mtMat(Attention,2));
%     modIndex.v1Dec = (illvsNoIll.v1Mat(Decision,1)-illvsNoIll.v1Mat(Decision,2)) ./ (illvsNoIll.v1Mat(Decision,1)+illvsNoIll.v1Mat(Decision,2));
%     modIndex.v1Attn = (illvsNoIll.v1Mat(Attention,1)-illvsNoIll.v1Mat(Attention,2)) ./ (illvsNoIll.v1Mat(Attention,1)+illvsNoIll.v1Mat(Attention,2));
%     modIndex.ipsDec = (illvsNoIll.ipsMat(Decision,1)-illvsNoIll.ipsMat(Decision,2)) ./ (illvsNoIll.ipsMat(Decision,1)+illvsNoIll.ipsMat(Decision,2));
%     modIndex.ipsAttn = (illvsNoIll.ipsMat(Attention,1)-illvsNoIll.ipsMat(Attention,2)) ./ (illvsNoIll.ipsMat(Attention,1)+illvsNoIll.ipsMat(Attention,2));
%
%     figure
%     h = mybar([mean(modIndex.v1Dec), mean(modIndex.v1Attn); ...
%         mean(modIndex.mtDec), mean(modIndex.mtAttn); ...
%         mean(modIndex.ipsDec), nanmean(modIndex.ipsAttn)], ...
%         'yError', [std(modIndex.v1Dec)/length(modIndex.v1Dec), std(modIndex.v1Attn)/length(modIndex.v1Attn); ...
%         std(modIndex.mtDec)/length(modIndex.mtDec), std(modIndex.mtAttn)/length(modIndex.mtAttn); ...
%         std(modIndex.ipsDec)/length(modIndex.ipsDec), nanstd(modIndex.ipsAttn)/length(modIndex.ipsAttn)], ...
%         'groupLabels', {'v1', 'mt', 'ips'}, ...
%         'yLabelText', 'Modulation index', ...
%         'withinGroupColors', {Colors(4,:), Colors(5,:)}, ...
%         'yAxisMin', -.1, 'yAxisMax', .3);
%
%     keyboard
% elseif plotLvsR == 1
%     keyboard
%     if strcmp(Analysis, 'indiv')
%         figure
%         h = mybar([mean(DecodingLvsR.mt)*100; mean(DecodingLvsR.v1)*100; nanmean(DecodingLvsR.ips)*100], ...
%             'yError', [(std(DecodingLvsR.mt)/sqrt(length(DecodingLvsR.mt)))*100; ...
%             (std(DecodingLvsR.v1)/sqrt(length(DecodingLvsR.mt)))*100; ...
%             (nanstd(DecodingLvsR.ips)/sqrt(length(DecodingLvsR.mt)))*100], ...
%             'groupLabels', {'v1', 'mt', 'ips'}, ...
%             'xLabelText', 'visual areas', ...
%             'yLabelText', 'decoding accuracy (%)', ...
%             'withinGroupColors', {Colors(3,:)}, ...
%             'yAxisMin', 0, 'yAxisMax', 100);
%
%         %hline(50); %chance decoding
%         %hline(prctile(DecodingLvsR.mtBoot(:), 95)*100); %95th percentile decoding
%
%         %leg = legend; leg.FontSize = 16;
%         axis square
%         h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
%         hFig1 = figure(n); %handle of the most recent figure
%         set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
%         keyboard
%     elseif strcmp(Analysis, 'group')
%         keyboard
%         DecodingLvsR.mt(DecodingLvsR.mt==0)=NaN;
%         DecodingLvsR.sts(DecodingLvsR.sts==0)=NaN;
%         DecodingLvsR.ips(DecodingLvsR.ips==0)=NaN;
%         DecodingLvsR.po(DecodingLvsR.po==0)=NaN;
%         DecodingLvsR.spl(DecodingLvsR.spl==0)=NaN;
%
%         violinplot([DecodingLvsR.mt , DecodingLvsR.sts , ...
%             DecodingLvsR.ips , DecodingLvsR.po , DecodingLvsR.spl])
%
%
%         keyboard
%
%         figure
%         h = mybar([nanmean(DecodingLvsR.mt)*100; nanmean(DecodingLvsR.sts)*100; ...
%             nanmean(DecodingLvsR.ips)*100; nanmean(DecodingLvsR.po)*100; ...
%             nanmean(DecodingLvsR.spl)*100], ...
%             'groupLabels', ROInamesShort', ...
%             'xLabelText', 'visual area', ...
%             'yLabelText', 'fMRI response (% change img intensity)', ...
%             'withinGroupColors', {Colors(1,:), Colors(2,:)}, ...
%             'yAxisMin', 0, 'yAxisMax', 4, ...
%             'yTick', [0 2 4], ...
%             'withinGroupLabels', {'illusion', 'no illusion'});
%         leg = legend; leg.FontSize = 16;
%         axis square
%         h =  findobj('type','figure'); n = length(h); if n == 0; n=1; end; %find how many figures there are
%         hFig1 = figure(n); %handle of the most recent figure
%         set(hFig1, 'Position', [5+(10*(n-1)) 23 30 23])
%     end
% elseif decodeAllTogether == 1
%     %new way to get nTrials - is this right?
%     nTrials = size(DecodingLvsR.mt.allEhdr,2) / 3;
%
%     for iROI = 1:length(ROInames)
%         clear betas; clear accuracy; clear bootAccuracy;
%         switch iROI
%             case 1
%                 goodVox.mt = DecodingLvsR.mt.r2>prctile(DecodingLvsR.mt.r2,50);
%                 betas = DecodingLvsR.mt.allEhdr(goodVox.mt,[1:nTrials (nTrials*2)+1:nTrials*3])';
%             case 2
%                 goodVox.v1 = DecodingLvsR.v1.r2>prctile(DecodingLvsR.v1.r2,50);
%                 betas = DecodingLvsR.v1.allEhdr(goodVox.v1,[1:nTrials (nTrials*2)+1:nTrials*3])';
%
%             case 3
%                 goodVox.ips = DecodingLvsR.ips.r2>prctile(DecodingLvsR.ips.r2,50);
%                 betas = DecodingLvsR.ips.allEhdr(goodVox.ips,[1:nTrials (nTrials*2)+1:nTrials*3])';
%         end
%         keyboard
%         % To compare conditions 1 vs 3 (left vs right)
%         % betas = DecodingLvsR.mt.allEhdr(goodVox.mt,[1:nTrials (nTrials*2)+1:nTrials*3])';
%
%         betas = zscore(betas);
%         % create a grouping variable
%         group = cat(1, repmat('l', nTrials, 1), repmat('r', nTrials, 1));
%         disp('  '); disp(sprintf('ROI = %s', ROInames{iROI})); disp('  ');
%         disp(sprintf('(plotROI): Computing accuracy for L vs R'));
%         % create partitioning scheme
%         CVO = cvpartition(group, 'KFold', nTrials*2);
%         % %CVO = cvpartition(group, 'KFold', 16); %ONLY IF YOU WANT TO DO LEAVE-1-RUN-OUT (NOTE: only 4 trials, not all 8)
%         %CVO = cvpartition(group, 'LeaveOut', 1);
%         %CVO = cvpartition(group, 'HoldOut', 0.5);
%
%         err = zeros(CVO.NumTestSets,1);
%         for i = 1:CVO.NumTestSets
%             trIdx = CVO.training(i);
%             teIdx = CVO.test(i);
%             ytest = classify(betas(teIdx,:),betas(trIdx,:),group(trIdx), 'diagLinear');
%             err(i) = sum(~strcmp(ytest,group(teIdx)));
%         end
%         accuracy = 1 - sum(err)/sum(CVO.TestSize);
%         disp(sprintf('Classification accuracy: %f', accuracy));
%
%         disp(sprintf('(plotROI): Computing bootstrap for L vs R'));
%         nBoot = 1000;
%         bootAccuracy = zeros(nBoot, 1);
%         for iBoot=1:nBoot;
%             gNull = group(randperm(size(group,1)));
%             for i = 1:CVO.NumTestSets
%                 trIdx = CVO.training(i);
%                 teIdx = CVO.test(i);
%                 ytest = classify(betas(teIdx,:),betas(trIdx,:),gNull(trIdx), 'diagLinear');
%                 err(i) = sum(~strcmp(ytest,gNull(teIdx)));
%             end
%             bootAccuracy(iBoot) = 1 - sum(err)/sum(CVO.TestSize);
%         end
%         disp(sprintf('Bootstrapped threshold is: %f', prctile(bootAccuracy, 95)));
%     end
%     keyboard
%     clear randMat
%     for i=1:10000
%         for j=1:19
%             randMat(i,j) = DecodingLvsR.mtBoot(j,randperm(1000,1));
%         end
%     end
%     foo = mean(randMat,2);
%     prctile(foo, 95)
%
%     keyboard
%     % THEN DO THIS FOR V1 AND IPS
% elseif plotScatter == 1
%     keyboard
%
%
% end
%
% % keyboard
% % illvsNoIll.Idx.mtDec = (illvsNoIll.mtMat(Decision,1)-illvsNoIll.mtMat(Decision,2)) ./  (illvsNoIll.mtMat(Decision,1)+illvsNoIll.mtMat(Decision,2));
% % illvsNoIll.Idx.mtAttn = (illvsNoIll.mtMat(Attention,1)-illvsNoIll.mtMat(Attention,2)) ./  (illvsNoIll.mtMat(Attention,1)+illvsNoIll.mtMat(Attention,2));
% % illvsNoIll.Idx.v1Dec = (illvsNoIll.v1Mat(Decision,1)-illvsNoIll.v1Mat(Decision,2)) ./  (illvsNoIll.v1Mat(Decision,1)+illvsNoIll.v1Mat(Decision,2));
% % illvsNoIll.Idx.v1Attn = (illvsNoIll.v1Mat(Attention,1)-illvsNoIll.v1Mat(Attention,2)) ./  (illvsNoIll.v1Mat(Attention,1)+illvsNoIll.v1Mat(Attention,2));
% % illvsNoIll.Idx.ipsDec = (illvsNoIll.ipsMat(Decision,1)-illvsNoIll.ipsMat(Decision,2)) ./  (illvsNoIll.ipsMat(Decision,1)+illvsNoIll.ipsMat(Decision,2));
% % illvsNoIll.Idx.ipsAttn = (illvsNoIll.ipsMat(Attention,1)-illvsNoIll.ipsMat(Attention,2)) ./  (illvsNoIll.ipsMat(Attention,1)+illvsNoIll.ipsMat(Attention,2));
%
% clear randMat
% for i=1:10000
%     for j=1:19
%         randMat(i,j) = DecodingLvsR.mtBoot(j,randperm(1000,1));
%     end
% end
% pct95 = mean(randMat,2);
%         disp(prctile(pct95, 95))
%
% clear randMat
% for i=1:10000
%     for j=1:19
%         randMat(i,j) = DecodingLvsR.v1Boot(j,randperm(1000,1));
%     end
% end
% pct95 = mean(randMat,2);
%         disp(prctile(pct95, 95))
%
% clear randMat
% for i=1:10000
%     for j=1:18
%         randMat(i,j) = DecodingLvsR.ipsBoot(j,randperm(1000,1));
%     end
% end
% pct95 = mean(randMat,2);
%         disp(prctile(pct95, 95))
%