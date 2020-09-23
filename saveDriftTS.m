mrQuit
close all
clear all
tic
if exist('v','var')
    deleteView(v);
end
% directionOrIllusion = 2;%1 for direction, 2 for illusion
expNum=2;
% crossValidation = 'trial';%'trial','run'
% locThresh = 0.2; %only include voxels with coherence>locThresh, for each of the two localizers.
% thrsh = num2str(locThresh,'%.2f');
% locStr = thrsh([1 3:4]);

% ROIs = {'lV1','lV2','lV3','lhV4','lLO1','lLO2','lTO1','lTO2','lV3a','lV3b',...
%     'rV1','rV2','rV3','rhV4','rLO1','rLO2','rTO1','rTO2','rV3a','rV3b'};

ROIs = {'lV1','lV2','lV3','lhV4','lLO1','lLO2','lTO1','lTO2','lV3a','lV3b','lV01','lV02',...
    'rV1','rV2','rV3','rhV4','rLO1','rLO2','rTO1','rTO2','rV3a','rV3b','rV01','rV02'};
% ROIs = {'lV02','rV02'};
% ROIs = {'lV1','lV2','lV3','lhV4','lLO1','lLO2','lTO1','lTO2','lV3a','lV3b','lV01',...
%     'rV1','rV2','rV3','rhV4','rLO1','rLO2','rTO1','rTO2','rV3a','rV3b','rV01'};
% ROIs = {'lV1','rV1'};
switch expNum
    case 1
        expName = '3conds';
        dataFolder = '/Volumes/Drift/Drift_Decision_3_Conds/';
        subNames = {'s002920180301','s004020180427','s004720180830','s005020180430',...
            's005920181115','s006320181015','s003320180319','s004720180628',...
            's005120180524','s005720180906','s006020181001','s007320181116'};
        %         subNames = {'s005120180524','s005720180906','s006020181001','s007320181116'};
        %didn't work: s004920180329, s005520180628
        % bad coverage: s004720180628, s005720180906
        %why doesn't driftSubjDir.m have s0060?
        %where is s006220190107?
    case 2
        expName = 'attn';
        dataFolder = '/Volumes/Drift/DriftIllusion_Attention/';
        subNames = {'s007420190128','s007420190207','s008820190422','s005920190129',...
            's009320190422','s006220190107','s007720190204','s006320190211','s007920190114',...
            's008020190201'};
        %         subNames = {'s006320190211','s007920190114','s007420190128',...
        %             's008020190201','s007420190207','s008820190422','s005920190129',...
        %             's009320190422','s006220190107','s007720190204'};
        %doesn't have averages: s007720190214, s007320190227, s007420190213
        %didn't work: s007420190128 (no coverage of lV02)
        %missing second averages (or just stimfile)? s005720190222 - has no Occm control!
        %missing both controls: s007720190214,
        %no controls, no motion comp, no averages: s007320190227, s007420190213
    case 3
        %doesn't have the localizers??
        expName = '4conds';
        dataFolder = '/Volumes/Drift/DriftIllusion_Attn_4_Conds/';
        subNames = {'s005520190401','s008620190415','s008720190405','s008820190318','s009320190415'};
    case 4
        expName = '3tdata';
        dataFolder = '/Volumes/Drift/DriftIllusion_Attention/3TData/';
        subNames = {'s007320190227_3T','s007620190320_3T','s008620190417_3T',...
            's008720190313_3T','s009220190327_3T'};%no localizers/averages
end

numRois = length(ROIs);
% subNames = {'s005920190129','s006220190107','s007720190204'};
% subNames = {'s005920190129','s006220190107'};
numSubs = length(subNames);
% dataFolder = '/Volumes/ZviDrive/noah/drift_analysis_attn/';
saveFolder = '~/noah/';
cd(dataFolder);

nhdr = 3;
hdrlen = 1;
matchScan=1;
for isub=1:numSubs
    cd(subNames{isub});
    
    v=newView;
    
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
    v = viewSet(v,'curGroup','Concatenation');
    erFilename =  dir('Concatenation/erAnal/erAnal.mat');
    if isempty(erFilename)
        erFilename =  dir('Concatenation/erAnal/erAnal_24*.mat');
    end
    clear d
    d = viewGet(v, 'd');
    if isempty(d)
        v=loadAnalysis(v,erFilename.name,fullfile(pwd,'/Concatenation/erAnal/'));
        d = viewGet(v, 'd');
    end
    dd=d;
    %first use deconvolution to model HRF as mean response
    [stimvol stimNames var] = getStimvol(v,{'groupNum=[1 3]','groupNum=2'},'taskNum=1','phaseNum=1');
    d.stimvol = stimvol;
    d = makescm(d, 16, 0, stimvol);%design matrix convolved with a 16 TR canonical HRF, with 2 conditions - 1+3 and 2
    scm = makescm(v, 1, 0, dd.stimvol);%design matrix with delta functions, with 3 conditions
    % Find how many scans there are by the name of the Concatenation group
    concatInfo = viewGet(v, 'concatInfo');
    subNumScans(isub) = concatInfo.n;
    
    % construct an scm with individual trials
    trialscm = [];
    nTrials = length(dd.stimvol{1});
    subNumTrials(isub) = nTrials;
    for iCond=1:3
        for iTrial=1:nTrials
            temp = dd;
            temp.stimvol{iCond} = temp.stimvol{iCond}(iTrial);
            foo = makescm(v, 1, 0, temp.stimvol);
            trialscm = cat(2, trialscm, foo(:,iCond));
        end
    end
    
    %     keyboard
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
            
            
            %         % convolve scm with hemodynamic response
            %         scmHdr = conv2(scm, hdr');
            %         scmHdr = scmHdr(1:size(dd.scm,1),:);
            
            
            % trial-wise-scm, convolved with cannonical hdr
            trialscmhdr = conv2(trialscm,hdr');
            trialscmhdr = trialscmhdr(1:size(rois.tSeries,2),:);
            %get betas
            driftBetas{isub,iroi} = getr2timecourse(rois.tSeries, size(trialscmhdr,2), hdrlen, trialscmhdr);
            
            condscm = conv2(scm,hdr');
            condscm = condscm(1:size(rois.tSeries,2),:);
            meanRoiBetas{isub,iroi} = getr2timecourse(mean(rois.tSeries),nhdr,hdrlen,condscm);
        else
            driftBetas{isub,iroi} = [];
            meanRoiBetas{isub,iroi}=[];
            
        end
        
        %get localizer timeseries
        %     matchScan = viewGet(v,'curScan');
        stimlocCorrData = loadROIcoranalMatching(v, roiname, stimlocScan, avgGroupNum, matchScan, concatGroupNum);
        eyelocCorrData = loadROIcoranalMatching(v, roiname, eyelocScan, avgGroupNum, matchScan, concatGroupNum);
        stimLoc{isub,iroi} = stimlocCorrData{1};
        eyeLoc{isub,iroi} = eyelocCorrData{1};
        
        subRoiNans{isub,iroi} = roiNans;
        %         temp = trialReg.r2>prctile(trialReg.r2,50) & stimlocCorrData{1}.co(~roiNans)>locThresh & eyelocCorrData{1}.co(~roiNans)>locThresh;
        %         res{1}.sharedVox = sum(temp);
        %         res{1}.totalVox = length(temp);
        
    end
    deleteView(v);
    cd ..
end
%%
varinfo=whos('driftBetas');
saveopt='';
if varinfo.bytes >= 2^31
    saveopt='-v7.3';
end
save([saveFolder 'driftTS_' expName '.mat'],'ROIs','subNames','dataFolder',...
    'stimLoc','eyeLoc','driftBetas','meanRoiBetas','subRoiNans','subNumTrials','subNumScans',saveopt);

toc


