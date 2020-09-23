mrQuit
close all
clear all
tic
if exist('v','var')
    deleteView(v);
end
directionOrIllusion = 2;%1 for direction, 2 for illusion
expNum=1;
crossValidation = 'trial';%'trial','run'
locThresh = 0.2; %only include voxels with coherence>locThresh, for each of the two localizers.
thrsh = num2str(locThresh,'%.2f');
locStr = thrsh([1 3:4]);
ROIs = {'lV1','rV1'};
ROIs = {'lV1','lV2','lV3','lhV4','lLO1','lLO2','lTO1','lTO2','lV3a','lV3b',...
    'rV1','rV2','rV3','rhV4','rLO1','rLO2','rTO1','rTO2','rV3a','rV3b'};
% ROIs = {'lV1','lV2','lV3','lhV4','lLO1','lLO2','lTO1','lTO2','lV3a','lV3b','lV01','lV02',...
%     'rV1','rV2','rV3','rhV4','rLO1','rLO2','rTO1','rTO2','rV3a','rV3b','rV01','rV02'};
% ROIs = {'lV1','lV2','lV3','lhV4','lLO1','lLO2','lTO1','lTO2','lV3a','lV3b','lV01',...
%     'rV1','rV2','rV3','rhV4','rLO1','rLO2','rTO1','rTO2','rV3a','rV3b','rV01'};
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
            's008720190313_3T','s009220190327_3T'};
end

numRois = length(ROIs);
% subNames = {'s005920190129','s006220190107','s007720190204'};
% subNames = {'s005920190129','s006220190107'};
numSubs = length(subNames);
% dataFolder = '/Volumes/ZviDrive/noah/drift_analysis_attn/';
saveFolder = '~/noah/';
cd(dataFolder);
for isub=1:numSubs
    cd(subNames{isub});
    
    v=newView;
    for iroi=1:numRois
        switch directionOrIllusion
            case 1
                [ res ] = plotROI_localizer(v, ROIs{iroi}, crossValidation,locThresh);
                
                for iloc=1:2
                    decodingAcc(isub,iroi,iloc) = res{iloc}.accuracy;
                    randAccDist(isub,iroi,iloc,:) = res{iloc}.randAccuracy;
                    rand95acc(isub,iroi,iloc) = res{iloc}.rand95accuracy;
                    numVox(isub,iroi,iloc) = res{iloc}.numVox;
                end
                sharedVox(isub,iroi) = res{1}.sharedVox;
                totalVox(isub,iroi) = res{1}.totalVox;
            case 2
                [ res ] = plotROI_localizer_illusion(v, ROIs{iroi}, crossValidation,locThresh);
                
                for iloc=1:2
                    decodingAccLR(isub,iroi,iloc,:) = res{iloc}.accuracy;
                    randAccDistLR(isub,iroi,iloc,:,:) = res{iloc}.randAccuracy;
                    rand95accLR(isub,iroi,iloc,:) = res{iloc}.rand95accuracy;
                    
                    decodingAcc(isub,iroi,iloc) = mean(res{iloc}.accuracy);
                    randAccDist(isub,iroi,iloc,:) = mean(res{iloc}.randAccuracy,2);
                    rand95acc(isub,iroi,iloc) = mean(res{iloc}.rand95accuracy);
                    
                    numVox(isub,iroi,iloc) = res{iloc}.numVox;
                end
                sharedVox(isub,iroi) = res{1}.sharedVox;
                totalVox(isub,iroi) = res{1}.totalVox;
                
        end
    end
    deleteView(v);
    cd ..
end
%%
varinfo=whos('randAccDist');
saveopt='';
if varinfo.bytes >= 2^31
    saveopt='-v7.3';
end
if directionOrIllusion==1
    save([saveFolder 'decodeAccuracy_' expName '_' crossValidation '_thresh' locStr '.mat'],'crossValidation','locThresh','ROIs','subNames','dataFolder',...
        'decodingAcc','randAccDist','rand95acc','numVox','sharedVox','totalVox',saveopt);
else
    varinfo=whos('randAccDistLR');
    saveopt='';
    if varinfo.bytes >= 2^31
        saveopt='-v7.3';
    end
    save([saveFolder 'decodeAccuracyIllusion_' expName '_' crossValidation '_thresh' locStr '.mat'],'crossValidation','locThresh','ROIs','subNames','dataFolder',...
        'decodingAcc','randAccDist','rand95acc','numVox','sharedVox','totalVox',...
        'decodingAccLR','randAccDistLR', 'rand95accLR', saveopt);
    
end

toc