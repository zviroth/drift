close all
clear all
runDecoding=1;
pursuitError=1;
dataFolder='/Volumes/Drift/DriftIllusion_Attention/3TData/';
subNames = {'s007320190227_3T','s007620190320_3T','s008620190417_3T','s008720190313_3T','s009220190327_3T'};
% subNames = {'s007320190227_3T','s007620190320_3T','s008720190313_3T','s009220190327_3T'};
% subNames = {'s008720190313_3T','s009220190327_3T'};
% subNames = {'s007320190227_3T'};
% subNames = {'s008620190417_3T'};
pursuitErrorStr='';
if pursuitError
    pursuitErrorStr = '_pursuitError';
end
for isub=1:length(subNames)
    
    if ~strcmp(subNames{isub},'s009220190327_3T')
        filename = fullfile(dataFolder,subNames{isub},'Etc',['driftEye_right' pursuitErrorStr '_ZR.mat']);
%         filename = fullfile(dataFolder,subNames{isub},'Etc',['driftEye_right' pursuitErrorStr '.mat']);
    else
        filename = fullfile(dataFolder,subNames{isub},'Etc',['driftEye_left' pursuitErrorStr '_ZR.mat']);
%         filename = fullfile(dataFolder,subNames{isub},'Etc',['driftEye_left' pursuitErrorStr '.mat']);
    end
    if runDecoding
        cd(dataFolder)
        cd(subNames{isub})
        driftEye3T_ZR(pursuitError, 0);
        cd ..
    end
    load(filename,'eyeDecoding');

    
    subAccIllusion(isub) = eyeDecoding.IllvsNoIll.accuracy_runwise;
    subRandAccIllusion(isub,:)=eyeDecoding.IllvsNoIll.bootAccuracy;
    
    subAccDirection(isub) = eyeDecoding.LvR.accuracy_runwise;
    subRandAccDirection(isub,:)=eyeDecoding.LvR.bootAccuracy;
    
end

meanAccIllusion = mean(subAccIllusion);
meanRandIllusion = mean(subRandAccIllusion);
rand95illusion = prctile(meanRandIllusion,95);

meanAccDirection = mean(subAccDirection);
meanRandDirection = mean(subRandAccDirection);
rand95direction = prctile(meanRandDirection,95);


%%
figure(1)
rows=2;
cols=length(subNames);
for isub=1:length(subNames)
   subplot(rows,cols, isub);
   hist(subRandAccIllusion(isub,:));
   vline(subAccIllusion(isub));
      subplot(rows,cols, cols+isub);
   hist(subRandAccDirection(isub,:));
   vline(subAccDirection(isub));
   
end

figure(2)
rows=1;
cols=2;
subplot(rows,cols,1)
hist(meanRandIllusion)
vline(meanAccIllusion)

subplot(rows,cols,2)
hist(meanRandDirection)
vline(meanAccDirection)

[meanAccIllusion rand95illusion]
[meanAccDirection rand95direction]