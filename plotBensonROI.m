mrQuit
close all
clear all
tic
if exist('v','var')
    deleteView(v);
end
% directionOrIllusion = 2;%1 for direction, 2 for illusion
expNum=2;
illusion=0;
ROIs = {'lV1','lV2','lV3','lhV4','lLO1','lLO2','lTO1','lTO2','lV3a','lV3b','lV01','lV02',...
    'rV1','rV2','rV3','rhV4','rLO1','rLO2','rTO1','rTO2','rV3a','rV3b','rV01','rV02'};

switch expNum
    case 1
        expName = '3conds';
        dataFolder = '/Volumes/Drift/Drift_Decision_3_Conds/';
        subNames = {'s002920180301','s004020180427','s004720180830','s005020180430',...
            's005920181115','s006320181015','s003320180319','s004720180628',...
            's005120180524','s005720180906','s006020181001','s007320181116'};
    case 2
        expName = 'attn';
        dataFolder = '/Volumes/Drift/DriftIllusion_Attention/';
        subNames = {'s007420190128','s007420190207','s008820190422','s005920190129',...
            's009320190422','s006220190107','s007720190204','s006320190211','s007920190114',...
            's008020190201'};
            subNames = {'s007720190204'};
%             subNames = {'s007420190128','s007420190207','s008820190422','s005920190129',...
%             's009320190422','s006220190107','s006320190211','s007920190114',...
%             's008020190201'};
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
    if illusion==1
    [stimvol stimNames var] = getStimvol(v,{'groupNum=[1 3]','groupNum=2'},'taskNum=1','phaseNum=1');
    elseif illusion==0%direction
        [stimvol stimNames var] = getStimvol(v,{'groupNum=[1]','groupNum=3'},'taskNum=1','phaseNum=1');
    end
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

            t = 0:1.5:24-1.5;
            Colors = [[0/255 171/255 75/255]; [0/255 126/255 200/255]];
            
            h = figure(1);
            clf
            hold on
            %plot it so that the no illusion (blue) is under the illusion (green)
            markersize = 8;
            linewidth=10;
            i=2;
            myerrorbar(t, d.ehdr(i,:)', 'yError', d.ehdrste(i,:)', 'Color', Colors(i,:), 'MarkerSize', markersize, 'LineWidth', linewidth)
            i=1;
            myerrorbar(t, d.ehdr(i,:)', 'yError', d.ehdrste(i,:)', 'Color', Colors(i,:), 'MarkerSize', markersize, 'LineWidth', linewidth)
            
            h1=plot(0,d.ehdr(1,1), 'Color', Colors(1,:), 'MarkerSize', 1);
            h2=plot(0,d.ehdr(2,1), 'Color', Colors(2,:), 'MarkerSize', 1);
            
           
            
            axis square; axis tight
            YLabel = sprintf('fMRI response \n(%%change image intensity)');
            drawPublishAxis('xLabelOffset', -4/64,'yLabelOffset', -9/64, 'xAxisMargin', 2/64, 'yAxisMargin', 2/64,...
                'xTick', [0 24], 'yTick', [-1.5 0 1.5], 'xLabel', 'Time (sec)', ...
                'yLabel',  YLabel,  ...
                'figSize', [17.6], 'labelFontSize', 25);
            legend off
            if illusion==1
            legend([h1 h2],{'Illusion','No Illusion'});
            elseif illusion==0
                legend([h1 h2],{'Left','Right'});
            end
            leg=legend;
            leg.FontSize = 19;
            leg.Box='off';
            hl = findobj(gca,'type','line');
            set(hl,'LineWidth',2.5);
            
            h.Position=[5   5   22.0839   17.6389];
            if illusion==1
            filename = sprintf('Figures/tsPlot_%s.pdf', roiname);
            elseif illusion==0
                filename = sprintf('Figures/tsPlot_direction_%s.pdf', roiname);
            end
            
            if ~(exist('Figures')==7)
                mkdir('Figures')
            end
            savepdf(h, filename)
            
        end
    end
    deleteView(v);
    cd ..
end


%print('-painters','-dpdf',['/Volumes/Drift/DriftIllusion_Attention/s009320190422/flatPatchROIs.pdf']);
%print('-painters','-dpdf',['~/Documents/MATLAB/min/figures/fig2_circular_cmap.pdf']);

% print('-painters','-dpdf',['flatPatchROIs_black.pdf']);