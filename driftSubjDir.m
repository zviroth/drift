function [ SubjDir ] = driftSubjDir( Analysis )
% Given a string, this simple function will return a list of directories to
% do a particular analysis on.
%   Detailed explanation goes here

if strcmp(Analysis, 'LeftvsNI')
    % This subject list has been used for the following analyses:
    % Left vs No Illusion
    % Right vs No Illusion
    % Univariate test (GLM 1 and 3 vs 2 [Left and Right vs NI])
    
    % 4/24/19: 15 good subjects, 6 removed for documented reasons
    
    
    SubjDir{1} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s002920180301/';
    SubjDir{2} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s003320180319/';
    %SubjDir{3} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s004720180628/'; %REMOVED BECAUSE BAD MOTION COMPARED TO OTHER SCAN
    SubjDir{4} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s004720180830/';
    SubjDir{5} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005020180430/';
    %SubjDir{6} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005520180628/'; %REMOVED BECAUSE INFERIOR COVERAGE, DOESN'T ALIGN WELL WITH OTHER AVERAGES
    %SubjDir{7} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005720180906/'; % REMOVED BECAUSE INFERIOR COVERAGE, THIS ONE HAS POOR MOTION (COULDN'T GET MOTION FROM OTHER TO COMPARE THO)
    SubjDir{8} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005920181115/';
    SubjDir{9} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s006320181015/';
    SubjDir{10} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s007320181116/';
    SubjDir{11} = '/Users/steinbergnj/data/DriftIllusion_Attention/s005720190222/';
    %SubjDir{12} = '/Users/steinbergnj/data/DriftIllusion_Attention/s005920190129/'; % REMOVED BECAUSE (SLIGHTLY) STRANGER ACCURACY MAP. tSNR was comperable, couldn't get plotZest, accuracy and attn are high
    SubjDir{13} = '/Users/steinbergnj/data/DriftIllusion_Attention/s006220190107/';
    SubjDir{14} = '/Users/steinbergnj/data/DriftIllusion_Attention/s006320190211/';
    %SubjDir{15} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007420190128/'; % REMOVED FOR BAD COVERAGE
    %SubjDir{16} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007420190207/'; % REMOVED FOR BAD COVERAGE
    SubjDir{17} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007720190204/'; 
    SubjDir{18} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007720190214/'; 
    SubjDir{19} = '/Users/steinbergnj/data/DriftIllusion_Attention/s008020190201/';
    SubjDir{20} = '/Users/steinbergnj/data/DriftIllusion_Attention/s009320190422/';
    SubjDir{21} = '/Users/steinbergnj/data/DriftIllusion_Attention/s008820190422/';
    
elseif strcmp(Analysis, 'Localizers')
    
    % This subject list has been used for the following analyses:
    % Coherence and phase of stim localizer
    % Coherence and phase of eye movement localizer
    
    % 4/24/19: 14 good subjects, 5 removed for documented reasons 
    
    SubjDir{1} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s002920180301/';
    SubjDir{2} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s003320180319/';
    % SubjDir{3} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s004720180628/'; % Data was NOT grabbed because I didn't want to duplicate this subject
    SubjDir{4} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s004720180830/';
    SubjDir{5} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005020180430/';
    SubjDir{6} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005520180628/';
    SubjDir{7} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005720180906/';
    SubjDir{8} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005920181115/';
    % SubjDir{9} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s006320181015/'; % Data was NOT grabbed because only 1 run instead of 2 or 3 (like other subjects) 
    % SubjDir{10} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s007320181116/';% Data was NOT grabbed because I didn't want to duplicate this subject
    SubjDir{12} = '/Users/steinbergnj/data/DriftIllusion_Attention/s005920190129/';
    SubjDir{13} = '/Users/steinbergnj/data/DriftIllusion_Attention/s006220190107/';
    SubjDir{14} = '/Users/steinbergnj/data/DriftIllusion_Attention/s006320190211/';
    %SubjDir{15} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007420190128/'; % Data was NOT grabbed because weirdly shaped brain!
    %SubjDir{16} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007420190207/'; % Data was NOT grabbed because weirdly shaped brain!
    SubjDir{17} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007720190204/';
    SubjDir{19} = '/Users/steinbergnj/data/DriftIllusion_Attention/s008020190201/';
    SubjDir{20} = '/Users/steinbergnj/data/DriftIllusion_Attention/s009320190422/';
    SubjDir{21} = '/Users/steinbergnj/data/DriftIllusion_Attention/s008820190422/';
    
elseif strcmp(Analysis, 'ROIs')
    
    % This subject list has been used for the following analyses:
    % plotROI 
    
    % 4/25/19: 18 good subjects 
    SubjDir{1} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s002920180301/';
    SubjDir{2} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s003320180319/';
    SubjDir{3} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s004720180628/';
    SubjDir{4} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s004720180830/';
    SubjDir{5} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005020180430/';
    SubjDir{6} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005520180628/';
    SubjDir{7} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005720180906/';
    SubjDir{8} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005920181115/';
    SubjDir{9} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s006320181015/';
    SubjDir{10} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s007320181116/';
    SubjDir{11} = '/Users/steinbergnj/data/DriftIllusion_Attention/s005720190222/';
    SubjDir{12} = '/Users/steinbergnj/data/DriftIllusion_Attention/s005920190129/';
    SubjDir{13} = '/Users/steinbergnj/data/DriftIllusion_Attention/s006220190107/';
    SubjDir{14} = '/Users/steinbergnj/data/DriftIllusion_Attention/s006320190211/';
    SubjDir{15} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007720190204/';
    SubjDir{16} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007720190214/';
    SubjDir{17} = '/Users/steinbergnj/data/DriftIllusion_Attention/s008020190201/';
    SubjDir{18} = '/Users/steinbergnj/data/DriftIllusion_Attention/s009320190422/';
    SubjDir{19} = '/Users/steinbergnj/data/DriftIllusion_Attention/s008820190422/';
    
elseif strcmp(Analysis, 'Attn_4_Conds')
    
    % This subject list has been used for the following analyses:
    % Comparison of Left vs Right, with Local Left vs Local Right
    
    % 4/26/19: 5 good subjects
    SubjDir{1} = '/Users/steinbergnj/data/DriftIllusion_Attn_4_Conds/s008820190318/';
    SubjDir{2} = '/Users/steinbergnj/data/DriftIllusion_Attn_4_Conds/s005520190401/';
    SubjDir{3} = '/Users/steinbergnj/data/DriftIllusion_Attn_4_Conds/s008720190405/';
    SubjDir{4} = '/Users/steinbergnj/data/DriftIllusion_Attn_4_Conds/s009320190415/';
    SubjDir{5} = '/Users/steinbergnj/data/DriftIllusion_Attn_4_Conds/s008620190415/';
    
elseif strcmp(Analysis, 'LvR')
    
    % This subject list has been used for the following analyses:
    % Left vs Right
    
    % 4/26/19: 26 good subjects, 8 removed for documented reasons
    
    SubjDir{1} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s002920180301/';
    SubjDir{2} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s003320180319/';
    %SubjDir{3} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s004720180628/'; %REMOVED BECAUSE BAD MOTION COMPARED TO OTHER SCAN
    SubjDir{4} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s004720180830/';
    SubjDir{5} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005020180430/';
    %SubjDir{6} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005520180628/'; %REMOVED BECAUSE INFERIOR COVERAGE, DOESN'T ALIGN WELL WITH OTHER AVERAGES
    %SubjDir{7} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005720180906/'; % REMOVED BECAUSE INFERIOR COVERAGE, THIS ONE HAS POOR MOTION (COULDN'T GET MOTION FROM OTHER TO COMPARE THO)
    SubjDir{8} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s005920181115/';
    SubjDir{9} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s006320181015/';
    SubjDir{10} = '/Users/steinbergnj/data/Drift_Decision_3_Conds/s007320181116/';
    SubjDir{11} = '/Users/steinbergnj/data/DriftIllusion_Attention/s005720190222/';
    %SubjDir{12} = '/Users/steinbergnj/data/DriftIllusion_Attention/s005920190129/'; % REMOVED BECAUSE (SLIGHTLY) STRANGER ACCURACY MAP. tSNR was comperable, couldn't get plotZest, accuracy and attn are high
    SubjDir{13} = '/Users/steinbergnj/data/DriftIllusion_Attention/s006220190107/';
    SubjDir{14} = '/Users/steinbergnj/data/DriftIllusion_Attention/s006320190211/';
    %SubjDir{15} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007420190128/'; % REMOVED FOR BAD COVERAGE, weird brain
    %SubjDir{16} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007420190207/'; % REMOVED FOR BAD COVERAGE, weird brain
    SubjDir{17} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007720190204/'; % Don't have a good sense of which one to eliminate, keeping both right now
    SubjDir{18} = '/Users/steinbergnj/data/DriftIllusion_Attention/s007720190214/'; % Don't have a good sense of which one to eliminate, keeping both right now
    SubjDir{19} = '/Users/steinbergnj/data/DriftIllusion_Attention/s008020190201/';
    %SubjDir{20} = '/Users/steinbergnj/data/DriftIllusion_Attention/s009320190422/'; %removed because other scan below is better (?)
    %SubjDir{21} = '/Users/steinbergnj/data/DriftIllusion_Attention/s008820190422/'; %removed because other scan below is better (?)
    SubjDir{22} = '/Users/steinbergnj/data/DriftIllusion_Attn_4_Conds/s008820190318/';
    SubjDir{23} = '/Users/steinbergnj/data/DriftIllusion_Attn_4_Conds/s005520190401/';
    SubjDir{24} = '/Users/steinbergnj/data/DriftIllusion_Attn_4_Conds/s008720190405/';
    SubjDir{25} = '/Users/steinbergnj/data/DriftIllusion_Attn_4_Conds/s009320190415/';
    SubjDir{26} = '/Users/steinbergnj/data/DriftIllusion_Attn_4_Conds/s008620190415/';
    
elseif strcmp(Analysis, '3T')
%     SubjDir{1} = '/Users/steinbergnj/data/DriftIllusion_Attention/3TData/s007620190320_3T/';
%     SubjDir{2} = '/Users/steinbergnj/data/DriftIllusion_Attention/3TData/s009220190327_3T/';
    SubjDir{3} = '/Users/steinbergnj/data/DriftIllusion_Attention/3TData/s007320190227_3T/';
%     SubjDir{4} = '/Users/steinbergnj/data/DriftIllusion_Attention/3TData/s008720190313_3T/';
%     SubjDir{5} = '/Users/steinbergnj/data/DriftIllusion_Attention/3TData/s008620190417_3T/';
    
else
    keyboard
end


end

