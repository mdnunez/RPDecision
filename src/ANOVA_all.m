%% this script computes the ANOVA for erp onset/peaks/behavioral results
%% ANOVA FOR N200
load('n200peak_svd.mat');

%organize onset times data for anova - resp onset
data = [n200peak.as.all'; n200peak.eo.all'];
task = [ones(14,1); 2*ones(14,1)];

subj = [[1:14]';[1:14]'];

[P,T,STATS,TERMS]=anovan(data,{task;subj;},'random',[2],'model','interaction','varnames',{'Task';'Subject';});
%% ANOVA for P300 WITHIN tasks, on P300 onset times baed on the tertile split
load('n200peak_svd.mat');

%organize onset times data for anova - resp onset
data = [n200peak.as.all'; n200peak.eo.all'];
task = [ones(14,1); 2*ones(14,1)];

subj = [[1:14]';[1:14]'];

[P,T,STATS,TERMS]=anovan(data,{task;subj;},'random',[2],'model','interaction','varnames',{'Task';'Subject';});

%% ANOVA for P300 onset for stimulus-locked and response-locked onset
load('p300peak_svd.mat');
%organize onset times data for anova - rtthirds - stim onset
data = [p300peak.as.fastest'; p300peak.as.medfast'; p300peak.as.slowest'; p300peak.eo.fastest'; 
    p300peak.eo.medfast'; p300peak.eo.slowest'];
task = [ones(42,1); 2*ones(42,1)];
thirds = [ones(14,1); 2*ones(14,1); 3*ones(14,1); ones(14,1); 2*ones(14,1); 3*ones(14,1)];
subj = [[1:14]';[1:14]';[1:14]';[1:14]';[1:14]';[1:14]'];

[P,T,STATS,TERMS]=anovan(data,{task;thirds;subj},'random',[3],'model','interaction','varnames',{'Task';'Thirds';'Subject'});

%organize onset times data for anova - resp onset
data = [p300peak.as.all'; p300peak.eo.all'];
task = [ones(14,1); 2*ones(14,1)];
subj = [[1:14]';[1:14]'];

[P,T,STATS,TERMS]=anovan(data,{task;subj;},'random',[2],'model','interaction','varnames',{'Task';'Subject';});

%% ANOVA for P300 peak latency (stimulus-locked)
load p300slope.mat

data = [p300peakvolt.as'; p300peakvolt.eo'];
% data = [rpvolt.as.left'; rpvolt.as.right'; rpvolt.eo.left'; rpvolt.eo.right'];
task = [ones(14,1); 2*ones(14,1)];
%side = [ones(14,1); 2*ones(14,1); ones(14,1); 2*ones(14,1)];
subj = [[1:14]';[1:14]'];

[P,T,STATS,TERMS]=anovan(data,{task;subj;},'random',[2],'model','interaction','varnames',{'Task';'Subject';});

%% ANOVA for RP onset times - stimulus-locked and response-locked and when split into tertiles
load rp_onsets.mat

%organize onset times data for anova - resp onset
data = [rponset.as.all'; rponset.eo.all'];
% data = [rponset.as.left'; rponset.as.right'; rponset.eo.left'; rponset.eo.right'];
task = [ones(14,1); 2*ones(14,1)];
%side = [ones(14,1); 2*ones(14,1); ones(14,1); 2*ones(14,1)];
subj = [[1:14]';[1:14]'];

[P,T,STATS,TERMS]=anovan(data,{task;subj;},'random',[2],'model','interaction','varnames',{'Task';'Subject';});

%organize onset times data for anova - stim onsets
data = [stimrponsets.as.all'; stimrponsets.eo.all'];
task = [ones(14,1); 2*ones(14,1)];
%side = [ones(14,1); 2*ones(14,1); ones(14,1); 2*ones(14,1)];
subj = [[1:14]';[1:14]'];

[P,T,STATS,TERMS]=anovan(data,{task;subj},'random',[2],'model','interaction','varnames',{'Task';'Subject'});


%organize onset times data for anova - rtthirds - resp onset
data = [rponset.as.fastest'; rponset.as.medfast'; rponset.as.slowest'; rponset.eo.fastest'; 
    rponset.eo.medfast'; rponset.eo.slowest'];
task = [ones(42,1); 2*ones(42,1)];
thirds = [ones(14,1); 2*ones(14,1); 3*ones(14,1); ones(14,1); 2*ones(14,1); 3*ones(14,1)];
subj = [[1:14]';[1:14]';[1:14]';[1:14]';[1:14]';[1:14]'];

[P,T,STATS,TERMS]=anovan(data,{task;thirds;subj},'random',[3],'model','interaction','varnames',{'Task';'Thirds';'Subject'});


%organize onset times data for anova - rtthirds - stim onset
data = [stimrponsets.as.fastest'; stimrponsets.as.medfast'; stimrponsets.as.slowest'; stimrponsets.eo.fastest'; 
    stimrponsets.eo.medfast'; stimrponsets.eo.slowest'];
task = [ones(42,1); 2*ones(42,1)];
thirds = [ones(14,1); 2*ones(14,1); 3*ones(14,1); ones(14,1); 2*ones(14,1); 3*ones(14,1)];
subj = [[1:14]';[1:14]';[1:14]';[1:14]';[1:14]';[1:14]'];

[P,T,STATS,TERMS]=anovan(data,{task;thirds;subj},'random',[3],'model','interaction','varnames',{'Task';'Thirds';'Subject'});

%% ANOVA for RTs
load young_rt.mat

%organize onset times data for anova - resp onset
%data = [allavgrt.as.all'; allavgrt.eo.all'];
data = [allavgrt.as.left'; allavgrt.as.right'; allavgrt.eo.left'; allavgrt.eo.right'];
task = [ones(28,1); 2*ones(28,1)];
side = [ones(14,1); 2*ones(14,1); ones(14,1); 2*ones(14,1)];
subj = [[1:14]';[1:14]';[1:14]';[1:14]'];

[P,T,STATS,TERMS]=anovan(data,{task;side;subj;},'random',[2 3],'model','interaction','varnames',{'Task';'Rotation';'Subject';});
