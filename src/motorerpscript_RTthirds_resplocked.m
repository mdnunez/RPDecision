%% getting the P300 and RP based on the tertile split (Figure 5a)
clear all
erpwindowslow = [-1400:100];
erpwindowmed = [-1200:100];
erpwindowfast = [-1200:100];
erpbaseline = [1:200];

whichtrials = 'all';  %or 'correct';
subid = '';
eval(['load ' subid '/' subid '_cleaned.mat']);
eval(['load ' subid '/' subid '_RTthirds.mat']);
%trials = organizetrials(artifact,blockorder,correctresponses,responseside,responsetime);
%have to do this separately for each trial and condition
%as fast
[erp.as.fastest,ntrials.as.fastest] = motorERP_RTsplit_resplocked(data,trials.as.fastest,erpwindowfast,responsetime);
[erp.as.medfast,ntrials.as.medfast] = motorERP_RTsplit_resplocked(data,trials.as.medfast,erpwindowmed,responsetime);
[erp.as.slowest,ntrials.as.slowest] = motorERP_RTsplit_resplocked(data,trials.as.slowest,erpwindowslow,responsetime);

[erp.eo.fastest,ntrials.eo.fastest] = motorERP_RTsplit_resplocked(data,trials.eo.fastest,erpwindowfast,responsetime);
[erp.eo.medfast,ntrials.eo.medfast] = motorERP_RTsplit_resplocked(data,trials.eo.medfast,erpwindowmed,responsetime);
[erp.eo.slowest,ntrials.eo.slowest] = motorERP_RTsplit_resplocked(data,trials.eo.slowest,erpwindowslow,responsetime);



%baseline correct

erp.as.fastest = baselinecorrect(erp.as.fastest,erpbaseline);
erp.as.medfast = baselinecorrect(erp.as.medfast,erpbaseline);
erp.as.slowest = baselinecorrect(erp.as.slowest,erpbaseline);
erp.eo.fastest = baselinecorrect(erp.eo.fastest,erpbaseline);
erp.eo.medfast = baselinecorrect(erp.eo.medfast,erpbaseline);
erp.eo.slowest = baselinecorrect(erp.eo.slowest,erpbaseline);
%Filter data
% Lowpass the data
Fpass = 4;          % Passband Frequency
Fstop = 8;          % Stopband Frequency
Apass = 1;           % Passband Ripple (dB)
Astop = 10;          % Stopband Attenuation (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, 1000);
Hd = design(h, 'butter', 'MatchExactly', match);

% Carry out the filtering
erp.as.fastest=filtfilthd(Hd,erp.as.fastest);
erp.as.medfast=filtfilthd(Hd,erp.as.medfast);
erp.as.slowest = filtfilthd(Hd,erp.as.slowest);
erp.eo.fastest=filtfilthd(Hd,erp.eo.fastest);
erp.eo.medfast=filtfilthd(Hd,erp.eo.medfast);
erp.eo.slowest = filtfilthd(Hd,erp.eo.slowest);
%carry along some information 
erp.window.slow = erpwindowslow;
erp.window.fast = erpwindowfast;
erp.window.med = erpwindowmed;
erp.baseline = erpbaseline;
erp.filter.Hd = Hd;
erp.filter.Fpass = 15;
erp.filter.Fstop = 20;
erp.filter.Apass = 1;
erp.filter.Astop = 10;
erp.filter.match = match;


%save to same directory 
eval(['save -v7.3 ' subid '/' subid '_motorerp_RTthirds_new_resplocked.mat' ]);