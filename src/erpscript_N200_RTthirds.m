%% this script is for getting the N200 for each subject with the RT tertile split
clear all
erpwindow = [-400:1500];
erpbaseline = [300:400];

whichtrials = 'all';  %or 'correct';
subid = ''; %put subject ID in 
eval(['load ' subid '/' subid '_cleaned.mat']);
eval(['load ' subid '/' subid '_RTthirds.mat']);
%as fast
[erp.as.fastest,ntrials.as.fastest] = motorERP_RTsplit_stimlocked(data,trials.as.fastest,erpwindow,responsetime);
[erp.as.medfast,ntrials.as.medfast] = motorERP_RTsplit_stimlocked(data,trials.as.medfast,erpwindow,responsetime);
[erp.as.slowest,ntrials.as.slowest] = motorERP_RTsplit_stimlocked(data,trials.as.slowest,erpwindow,responsetime);

[erp.eo.fastest,ntrials.eo.fastest] = motorERP_RTsplit_stimlocked(data,trials.eo.fastest,erpwindow,responsetime);
[erp.eo.medfast,ntrials.eo.medfast] = motorERP_RTsplit_stimlocked(data,trials.eo.medfast,erpwindow,responsetime);
[erp.eo.slowest,ntrials.eo.slowest] = motorERP_RTsplit_stimlocked(data,trials.eo.slowest,erpwindow,responsetime);



%baseline correct

erp.as.fastest = baselinecorrect(erp.as.fastest,erpbaseline);
erp.as.medfast = baselinecorrect(erp.as.medfast,erpbaseline);
erp.as.slowest = baselinecorrect(erp.as.slowest,erpbaseline);
erp.eo.fastest = baselinecorrect(erp.eo.fastest,erpbaseline);
erp.eo.medfast = baselinecorrect(erp.eo.medfast,erpbaseline);
erp.eo.slowest = baselinecorrect(erp.eo.slowest,erpbaseline);
%Filter data
% Lowpass the data
Fpass = 10;          % Passband Frequency
Fstop = 20;          % Stopband Frequency
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
erp.window = erpwindow;
erp.baseline = erpbaseline;
erp.filter.Hd = Hd;
erp.filter.Fpass = 10;
erp.filter.Fstop = 20;
erp.filter.Apass = 1;
erp.filter.Astop = 10;
erp.filter.match = match;


%save to same directory 
eval(['save -v7.3 ' subid '/' subid '_N200_RTthirds.mat erp  whichtrials' ]);