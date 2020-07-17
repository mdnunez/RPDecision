%% this script is for getting the N200 for each subject with the RT tertile split (figure 3a)
clear all
erpwindow = [-400:1500];
erpbaseline = [300:400];

whichtrials = 'all';  %or 'correct';
subid = '';
eval(['load ' subid '/' subid '_cleaned.mat']);
trials = organizetrials(artifact,blockorder,correctresponses,responseside,responsetime);
[erp,ntrials] = motorERP_stimlocked(data,responsetime,trials,erpwindow,whichtrials);
erp.as.left = baselinecorrect(erp.as.left,erpbaseline);
erp.as.right = baselinecorrect(erp.as.right,erpbaseline);
erp.eo.left = baselinecorrect(erp.eo.left,erpbaseline);
erp.eo.right = baselinecorrect(erp.eo.right,erpbaseline);
%Filter data

% Highpass the data
Fstop = 0.25;        % Stopband Frequency
Fpass = 1;           % Passband Frequency
Astop = 10;          % Stopband Attenuation (dB)
Apass = 1;           % Passband Ripple (dB)
match = 'passband';  % Band to match exactly

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.highpass(Fstop, Fpass, Astop, Apass, 1000);
Hd = design(h, 'butter', 'MatchExactly', match);

% Carry out the filtering
erp.as.left =filtfilthd(Hd,erp.as.left);
erp.as.right =filtfilthd(Hd,erp.as.right);
erp.eo.left =filtfilthd(Hd,erp.eo.left);
erp.eo.right =filtfilthd(Hd,erp.eo.right);

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
erp.as.left=filtfilthd(Hd,erp.as.left);
erp.as.right=filtfilthd(Hd,erp.as.right);
erp.eo.left=filtfilthd(Hd,erp.eo.left);
erp.eo.right=filtfilthd(Hd,erp.eo.right);

%carry along some information 
erp.window = erpwindow;
erp.baseline = erpbaseline;
erp.filter.Hd = Hd;
erp.filter.Fpass = Fpass;
erp.filter.Fstop = Fstop;
erp.filter.Apass = 1;
erp.filter.Astop = 10;
erp.filter.match = match;
erp.whichtrials = whichtrials;

%save to same directory 
eval(['save -v7.3 ' subid '/' subid '_erp_N200.mat' ]);

