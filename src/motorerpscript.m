%% this script gets you the RP and P300 response-locked for each subject (Figure 3c)

clear all
erpwindow = [-1200:100];
erpbaseline = [1:200];

whichtrials = 'all';  %or 'correct';
subid = ''; 
eval(['load ' subid '/' subid '_cleaned.mat']);
trials = organizetrials(artifact,blockorder,correctresponses,responseside,responsetime);
[erp,ntrials] = motorERP(data,responsetime,trials,erpwindow,whichtrials);
erp.as.left = baselinecorrect(erp.as.left,erpbaseline);
erp.as.right = baselinecorrect(erp.as.right,erpbaseline);
erp.eo.left = baselinecorrect(erp.eo.left,erpbaseline);
erp.eo.right = baselinecorrect(erp.eo.right,erpbaseline);
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
erp.as.left=filtfilthd(Hd,erp.as.left);
erp.as.right=filtfilthd(Hd,erp.as.right);
erp.eo.left=filtfilthd(Hd,erp.eo.left);
erp.eo.right=filtfilthd(Hd,erp.eo.right);

%carry along some information 
erp.window = erpwindow;
erp.baseline = erpbaseline;
erp.filter.Hd = Hd;
erp.filter.Fpass = 15;
erp.filter.Fstop = 20;
erp.filter.Apass = 1;
erp.filter.Astop = 10;
erp.filter.match = match;
erp.whichtrials = whichtrials;

%save to same directory 
eval(['save -v7.3 ' subid '/' subid '_motorerp.mat' ]);

