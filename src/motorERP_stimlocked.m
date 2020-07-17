function [erp,ntrials] = motorERP_stimlocked(data,responsetime,trials,window,whichtrials);
%this function calculates motor erp for stroke2AFC data set.
%there are 4 ERP averages Execution Only -Left, Execution Only-Right, 
%Action Selection-Left, Action-Selection - Right
%INPUT: data - EEG data
%	responsetime - RT must add 1000 from start of trial
%	trials - output of organize trials
%	window - window of data around response, e.g. [-1000:200];
%	whichtrials - 'all' for all trials, 'correct' for correct trials only	
%OUTPUT: erp - structure with motor-locked ERP. 
if nargin < 5
error('missing inputs')
end;
erp.eo.left = zeros(length(window),256);
erp.eo.right = erp.eo.left;
erp.as.left = erp.eo.left;
erp.as.right = erp.as.left;
if strcmp(whichtrials,'all')
%EXECUTION ONLY
ntrials.eo.left = length(trials.eo.left);
for j = 1:ntrials.eo.left;
thetrial = trials.eo.left(j);
erp.eo.left = erp.eo.left + squeeze(data(window+1000,:,thetrial));	
end;
erp.eo.left = erp.eo.left/ntrials.eo.left;
ntrials.eo.right = length(trials.eo.right);
for j = 1:ntrials.eo.right;
thetrial = trials.eo.right(j);
erp.eo.right = erp.eo.right + squeeze(data(window+1000,:,thetrial));	
end;
erp.eo.right = erp.eo.right/ntrials.eo.right;
%ACTION SELECTION
ntrials.as.left = length(trials.as.left);
for j = 1:ntrials.as.left;
thetrial = trials.as.left(j);
erp.as.left = erp.as.left + squeeze(data(window+1000,:,thetrial));	
end;
erp.as.left = erp.as.left/ntrials.as.left;
ntrials.as.right = length(trials.as.right);
for j = 1:ntrials.as.right;
thetrial = trials.as.right(j);
erp.as.right = erp.as.right + squeeze(data(window+1000,:,thetrial));	
end;
erp.as.right = erp.as.right/ntrials.as.right;
elseif strcmp(whichtrials,'correct');
%EXECUTION ONLY
ntrials.eo.left = length(trials.correct.eo.left);
for j = 1:ntrials.eo.left;
thetrial = trials.correct.eo.left(j);
erp.eo.left = erp.eo.left + squeeze(data(window+1000,:,thetrial));	
end;
erp.eo.left = erp.eo.left/ntrials.eo.left;
ntrials.eo.right = length(trials.correct.eo.right);
for j = 1:ntrials.eo.right;
thetrial = trials.correct.eo.right(j);
erp.eo.right = erp.eo.right + squeeze(data(window+1000,:,thetrial));	
end;
erp.eo.right = erp.eo.right/ntrials.eo.right;
%ACTION SELECTION
ntrials.as.left = length(trials.correct.as.left);
for j = 1:ntrials.as.left;
thetrial = trials.correct.as.left(j);
erp.as.left = erp.as.left + squeeze(data(window+1000,:,thetrial));	
end;
erp.as.left = erp.as.left/ntrials.as.left;
ntrials.as.right = length(trials.correct.as.right);
for j = 1:ntrials.as.right;
thetrial = trials.correct.as.right(j);
erp.as.right = erp.as.right + squeeze(data(window+1000,:,thetrial));	
end;
erp.as.right = erp.as.right/ntrials.as.right;
else
exit('Unknown trial types')
end;


 

