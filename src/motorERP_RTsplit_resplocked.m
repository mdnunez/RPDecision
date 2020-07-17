function [erp,ntrials] = motorERP_RTsplit(data,trials,erpwindow,responsetime);
%this function calculates response-locked erp for  data set for when the RT is split into
%tertiles 
%INPUT: data - EEG data
%	responsetime - RT must add 1000 from start of trial
%	trials - the trials to average
%	erpwindow - window of data around response, e.g. [-1000:200];	
%OUTPUT: erp - response-locked ERP. 
if nargin < 4
error('missing inputs')
end;
%create vectors of erp - 
erp = zeros(length(erpwindow),256);


%ACTION SELECTION
%fastest
ntrials = length(trials);
for j = 1:ntrials
thetrial = trials(j);
erp = erp + squeeze(data(erpwindow+responsetime(thetrial)+1000,:,thetrial));	
end;
erp  = erp/ntrials;
